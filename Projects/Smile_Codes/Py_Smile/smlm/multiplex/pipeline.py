# smlm/multiplex/pipeline.py
"""一站式配准流程: 依次载入多个文件 -> 各自配对 -> 合并拟合 -> 报告。

典型用法::

    from smlm.multiplex import calibrate_dual_view, DualViewROI
    calib, report = calibrate_dual_view(
        [r'...\\std2.sif', r'...\\std.sif'],
        roi=DualViewROI(split_x=256, row_range=(0, 512)),
        order='auto', cap_per_file=2000)
    calib.save(r'...\\dualcal_v2.npz')

For general multichannel inputs, ``pair_multichannel_files(files, roi)`` returns
one merged ``PairSet`` per moving channel, all referenced to the selected
channel (channel 0 by default). ``fit_multichannel_pairs(pairs, roi)`` fits
those pairs and returns channel-aligned calibrations with the reference slot
set to ``None``.

要点:
  - 每个文件**独立配对**(不能跨文件配), 拟合时再合并 —— 因为不同文件的
    粗偏移/漂移状态不一样, 混在一起配会互相干扰。
  - 合并时按文件配额抽样, 避免亮文件主导。
  - ``localize_channel`` 可显式读取或覆盖包含过滤后定位的缓存。
"""
from __future__ import annotations

import os
import warnings
import numpy as np

from ..datatypes import ImageStack, MultiImageStack
from ..camera import EMCCD
from ..localize import (ClassicalPipeline, LocalMaxDetector, MLELocalizer,
                        localize_stack)
from ..io import load_raw

from .transform import DualViewCalibration, NPZ_VERSION
from .pair import PairSet, pair_channels, pair_multichannels, merge_pairs
from .channels import (DualViewROI, split_dual_view, cut_pair,
                       refine_coarse_offset, check_drift_margin,
                       seam_windows, split_from_crops, windows_from_crops,
                       DEFAULT_SEAM)
from .fit import fit_transform, residual_report
from .roi import MultichannelROI
from .bundle import MultiChannelCalibrationBundle

# ---------------------------------------------------------------------------
#  标定用的默认参数 (分子密度低、要求干净多过要求灵敏)
# ---------------------------------------------------------------------------
CALIB_CAM = dict(conv_factor=5.75, em_gain=300)
CALIB_DET_KW = dict(edge=5, f_em=np.sqrt(2), bg_window=15,
                    max_size=7, min_neighbors=2)
CALIB_LOC_KW = dict(edge=5, F_EM=np.sqrt(2), use_pixel=False,
                    use_gradient=True, n_jobs=1, verbose=False)
CALIB_FILT_KW = dict(photons=(15, 1e8), bg=(0, 1e6), sigma=(0.6, 5.0), lp=None)
CALIB_K_SIGMA = dict(green=1.8, red=1.3)
CALIB_PAIR_TOLS = (12.0, 6.0, 3.0, 1.5)


# ---------------------------------------------------------------------------
#  单通道定位 (可缓存)
# ---------------------------------------------------------------------------
def mean_frame(stack, frames=None):
    """把 (T,H,W) 压成单帧 ImageStack(时间平均, 标准样品固定不动时 SNR 更高)。

    ``frames=(a, b)`` 可选, 只平均 [a, b) 这段。
    """
    n = int(stack.shape[0])
    if frames is not None:
        a, b = max(0, int(frames[0])), min(n, int(frames[1]))
        if b <= a:
            raise ValueError(f'frames={frames} 超出帧数 {n}')
        sub, tag = stack.data[a:b], f'{a}-{b}'
    else:
        sub, tag = stack.data, 'mean'
    data = sub.mean(axis=0, keepdims=True)
    meta = {**dict(stack.metadata), 'frame_mode': 'mean', 'frame_window': tag}
    return ImageStack(data=data, pixel_size=stack.pixel_size, metadata=meta), tag


def _mean_stack(stack, frame_mode='mean', frame_window=None):
    """按 frame_mode 决定是否只用某段帧, 然后压成单帧。"""
    if frame_mode == 'window' and frame_window is not None:
        return mean_frame(stack, frame_window)
    return mean_frame(stack)



def localize_channel(stack, *, k_sigma=1.5, det_kw=None, loc_kw=None,
                     filt_kw=None, src_path=None, read_cache=False,
                     verbose=False, tag=''):
    """Compatibility wrapper using calibration defaults for one stack."""
    det_kw = {**CALIB_DET_KW, 'verbose': verbose, **(det_kw or {})}
    loc_kw = {**CALIB_LOC_KW, **(loc_kw or {})}
    filt_kw = {**CALIB_FILT_KW, **(filt_kw or {})}
    return localize_stack(
        stack,
        k_sigma=k_sigma,
        det_kw=det_kw,
        loc_kw=loc_kw,
        filt_kw=filt_kw,
        src_path=src_path,
        read_cache=read_cache,
        verbose=verbose,
        tag=tag,
    )


# ---------------------------------------------------------------------------
#  单文件: 载入 -> 切通道 -> 定位 -> 配对
# ---------------------------------------------------------------------------
def calib_one_file(path, *, roi=None, split_x=None, row_range=None,
                   cam=None, cam_kw=None, frame_mode='mean',
                   frame_window=None, pixel_size=None,
                   k_sigma=None, det_kw=None, loc_kw=None, filt_kw=None,
                   tols=CALIB_PAIR_TOLS, seam=DEFAULT_SEAM, crops=None,
                   read_cache=False, verbose=True):
    """处理一个文件, 返回 (PairSet, DualViewROI, 实际窗口 dict)。

    Parameters
    ----------
    crops : dict, optional
        **人工指定**的两个框, 形如
        ``{'g': dict(y=(y0,y1), x=(x0,x1)), 'r': {...}}``。
        给了就直接用它 —— 取框是人的判断, 见 ``viz.plot_cut_windows``。
        与 ``coords.geometry_from_crops`` 用同一种格式, 两边可以照抄。
    seam : int
        仅在没给 ``crops`` 时生效: 两通道各让出多少 px。默认 0,
        因为"该让多少"取决于光路和焦面, 不该由模块拍一个常数。
    read_cache : bool
        True 时严格读取并验证已有定位缓存; False 时重新定位并覆盖缓存。
    """
    cam_kw = {**CALIB_CAM, **(cam_kw or {})}
    cam = cam or EMCCD(**cam_kw)
    name = os.path.basename(path)

    raw = load_raw(path)
    img = cam.process(raw, pixel_size=pixel_size) if pixel_size is not None \
        else cam.process(raw)
    if not isinstance(img, ImageStack):
        img = ImageStack(data=np.asarray(img.data), pixel_size=getattr(
            img, 'pixel_size', None), metadata=dict(getattr(img, 'metadata', {})))
    if verbose:
        print(f'[{name}] stack {tuple(img.shape)}')

    if split_x is None:
        split_x = roi.split_x if roi is not None else img.shape[2] // 2
    if row_range is None:
        row_range = roi.row_range if roi is not None else (0, img.shape[1])
    roi = roi or DualViewROI(split_x=int(split_x),
                             row_range=tuple(int(v) for v in row_range),
                             n_cols=int(img.shape[2]))

    # 取框: 人工给的 crops 优先; 否则按 seam 从半幅裁
    if crops is not None:
        stack_g, stack_r, windows = split_from_crops(img, crops, verbose=verbose)
    else:
        stack_g, stack_r, roi, _offsets = split_dual_view(
            img, roi, seam=seam, verbose=verbose)
        gx, rx = seam_windows(roi, seam)
        ry = roi.row_range
        hh = int(ry[1] - ry[0])
        windows = {'g': (int(ry[0]), int(gx[0]), hh, int(gx[1] - gx[0])),
                   'r': (int(ry[0]), int(rx[0]), hh, int(rx[1] - rx[0]))}

    # 标定用平均图定位: 光子数叠加后 SNR 高得多, 密度也足够配对
    stack_g, tag_g = _mean_stack(stack_g, frame_mode, frame_window)
    stack_r, tag_r = _mean_stack(stack_r, frame_mode, frame_window)

    ks = {'green': 1.8, 'red': 1.3, **(k_sigma or {})}
    loc_g = localize_channel(stack_g, k_sigma=ks['green'], det_kw=det_kw,
                             loc_kw=loc_kw, filt_kw=filt_kw,
                             src_path=f'{os.path.splitext(path)[0]}_green_{tag_g}',
                             read_cache=read_cache,
                             verbose=verbose, tag=f'{name} green')
    loc_r = localize_channel(stack_r, k_sigma=ks['red'], det_kw=det_kw,
                             loc_kw=loc_kw, filt_kw=filt_kw,
                             src_path=f'{os.path.splitext(path)[0]}_red_{tag_r}',
                             read_cache=read_cache,
                             verbose=verbose, tag=f'{name} red')

    pairs, shift = pair_channels(loc_r, loc_g, tols=tols, verbose=verbose,
                                 tag=f'{name} r->g')
    pairs.file = np.full(len(pairs), name, dtype=object)
    return pairs, roi, windows


def load_multi(path, *, roi=None, split_x=None, row_range=None, cam_kw=None,
               frame_mode='mean', frame_window=None, pixel_size=None,
               coarse=None, margin=0, verbose=True):
    """载入一个文件并切成**同一物理窗口**的两通道 -> MultiImageStack。

    这是下游分析(漂移/定位/轨迹)该用的入口: 两通道已经对齐在同一个物理视场,
    只剩几个 px 的亚像素残差需要坐标层的标定来修。
    """
    from .channels import build_multi
    cam_kw = {**CALIB_CAM, **(cam_kw or {})}
    cam = EMCCD(**cam_kw)
    raw = load_raw(path)
    img = cam.process(raw, pixel_size) if pixel_size is not None else cam.process(raw)
    if not isinstance(img, ImageStack):
        img = ImageStack(data=np.asarray(img.data),
                         pixel_size=getattr(img, 'pixel_size', None),
                         metadata=dict(getattr(img, 'metadata', {})))
    if split_x is None:
        split_x = roi.split_x if roi is not None else img.shape[2] // 2
    if row_range is None:
        row_range = roi.row_range if roi is not None else (0, img.shape[1])
    roi = roi or DualViewROI.from_stack(img, split_x, row_range)
    if coarse is None:
        coarse = (0, roi.split_x)
    return build_multi(img, roi, coarse=coarse, margin=margin,
                       pixel_size=img.pixel_size, verbose=verbose)


# ---------------------------------------------------------------------------
#  多文件配准
# ---------------------------------------------------------------------------
def pair_set_summary(ps):
    """单个文件配对集的摘要: 对数 / 位移中位数 / 去平移残差中位数。

    ``shift`` 是 mov -> ref 的常数平移(取逐点位移的中位数)。多文件合并前要先
    比较这个量: 各文件相差很小才说明它们坐标系相同、光路没动, 合并才有意义。
    """
    if ps is None or len(ps) == 0:
        return {'n_pairs': 0, 'shift_xy': None, 'disp_median': None,
                'resid_median': None}
    d = np.asarray(ps.ref, float) - np.asarray(ps.mov, float)
    sh = np.median(d, axis=0)
    r = d - sh
    return {'n_pairs': int(len(ps)),
            'shift_xy': [float(np.round(sh[0], 3)), float(np.round(sh[1], 3))],
            'disp_median': float(np.round(np.median(np.hypot(*d.T)), 3)),
            'resid_median': float(np.round(np.median(np.hypot(*r.T)), 3))}


def pair_multichannel_files(files, roi, *, reference=None, verbose=True):
    """Read, crop, localize, pair per file, then merge by moving channel.

    ``reference`` defaults to the ROI's first channel. The returned list is
    aligned with ``roi.channels`` and contains ``None`` at the reference slot.
    """
    if isinstance(files, (str, os.PathLike)):
        files = [files]
    else:
        files = list(files)
    if not files:
        raise ValueError('files is empty')
    if not isinstance(roi, MultichannelROI):
        raise TypeError('roi must be a MultichannelROI')
    if len(roi.channels) < 2:
        raise ValueError('roi must define at least two channels')

    channels = list(roi.channels)
    reference = channels[0] if reference is None else reference
    if reference not in channels:
        raise ValueError(f'reference channel {reference!r} is not in roi.channels')
    reference_index = channels.index(reference)
    pair_tols = (6.0,)
    min_pairs = 10
    cap_per_file = 2000
    max_shift_spread = 1.0
    det_kw = {**CALIB_DET_KW, 'verbose': False}
    loc_kw = dict(CALIB_LOC_KW)
    sigma_by_name = {
        'g': CALIB_K_SIGMA['green'],
        'green': CALIB_K_SIGMA['green'],
        'r': CALIB_K_SIGMA['red'],
        'red': CALIB_K_SIGMA['red'],
    }
    detectors = [
        LocalMaxDetector(
            k_sigma=sigma_by_name.get(channel.lower(), 1.5),
            **det_kw,
        )
        for channel in channels
    ]
    localizer = MLELocalizer(**loc_kw)
    camera = EMCCD(**CALIB_CAM)
    file_pairs = [[] for _ in channels]

    for file_index, path in enumerate(files, start=1):
        name = os.path.basename(os.fspath(path))
        if verbose:
            print(f'--- [{file_index}/{len(files)}] {name}')

        raw = load_raw(path)
        image = camera.process(raw)
        if not isinstance(image, ImageStack):
            image = ImageStack(
                data=np.asarray(image.data),
                pixel_size=getattr(image, 'pixel_size', None),
                metadata=dict(getattr(image, 'metadata', {})),
            )
        multi = roi(image)
        mean_stacks = [
            mean_frame(multi.channel(channel))[0]
            for channel in channels
        ]
        channel_locs = [
            localizer.localize(stack, detector.detect(stack))
            for stack, detector in zip(mean_stacks, detectors)
        ]
        if verbose:
            counts = ', '.join(
                f'{channel}={len(locs)}'
                for channel, locs in zip(channels, channel_locs))
            print(f'  localizations: {counts}')

        pairs_for_file = pair_multichannels(
            channel_locs,
            channels=channels,
            reference=reference,
            tols=pair_tols,
            min_pairs=min_pairs,
            verbose=verbose,
        )
        for channel_index, pairs in enumerate(pairs_for_file):
            if channel_index == reference_index:
                continue
            pairs.file[:] = name
            if len(pairs):
                file_pairs[channel_index].append(pairs)
            if verbose:
                print(f'  [{channels[channel_index]} -> {reference}] '
                      f'{len(pairs)} pairs')

    if not any(file_pairs[index] for index in range(len(channels))
               if index != reference_index):
        raise RuntimeError('No multichannel pairs found in any input file')

    merged_pairs = [None] * len(channels)
    for channel_index, channel in enumerate(channels):
        if channel_index == reference_index:
            continue
        sets = file_pairs[channel_index]
        shifts = np.asarray([pairs.shift for pairs in sets], dtype=float)
        if len(shifts) > 1:
            spread = np.ptp(shifts, axis=0)
            if np.any(spread >= max_shift_spread):
                warnings.warn(
                    f'{channel} -> {reference}: per-file shift spread '
                    f'({spread[0]:.2f}, {spread[1]:.2f}) px exceeds '
                    f'{max_shift_spread:.1f} px; merging may combine '
                    'inconsistent calibration sets',
                    RuntimeWarning,
                    stacklevel=2,
                )
        merged_pairs[channel_index] = merge_pairs(
            sets, cap_per_file=cap_per_file)
        if verbose:
            print(f'--- merged {channel} -> {reference}: '
                  f'{len(merged_pairs[channel_index])} pairs '
                  f'from {len(sets)} files')

    return merged_pairs


def fit_multichannel_pairs(multi_channel_pairs, roi, *, order='auto',
                           method='poly', verbose=True):
    """Fit channel-to-reference calibrations from aligned multichannel pairs.

    The ``None`` entry in ``multi_channel_pairs`` identifies the reference
    channel. The returned calibration list preserves that channel ordering.
    """
    if not isinstance(roi, MultichannelROI):
        raise TypeError('roi must be a MultichannelROI')
    if len(roi.channels) < 2:
        raise ValueError('roi must define at least two channels')

    pairs_by_channel = list(multi_channel_pairs)
    if len(pairs_by_channel) != len(roi.channels):
        raise ValueError('pairs length must match roi.channels')
    reference_slots = [
        index for index, pairs in enumerate(pairs_by_channel)
        if pairs is None
    ]
    if len(reference_slots) != 1:
        raise ValueError(
            f'Expected one reference slot, got {len(reference_slots)}')

    reference_index = reference_slots[0]
    reference = roi.channels[reference_index]
    reference_y, reference_x, _, _ = roi.windows[reference]
    calibrations = [None] * len(roi.channels)

    for channel_index, channel in enumerate(roi.channels):
        if channel_index == reference_index:
            continue
        pairs = pairs_by_channel[channel_index]
        if not isinstance(pairs, PairSet):
            raise TypeError(
                f'pairs for channel {channel!r} must be a PairSet')
        if len(pairs) == 0:
            raise ValueError(f'No pairs available for channel {channel!r}')

        moving_y, moving_x, _, _ = roi.windows[channel]
        calibration_roi = {
            'type': 'MultichannelROI',
            'reference': reference,
            'moving': channel,
            'channels': list(roi.channels),
            'y': roi.y,
            'x': roi.x,
            'offsets': {
                name: offset
                for name, offset in zip(roi.channels, roi.offsets)
            },
            'windows': roi.windows,
            'files': pairs.files(),
        }
        if verbose:
            print(f'--- fitting {channel} -> {reference}')
        calibrations[channel_index] = fit_transform(
            pairs,
            order=order,
            method=method,
            roi=calibration_roi,
            coarse_offset=(moving_y - reference_y,
                           moving_x - reference_x),
            reference=reference,
            moving=channel,
            verbose=verbose,
        )

    return calibrations


def multichannel_stack_from_calibration(stack, calibration_path, *,
                                        show_roi=False):
    """Crop an ``ImageStack`` with a saved bundle and attach its registrations.

    Set ``show_roi=True`` to display the bundle's ROI on the input image before
    cropping. Registration is stored in ``MultiImageStack.reg``; image pixels
    are not warped.
    """
    if not isinstance(stack, ImageStack):
        raise TypeError('stack must be an ImageStack')

    bundle = MultiChannelCalibrationBundle.load(calibration_path)
    roi_metadata = bundle.roi
    if tuple(roi_metadata.get('channels', ())) != bundle.channels:
        raise ValueError('calibration ROI channel order does not match bundle')
    if roi_metadata.get('reference') != bundle.reference:
        raise ValueError('calibration ROI reference does not match bundle')

    roi = MultichannelROI(
        y=roi_metadata['y'],
        x=roi_metadata['x'],
        offsets=[roi_metadata['offsets'][channel]
                 for channel in bundle.channels],
        channels=bundle.channels,
    )
    if show_roi:
        roi.show_roi(stack)
        import matplotlib.pyplot as plt
        plt.show()

    cropped = roi(stack)
    metadata = dict(cropped.metadata)
    metadata['multichannel_roi'] = {
        **metadata.get('multichannel_roi', {}),
        'reference': bundle.reference,
    }
    metadata.update({
        'roi': bundle.roi,
        'calibration_bundle': os.fspath(calibration_path),
        'calibration_source_files': bundle.source_files,
    })
    return MultiImageStack(
        data=cropped.data,
        channels=list(bundle.channels),
        reference=bundle.channels.index(bundle.reference),
        reg={
            channel: calibration
            for channel, calibration in zip(
                bundle.channels, bundle.calibrations)
            if calibration is not None
        },
        offsets=cropped.offsets,
        pixel_size=cropped.pixel_size,
        metadata=metadata,
    )


def calibrate_dual_view(files, *, roi=None, split_x=None, row_range=None,
                        reference='g', moving='r', order='auto',
                        cap_per_file=2000, method='poly',
                        robust_k=4.0, robust_iter=4,
                        cam_kw=None, k_sigma=None, det_kw=None, loc_kw=None,
                        filt_kw=None, tols=CALIB_PAIR_TOLS,
                        pixel_size=None, frame_mode='mean',
                        frame_window=None, save_path=None,
                        seam=DEFAULT_SEAM, crops=None,
                        min_consistent_shift=1.0, read_cache=False,
                        verbose=True):
    """依次处理多个标准样品文件, 合并配对后拟合一个统一变换。

    Parameters
    ----------
    files : list[str]
        标准样品(多色荧光小球)的 sif 文件。green 通道作参考。
    cap_per_file : int or None
        每个文件保留的配对数上限, 防止亮文件主导。
    order : int or 'auto'
        多项式阶数。'auto' 用留一交叉验证选阶(见 ``fit.cv_select_order``)。
    seam : int
        两通道在分界线那侧各让出多少 px, 避开双视角的重叠/鬼影带。
        **分析阶段要留同样多**(见 ``channels.cut_pair`` 的 margin /
        ``coords.geometry_from_crops``), 否则分析窗口会比标定覆盖范围更靠外。
    save_path : str, optional
        给了就把标定存成 npz。
    read_cache : bool
        True 时读取已存在的定位缓存; 缺失或不兼容时直接报错。

    Returns
    -------
    (DualViewCalibration, report dict)
    """
    files = [files] if isinstance(files, str) else list(files)
    if not files:
        raise ValueError('files 为空')

    sets, rois, all_windows, file_stats = [], [], [], {}
    for i, f in enumerate(files):
        if verbose:
            print(f'--- [{i + 1}/{len(files)}] {os.path.basename(f)}')
        try:
            ps, roi_i, win = calib_one_file(
                f, roi=roi, split_x=split_x, row_range=row_range,
                cam_kw=cam_kw, k_sigma=k_sigma, det_kw=det_kw, loc_kw=loc_kw,
                filt_kw=filt_kw, tols=tols, pixel_size=pixel_size,
                frame_mode=frame_mode, frame_window=frame_window,
                seam=seam, crops=crops, read_cache=read_cache,
                verbose=verbose)
        except Exception as exc:                     # 单文件失败不该毁掉整批
            if read_cache:
                raise
            print(f'  !! {os.path.basename(f)} 失败: {exc}')
            continue
        name = os.path.basename(f)
        file_stats[name] = pair_set_summary(ps)
        if verbose:
            s = file_stats[name]
            sh = s['shift_xy'] or [float('nan')] * 2
            print(f'  -> {s["n_pairs"]} pairs   '
                  f'shift ({sh[0]:+.2f}, {sh[1]:+.2f}) px')
        if len(ps):
            sets.append(ps)
            rois.append(roi_i)
            all_windows.append(win)

    if not sets:
        raise RuntimeError('所有文件都没配出点, 检查 crops / split_x / 阈值')

    # ---------- 合并前的一致性检查 ----------
    # 合并的合法性前提: 各文件坐标系相同(窗口布局一致) 且光路没变。若逐文件
    # median shift 相差明显, 说明不是同一批次 —— 这时合并, 多项式会把"批次差"
    # 当成形变拟合进去, 得到一个对谁都不准的中间态。
    _sh = np.array([s['shift_xy'] for s in file_stats.values()
                    if s['shift_xy'] is not None])
    if len(_sh) > 1:
        spread = np.ptp(_sh, axis=0)
        consistent = bool(np.all(spread < float(min_consistent_shift)))
    else:
        spread = np.zeros(2)
        consistent = True
    if not consistent:
        warnings.warn(
            f'各文件 median shift 极差 ({spread[0]:.2f}, {spread[1]:.2f}) px '
            f'超过 {min_consistent_shift} px —— 这些文件可能不是同一批次/光路被动过, '
            f'合并会把批次差当成形变拟合进去', RuntimeWarning, stacklevel=2)

    merged = merge_pairs(sets, cap_per_file=cap_per_file)
    if verbose:
        print(f'--- 合并: {len(merged)} pairs 来自 {len(sets)} 个文件'
              f'   shift 极差 ({spread[0]:.3f}, {spread[1]:.3f}) px'
              f'{"  一致" if consistent else "  [!] 不一致"}')

    roi_final = rois[-1] if rois else (roi or DualViewROI(split_x or 256))

    # 探测器中 mov 窗口相对 ref 窗口的位置 = 两者左上角之差 -> (dy, dx)
    if all_windows:
        w = all_windows[0]
        dy = int(w[moving][0] - w[reference][0])
        dx = int(w[moving][1] - w[reference][1])
        gy0, gx0, gh, gw = w[reference]
        ry0, rx0, rh, rw = w[moving]
    else:
        (gx0, gx1), (rx0, rx1) = seam_windows(roi_final, seam)
        gy0, ry0 = int(roi_final.row_range[0]), int(roi_final.row_range[0])
        gh = rh = int(roi_final.row_range[1] - roi_final.row_range[0])
        gw, rw = gx1 - gx0, rx1 - rx0
        dy, dx = ry0 - gy0, rx0 - gx0

    # 存进 roi 的 green_x/red_x 必须是**实际用于标定的**窗口, 因为
    # coords.DualViewGeometry 靠它还原标定窗口、做偏移记账。
    calib = fit_transform(
        merged, order=order, method=method, robust_k=robust_k,
        robust_iter=robust_iter,
        roi=roi_final.as_dict(green_x=(gx0, gx0 + gw), red_x=(rx0, rx0 + rw),
                              ref_rows=(gy0, gy0 + gh), mov_rows=(ry0, ry0 + rh),
                              seam=int(seam),
                              files=[os.path.basename(f) for f in files]),
        coarse_offset=(dy, dx),
        reference=str(reference), moving=str(moving), verbose=verbose)

    report = {
        'n_files': len(sets),
        'n_pairs': len(merged),
        'per_file_pairs': {s.files()[0] if len(s) else '?': len(s) for s in sets},
        'per_file': file_stats,          # 每文件的 对数/位移/残差
        'shift_spread': [float(np.round(spread[0], 3)),
                         float(np.round(spread[1], 3))],
        'shift_consistent': bool(consistent),
        'residual': residual_report(calib, merged, method=method),
        'pairs': merged,                 # 合并后的点, 给可视化直接画
        'pair_sets': sets,               # 每文件各自的点
    }
    if verbose:
        print('--- 分文件残差 (px)')
        for k, v in report['residual'].get('per_file', {}).items():
            print(f'    {k:28s} med {v["median"]:.3f}  p84 {v["p84"]:.3f}  '
                  f'max {v["max"]:.3f}')
    if save_path:
        calib.save(save_path)
        if verbose:
            print(f'--- 已保存 {save_path}')
    return calib, report
