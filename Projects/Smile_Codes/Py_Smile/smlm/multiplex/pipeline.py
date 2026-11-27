# smlm/multiplex/pipeline.py
"""多通道配准流程: 逐文件配对 -> 合并拟合 -> 落盘标定 bundle。

典型用法见 ``smTest_dual_calibration.ipynb``::

    rois  = pair_multichannel_files(files, roi, reference='g')
    calib = fit_multichannel_pairs(rois, roi)
    MultiChannelCalibrationBundle.from_fit(calib, roi, files).save(path)

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
from ..localize import LocalMaxDetector, MLELocalizer, localize_stack
from ..io import load_raw

from .pair import PairSet, pair_multichannels, merge_pairs
from .fit import fit_transform
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
# ---------------------------------------------------------------------------
#  多文件配准
# ---------------------------------------------------------------------------
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

