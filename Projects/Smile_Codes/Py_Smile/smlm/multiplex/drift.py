# smlm/multiplex/drift.py
"""多通道漂移: 先配准点云, 再估**一条**物理漂移, 然后按各通道自己的像素尺度校正。

为什么要先配准再合并
--------------------
漂移是样品相对探测器的运动, 两个通道看到的是**同一个物理位移** —— 只是各自成像
在不同的像素坐标里(像元尺度差 ~0.1%, 再加一点剪切)。所以正确顺序是::

    各通道定位 --(配准, 坐标层)--> 参考通道坐标系 --(合并)--> 一条点云 --> 一条漂移

**绝不能**各自估一条再合并: 那会把"通道间的固定 1.4 px 错位"和"真实的逐帧漂移"
混在一起, 得到两条互不一致的轨迹。

为什么只估一条、再换算回各通道
------------------------------
合并后点密度翻倍, 估计更稳; 也顺便解决红通道分子少、单独估不稳的问题。换算是纯
线性化的: 在视场中心取雅可比 J = d(ref)/d(ch), 则 ``d_ch = J^-1 · d_ref``。
双视角的形变只有 ~0.1 px 量级, 用中心 J 的误差在 1e-3 相对量, 可以忽略。

图像只做平移
------------
``DriftCorrector`` 做的是亚像素傅里叶平移, 光度中性, 不会改变 PSF 形状。空间变化
的形变一律留给坐标层(见 ``coords.py``)。
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import numpy as np

from ..datatypes import Localizations, DriftTrace
from ..imgprocess import DriftCorrector
from .coords import DualViewGeometry


# ---------------------------------------------------------------------------
#  点云合并 (必须先配准)
# ---------------------------------------------------------------------------
_CONCAT_FIELDS = ('frame', 'x', 'y', 'photons', 'bg', 'sigma', 'lp', 'z', 'z_lp')


def merge_registered(locs, geom: DualViewGeometry, *, reference='g',
                     verbose: bool = False):
    """把各通道定位搬到**参考分析坐标**后合并成一个点云。

    Parameters
    ----------
    locs : dict
        ``{通道: Localizations}``, 坐标是该通道**分析窗口**的 local 坐标。

    Returns
    -------
    (merged, parts)
        ``merged`` 是合并后的 Localizations(参考分析坐标);
        ``parts[ch]`` 是布尔掩码, 标记 ``merged`` 里哪些点来自该通道。
    """
    names = [reference] + [c for c in locs if c != reference]
    chunks, parts, cursor = [], {}, 0
    for ch in names:
        L = locs.get(ch)
        if L is None or len(L) == 0:
            parts[ch] = np.zeros(0, bool)
            continue
        xy = geom.to_reference(np.column_stack([L.x, L.y]), ch)
        n = len(L)
        m = np.zeros(n, bool)
        m[:] = True
        parts[ch] = m
        chunks.append((replace(L, x=xy[:, 0], y=xy[:, 1]), n))
        cursor += n
    if not chunks:
        return Localizations(frame=np.zeros(0, int), x=np.zeros(0),
                             y=np.zeros(0)), {c: np.zeros(0, bool) for c in locs}

    kw = {}
    for f in _CONCAT_FIELDS:
        vals = [getattr(c, f) for c, _ in chunks]
        kw[f] = (np.concatenate([np.asarray(v) for v in vals])
                 if all(v is not None for v in vals) else None)
    merged = Localizations(**kw)

    # 掩码改成按合并后的顺序
    out, pos = {}, 0
    for ch in names:
        n = len(locs[ch]) if (ch in locs and locs[ch] is not None) else 0
        m = np.zeros(len(merged), bool)
        m[pos:pos + n] = True
        out[ch] = m
        pos += n
    if verbose:
        cnt = ', '.join(f'{c}:{int(m.sum())}' for c, m in out.items())
        print(f'  merged cloud {len(merged)} locs  ({cnt})')
    return merged, out


# ---------------------------------------------------------------------------
#  估计
# ---------------------------------------------------------------------------
@dataclass
class MultiDrift:
    """多通道漂移结果。

    traces : {通道: DriftTrace}
        各通道**自己像素坐标**下的轨迹(可直接交给 ``DriftCorrector``)。
    trace_ref : DriftTrace
        在参考通道分析坐标里估计出来的那一条(源头)。
    merged / parts
        ``merge_registered`` 的输出, 便于事后检查合并点云。
    control : DriftTrace
        只用参考通道估的对照 —— 用来核对合并有没有带偏。
    checks : dict
        自检指标。
    """
    traces: dict
    trace_ref: DriftTrace
    merged: Localizations = None
    parts: dict = field(default_factory=dict)
    control: DriftTrace = None
    checks: dict = field(default_factory=dict)
    reference: str = 'g'

    def __getitem__(self, ch):
        return self.traces[ch]

    def __repr__(self):
        d = ', '.join(f'{c}:{t.total_drift:.2f}px' for c, t in
                      self.traces.items())
        return f'MultiDrift({d})'


def _drift_to_channel(d_ref_xy, geom, ch, reference):
    """参考分析坐标下的逐帧位移 -> ``ch`` 自己像素坐标下的位移。"""
    cal = geom.calibs.get(ch)
    if ch == reference or cal is None:
        return np.asarray(d_ref_xy, float)
    # 在该通道分析窗口中心取中心雅可比; 注意 jacobian_at 吃的是**标定 local**
    _, _, h, w = geom.windows[ch]
    c_ana = np.array([w / 2.0, h / 2.0])
    p_cal = geom.to_calib_frame(c_ana[None], ch)[0]
    J = cal.jacobian_at(p_cal)                 # d(ref)/d(ch)
    return np.asarray(d_ref_xy, float) @ np.linalg.inv(J).T


def _resolve_locs_geom(data, locs, n_frames, geom):
    """兼容两种调用: ``(multi, locs, n_frames=None)`` 或 ``(locs, geom, n_frames)``。"""
    if hasattr(data, 'channels'):              # MultiImageStack
        if locs is None:
            raise TypeError('传 MultiImageStack 时还要给 locs={通道: Localizations}')
        geom = geom if geom is not None else data.geometry(check_coverage=False)
        n_frames = int(data.n_frames) if n_frames is None else int(n_frames)
        return locs, geom, n_frames
    # 旧写法: 第一个参数就是 {通道: Localizations}
    if geom is None and hasattr(locs, 'windows'):      # (locs, geom, n)
        locs, geom = data, locs
    elif geom is None:
        geom, locs = locs, data
    else:
        locs = data
    if geom is None or not hasattr(geom, 'windows'):
        raise TypeError('estimate_drift_multi 需要 geom(DualViewGeometry) '
                        '或一个 MultiImageStack')
    if n_frames is None:
        raise TypeError('缺少 n_frames')
    return locs, geom, int(n_frames)


def estimate_drift_multi(data, locs=None, n_frames=None, *,
                         frames_per_seg=10, do_stage2=True,
                         uncertainty=None, reference='g', geom=None,
                         control=True, verbose=True):
    """多通道漂移估计: 配准 -> 合并 -> 一条轨迹 -> 换算回各通道。

    Parameters
    ----------
    data : MultiImageStack
        **推荐用法**: 窗口、偏移、各通道配准都在对象里, 不用另外传几何。
    locs : dict
        ``{通道: Localizations}``, 各通道分析窗口 local 坐标。
    n_frames : int, optional
        总帧数。用 MultiImageStack 时可省略。
    geom : DualViewGeometry, optional
        不给就从 ``data`` 现造。旧写法 ``estimate_drift_multi(locs, geom, n)``
        仍然可用。
    frames_per_seg : int
        分段帧数。``CloudDriftEstimator`` 推荐 10; 用大值(如 100)会让**首段
        完全没有被校正**(Stage 1 以第 0 段为参考, 段内漂移恒为 0), 别踩。
    control : bool
        额外只用参考通道估一条, 用于核对合并是否带偏。

    Returns
    -------
    MultiDrift
    """
    from ..drift import CloudDriftEstimator

    locs, geom, n_frames = _resolve_locs_geom(data, locs, n_frames, geom)
    if not geom.calibs:
        raise ValueError('geom 里没有配准: 多通道漂移必须先有配准才能合并点云')

    merged, parts = merge_registered(locs, geom, reference=reference,
                                     verbose=verbose)

    def _estimate(L, tag):
        unc = uncertainty
        if unc is None:
            unc = float(np.median(L.lp)) if getattr(L, 'lp', None) is not None else 0.1
        est = CloudDriftEstimator(frames_per_seg=frames_per_seg,
                                  uncertainty=unc, do_stage2=do_stage2,
                                  verbose=False)
        tr = est.estimate(L, n_frames=n_frames)
        if verbose:
            print(f'  [{tag}] uncertainty={unc:.3f} px   '
                  f'drift={tr.total_drift:.2f} px (参考分析坐标)')
        return tr, unc

    trace_ref, unc = _estimate(merged, 'merged')

    # 换算到各通道自己的像素坐标
    traces = {}
    d_ref = np.column_stack([trace_ref.dx, trace_ref.dy])
    for ch in geom.windows:
        if ch != reference and (ch not in parts or parts[ch].sum() == 0):
            continue
        d_ch = _drift_to_channel(d_ref, geom, ch, reference)
        traces[ch] = DriftTrace(dx=d_ch[:, 0], dy=d_ch[:, 1],
                                pixel_size=trace_ref.pixel_size,
                                metadata={'estimator': 'merged-multi',
                                          'frame': f'{ch}-analysis'})

    checks = {'uncertainty': float(unc), 'n_merged': int(len(merged)),
              'n_per_channel': {c: int(m.sum()) for c, m in parts.items()},
              'total_drift_ref': float(trace_ref.total_drift)}

    ctrl = None
    if control and reference in locs and locs[reference] is not None:
        ctrl, _ = _estimate(locs[reference], f'{reference} only (对照)')
        d = np.hypot(trace_ref.dx - ctrl.dx, trace_ref.dy - ctrl.dy)
        checks['control_max_diff_px'] = float(np.nanmax(d))
        if verbose:
            print(f'  合并 vs 仅{reference}: max |Δdrift| = '
                  f'{checks["control_max_diff_px"]:.3f} px (应是小量)')

    return MultiDrift(traces=traces, trace_ref=trace_ref, merged=merged,
                      parts=parts, control=ctrl, checks=checks,
                      reference=reference)


# ---------------------------------------------------------------------------
#  校正
# ---------------------------------------------------------------------------
def common_crop(traces, shape, pad=4, extra_trim=0):
    """各通道**共用**的居中裁切尺寸, 以及逐通道的漂移余量。

    ``DriftCorrector`` 的裁切是居中的: 只要各通道用**同一个** (H_out, W_out),
    校正后的图像就尺寸相同、原点相同, 可以直接叠成 ``MultiImageStack``。
    所以这里取各通道所需余量的**最大值**。
    """
    H, W = int(shape[0]), int(shape[1])
    ry = rx = 0
    for tr in traces:
        ry = max(ry, int(np.ceil(np.abs(np.asarray(tr.dy) - tr.dy[0]).max())))
        rx = max(rx, int(np.ceil(np.abs(np.asarray(tr.dx) - tr.dx[0]).max())))
    ry += int(pad) + int(extra_trim)
    rx += int(pad) + int(extra_trim)
    if H - 2 * ry <= 0 or W - 2 * rx <= 0:
        raise ValueError(f'漂移 {max(ry, rx)} px + pad {pad} + extra_trim '
                         f'{extra_trim} 超过窗口 {H}x{W}, 没有公共区域可用')
    return (H - 2 * ry, W - 2 * rx), (ry, rx)


def crop_for_trace(trace, shape, pad=4):
    """单通道: 按实际漂移量算"居中裁切"的尺寸, 留出平移缓冲。"""
    return common_crop([trace], shape, pad=pad)[0]


def correct_multi(multi, locs, md: MultiDrift, geom=None, *, pad=4, edge=5,
                  extra_trim=0, verbose=True):
    """按各通道自己的轨迹校正图像 + 定位, 并返回**跟着平移的几何对象**。

    统一裁切 (关键)
    ---------------
    所有通道用同一个 ``(H_out, W_out)``, 因此校正后:

    * 各通道图像**尺寸完全相同**;
    * 局部坐标原点相同(逐像素可比对);
    * 剩下的亚像素/形变差异仍由坐标层(``geom_c``)承载。

    ``extra_trim`` 是在"漂移余量 + pad"之外再往里收的像素(漂移校正后嫌边缘脏时用)。
    它被并入裁切尺寸而不是二次裁切, 所以定位与图像始终同步。

    Returns
    -------
    (multi_c, locs_c, geom_c)
        输入 ``MultiImageStack`` -> 返回 ``MultiImageStack``;
        输入 ``{通道: ImageStack}`` -> 返回 dict(旧行为)。
    """
    is_multi = hasattr(multi, 'channels')
    if is_multi:
        stacks = {c: multi.channel(c) for c in multi.channels}
        if geom is None:
            geom = multi.geometry(check_coverage=False)
    else:
        stacks = dict(multi)
        if geom is None:
            raise TypeError('传 {通道: ImageStack} 时必须给 geom')

    chs = [c for c in stacks if c in md.traces]
    H = min(int(np.asarray(stacks[c].data).shape[-2]) for c in chs)
    W = min(int(np.asarray(stacks[c].data).shape[-1]) for c in chs)
    crop, (ry, rx) = common_crop([md.traces[c] for c in chs], (H, W),
                                 pad=pad, extra_trim=extra_trim)

    images, out_locs, new_windows = {}, {}, {}
    for ch in chs:
        tr = md.traces[ch]
        corr = DriftCorrector(crop=crop, edge=edge)
        img_c, (r0, c0, H_out, W_out) = corr.correct_image(stacks[ch], tr)
        L_c, _ = corr.correct_locs(locs[ch], tr, window=(r0, c0, H_out, W_out))
        images[ch] = img_c
        out_locs[ch] = L_c

        y0, x0 = geom.windows[ch][:2]
        new_windows[ch] = (y0 + r0, x0 + c0, H_out, W_out)
        if verbose:
            print(f'  [{ch}] crop={crop} (+{r0},{c0})   '
                  f'locs kept {len(L_c)}/{len(locs[ch])}   '
                  f'total drift {tr.total_drift:.2f} px')

    if len({(w[2], w[3]) for w in new_windows.values()}) != 1:
        raise RuntimeError('校正后各通道尺寸不一致 —— 居中裁切的前提被打破了')

    geom_c = DualViewGeometry(
        windows=new_windows, calib=dict(geom.calibs),
        calib_windows=dict(geom.calib_windows or {}),
        reference=geom.reference, roi=geom.roi, margin=0,
        check_coverage=False)          # 漂移裁切只会更靠内, 不用再查外推

    if not is_multi:
        return images, out_locs, geom_c

    data = np.stack([np.asarray(images[c].data) for c in multi.channels], axis=0)
    offs = np.array([new_windows[c][:2] for c in multi.channels], dtype=int)
    meta = dict(multi.metadata)
    meta['drift_crop'] = {'size': [int(crop[0]), int(crop[1])],
                          'margin_rc': [int(ry), int(rx)],
                          'pad': int(pad), 'extra_trim': int(extra_trim)}
    return multi.with_data(data, offsets=offs, metadata=meta), out_locs, geom_c
