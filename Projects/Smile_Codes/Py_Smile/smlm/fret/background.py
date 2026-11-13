# smlm/fret/background.py
"""背景估计: 圆环 + 邻居核心挖洞 + sigma-clip 中位数。

沿用原 ``ApertureTraceExtractor`` 的口径(暂时不改), 但拆成纯函数, 并把
"背景到底估计成了多少"作为一等输出暴露出来 —— 复核界面要把它画进 trace。

为什么只挖邻居核心(``excl_radius=2``)
-------------------------------------
不挖 -> 邻居的外周抬高环内背景 -> 背景高估 -> trace 偏低 -> E 偏高。
挖太大 -> 把"最能代表本分子位置上背景水平"的那部分像素也删了 -> 背景低估。
所以只避开邻居的峰值核心, 外周仍然进入估计。
"""
from __future__ import annotations

import warnings

import numpy as np


# ---------------------------------------------------------------------------
#  环掩膜
# ---------------------------------------------------------------------------
def ring_mask(n, half, radius, spec, *, neighbours_rel=None, valid=None,
              center_frac=(0.0, 0.0)):
    """在 (n,n) 的块上取背景环, 返回 bool 掩膜。

    Parameters
    ----------
    n : int
        块边长(奇数), 中心像素在 ``[half, half]``。
    radius : int
        信号孔径半径。环从 ``radius + bg_gap`` 起算(严格 >)。
    neighbours_rel : (K,2) ndarray, optional
        邻居在**块坐标**里的 (row, col); 每个邻居周围 ``excl_radius`` 内不取。
    valid : (n,n) bool, optional
        图内有效像素(越界补 0 的地方要剔掉, 否则背景被往下拉)。
    center_frac : (fy, fx)
        中心相对 ``[half, half]`` 的小数偏移。旧实现用**浮点**中心算距离,
        这里保持一致, 否则环形会整体偏半个像素。
    """
    fy, fx = float(center_frac[0]), float(center_frac[1])
    cy, cx = half + fy, half + fx
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(yy - cy, xx - cx)
    r_in = float(radius + spec.bg_gap)
    r_out = float(radius + spec.bg_radius)
    m = (rr > r_in) & (rr <= r_out)

    if neighbours_rel is not None and len(neighbours_rel):
        nr = np.asarray(neighbours_rel, float).reshape(-1, 2)
        for (ry, cxx) in nr:
            m &= np.hypot(yy - ry, xx - cxx) > float(spec.excl_radius)

    if valid is not None:
        m &= np.asarray(valid, bool)
    return m


def valid_mask(half, r0, c0, image_shape):
    """块内哪些像素真的落在图里(crop_centered 越界会补 0)。"""
    H, W = int(image_shape[0]), int(image_shape[1])
    n = 2 * int(half) + 1
    ys = np.arange(r0, r0 + n)
    xs = np.arange(c0, c0 + n)
    return ((ys >= 0) & (ys < H))[:, None] & ((xs >= 0) & (xs < W))[None, :]


# ---------------------------------------------------------------------------
#  逐帧 sigma-clip 中位数
# ---------------------------------------------------------------------------
def clip_median(v, n_clip=3.0, n_iter=3):
    """沿 axis=1 做 sigma-clip 中位数。``v`` 是 (T, M), 允许 nan。

    逐步与旧 ``ApertureTraceExtractor._sigma_clip_rows`` 一致:
      * ``std < 1e-10`` 的帧不裁剪(否则常数列会被削成空);
      * 有效像素 < 5 的帧退回**未裁剪**的中位数。
    """
    v = np.asarray(v, float)
    if v.ndim != 2:
        raise ValueError(f'clip_median 需要 (T,M), 得到 {v.shape}')
    T, M = v.shape
    if M == 0:
        return np.full(T, np.nan)

    x = np.where(np.isfinite(v), v, np.nan)
    for _ in range(int(n_iter)):
        with np.errstate(all='ignore'):
            med = np.nanmedian(x, axis=1, keepdims=True)
            std = np.nanstd(x, axis=1, keepdims=True)
        std = np.where(std < 1e-10, np.inf, std)
        keep = np.abs(np.nan_to_num(x - med, nan=0.0)) <= float(n_clip) * std
        new = np.where(keep, x, np.nan)
        if np.array_equal(np.isnan(new), np.isnan(x)):
            break
        x = new

    with np.errstate(all='ignore'):
        n_keep = np.isfinite(x).sum(axis=1)
        out = np.nanmedian(x, axis=1)
        fb = n_keep < 5
        if fb.any():
            out[fb] = np.nanmedian(x[fb], axis=1) if x[fb].size else np.nan
    return out


# ---------------------------------------------------------------------------
#  一个分子的背景估计
# ---------------------------------------------------------------------------
def moving_average(x, w):
    """居中滑动平均, 边缘用**缩短的窗口**(不是补零), nan 自动跳过。"""
    x = np.asarray(x, float)
    w = int(w)
    if w <= 1 or x.size == 0:
        return x.copy()
    ok = np.isfinite(x)
    k = np.ones(w)
    num = np.convolve(np.where(ok, x, 0.0), k, mode='same')
    den = np.convolve(ok.astype(float), k, mode='same')
    with np.errstate(all='ignore'):
        return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)


def _drift(x, frac=0.1):
    """背景时间序列漂了多少: 首尾各 ``frac`` 段的中位。

    诊断用: 若 ``ratio = early/late`` 明显 >1, 说明**早期背景更高** —— 常见于
    前期还有很多分子亮着, 它们的 PSF 尾巴漏进了环里。此时"末端聚合常数"会比
    "逐帧"给出更干净的背景, 但也会漏掉真实的时间漂移。
    """
    x = np.asarray(x, float)
    v = x[np.isfinite(x)]
    if v.size < 20:
        return {}
    k = max(3, int(v.size * float(frac)))
    early, late = float(np.median(v[:k])), float(np.median(v[-k:]))
    return {'early': early, 'late': late,
            'ratio': (early / late) if late not in (0.0, np.nan) else np.nan}


def model_weights(ring_yx, neighbour_yx, sigma):
    """邻居 PSF 落在**每个环像素**上的通量占比。

    Parameters
    ----------
    ring_yx : (M,2)
        环像素在块坐标里的位置。
    neighbour_yx : (K,2)
        邻居在块坐标里的位置。
    sigma : (K,)
        邻居的 PSF 宽度。

    Returns
    -------
    (M, K) ndarray
        ``W[m, j]`` = 邻居 j 落在像素 m 上的光子占其总数的比例。
    """
    ring_yx = np.asarray(ring_yx, float).reshape(-1, 2)
    nb = np.asarray(neighbour_yx, float).reshape(-1, 2)
    if ring_yx.size == 0 or nb.size == 0:
        return np.zeros((len(ring_yx), len(nb)))
    s2 = np.maximum(np.asarray(sigma, float) ** 2, 0.25)      # σ 下限 0.5 px
    d2 = ((ring_yx[:, None, :] - nb[None, :, :]) ** 2).sum(-1)
    return np.exp(-d2 / (2.0 * s2[None, :])) / (2.0 * np.pi * s2)[None, :]


def emitter_frame_photons(emitters, locs, n_frames, n_expect=None):
    """每条 emitter 的**逐帧光子数**, 按 ``loc_idx`` 回查 ``locs``。

    Returns
    -------
    (T, K) ndarray; 同帧多条就累加, 没检出就是 0。
    """
    K = len(emitters) if n_expect is None else int(n_expect)
    out = np.zeros((int(n_frames), K), np.float32)
    if locs is None or not len(emitters):
        return out
    fr_all = np.asarray(locs.frame, int)
    ph_all = np.asarray(locs.photons, float)
    for j, e in enumerate(emitters):
        idx = np.asarray(getattr(e, 'loc_idx', []), int)
        if idx.size == 0:
            continue
        fr, ph = fr_all[idx], ph_all[idx]
        ok = (fr >= 0) & (fr < n_frames)
        if ok.any():
            np.add.at(out[:, j], fr[ok], ph[ok])
    return out


def emitter_sigma(emitters, locs, default=1.7):
    """每条 emitter 的 σ(取它自己定位点的中位 σ)。"""
    out = np.full(len(emitters), float(default))
    if locs is None or getattr(locs, 'sigma', None) is None:
        return out
    sg = np.asarray(locs.sigma, float)
    for j, e in enumerate(emitters):
        idx = np.asarray(getattr(e, 'loc_idx', []), int)
        if idx.size:
            v = sg[idx]
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                out[j] = float(np.median(v))
    return out


class AnnulusBackground:
    """某分子某通道的背景环。构造时定好掩膜, ``estimate`` 只做数值。"""

    def __init__(self, spec, radius, half, *, neighbours_rel=None, valid=None,
                 center_frac=(0.0, 0.0)):
        self.spec = spec
        self.radius = int(radius)
        self.half = int(half)
        n = 2 * self.half + 1
        self.mask = ring_mask(n, self.half, self.radius, spec,
                              neighbours_rel=neighbours_rel, valid=valid,
                              center_frac=center_frac)
        self.n_pixels = int(self.mask.sum())

    def estimate(self, box, subtract=None):
        """``box`` (T,n,n) -> (bg_per_pixel (T,), info dict)。

        顺序: (可选)减掉邻居模型 -> 逐帧 sigma-clip 中位数 -> (可选)滑动平均
        -> (可选)末端聚合成常数。

        Parameters
        ----------
        subtract : (T, M) ndarray, optional
            **邻居落在环像素上的逐帧光子数**。直接减在像素值上, 而不是把那些
            像素挖掉 —— 早期邻居还亮着时, 它们的 PSF 翅膀会把环抬高, 但那股光
            并不在本分子的孔径里, 逐帧扣就会过度扣除(把 I_D^pre 推成负数)。
        """
        box = np.asarray(box, float)
        if self.n_pixels == 0:
            T = box.shape[0]
            return np.full(T, np.nan), {'n_pixels': 0}
        vals = box[:, self.mask]                       # (T, M)
        n_model = 0
        sub_mean = 0.0
        sub_rel = 0.0
        if subtract is not None:
            sub = np.asarray(subtract, float)
            if sub.shape == vals.shape:
                vals = vals - sub
                n_model = int(np.isfinite(sub).any())
                # 扣除量: 每帧在环像素上的平均(photon/px), 以及相对环中位数的比例
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    sub_mean = float(np.nanmedian(np.nanmean(sub, axis=1)))
                    base = float(np.nanmedian(box[:, self.mask]))
                    sub_rel = sub_mean / base if abs(base) > 1e-9 else 0.0
        raw = clip_median(vals, self.spec.clip_n, self.spec.clip_iter)

        bg = raw
        if self.spec.bg_smooth and int(self.spec.bg_smooth) > 1:
            bg = moving_average(bg, int(self.spec.bg_smooth))

        info = {'n_pixels': self.n_pixels, 'ref': self.spec.bg_ref,
                'smooth': int(self.spec.bg_smooth), 'model': bool(n_model),
                'sub_mean': sub_mean, 'sub_rel': sub_rel}
        if self.spec.bg_ref == 'tail':
            k = int(min(self.spec.bg_tail, bg.size))
            const = float(np.nanmedian(bg[-k:])) if k > 0 else np.nan
            info['tail_value'] = const
            bg = np.full(bg.size, const)
        info['drift'] = _drift(raw)
        return bg, info

    def __repr__(self):
        return (f'AnnulusBackground(r={self.radius}, ring '
                f'{self.radius + self.spec.bg_gap}~'
                f'{self.radius + self.spec.bg_radius}, '
                f'{self.n_pixels} px)')


# ---------------------------------------------------------------------------
#  挖洞预算 (不读图像)
# ---------------------------------------------------------------------------
def ring_budget(centers_yx, spec, *, radius, mask_centers_yx=None,
                image_shape=None, excl_list=(2, 3, 4, 5), self_radius=1.5):
    """估算"挖洞半径 -> 还剩多少背景像素"。

    背景环本来就薄(默认 r+2..r+6, 只有 ~180 px), 挖洞半径一加大就会掉像素。
    在改 ``excl_radius`` 之前先看一眼这个表, 免得把环挖空。

    Parameters
    ----------
    centers_yx : (N,2)
        孔径中心 ``[row, col]``。
    radius : int
        信号半径(硬孔径就是固定值)。
    mask_centers_yx : (K,2), optional
        要挖掉的中心; ``None`` = 用 ``centers_yx`` 自己。

    Returns
    -------
    dict
        ``{excl_radius: (中位数, p10, 最小)}``, 单位像素。
    """
    from dataclasses import replace as _replace
    from scipy.spatial import cKDTree
    from .aperture import crop_centered  # noqa: F401  (保持依赖显式)

    c = np.asarray(centers_yx, float).reshape(-1, 2)
    self_mask = mask_centers_yx is None
    mc = c if self_mask else np.asarray(mask_centers_yx, float).reshape(-1, 2)
    tree = cKDTree(mc) if len(mc) else None

    half = int(radius) + int(spec.bg_radius)
    n = 2 * half + 1
    out = {}
    for excl in excl_list:
        sp = _replace(spec, excl_radius=int(excl))
        counts = []
        for (cy, cx) in c:
            ir, ic = int(round(cy)), int(round(cx))
            valid = None
            if image_shape is not None:
                valid = valid_mask(half, ir - half, ic - half, image_shape)
            nrel = None
            if tree is not None and len(mc):
                qr = half + float(excl) + 1.0
                idx = np.asarray(tree.query_ball_point([cy, cx], qr), int)
                near = mc[idx].astype(float) if idx.size else np.zeros((0, 2))
                if len(near):
                    d_self = np.hypot(near[:, 0] - cy, near[:, 1] - cx)
                    j = int(np.argmin(d_self))
                    if self_mask or d_self[j] <= self_radius:
                        near = np.delete(near, j, axis=0)
                if len(near):
                    nrel = near + np.array([half - ir, half - ic], float)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                m = ring_mask(n, half, int(radius), sp,
                              neighbours_rel=nrel, valid=valid,
                              center_frac=(cy - ir, cx - ic))
            counts.append(int(m.sum()))
        a = np.asarray(counts, int)
        out[int(excl)] = (float(np.median(a)) if a.size else 0.0,
                          float(np.percentile(a, 10)) if a.size else 0.0,
                          int(a.min()) if a.size else 0)
    return out
