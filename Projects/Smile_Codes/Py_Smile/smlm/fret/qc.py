# smlm/fret/qc.py
"""逐条 trace 的质量诊断。

核心指标是 **孔径流量 / MLE 流量**:

* 孔径积分依赖环形背景, 邻居密的时候会被污染;
* MLE 光子数是高斯拟合的直接输出, **不受环形背景影响**。

两者对不上, 就说明该分子的背景或孔径不可信, 该交给人工看清楚。
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def mle_track(locs, loc_idx, n_frames):
    """精确版: 直接摊开一条 emitter 的 ``loc_idx``, 返回 (n_frames,)。"""
    y = np.full(int(n_frames), np.nan)
    idx = np.asarray(loc_idx, int)
    if idx.size == 0:
        return y
    fr = np.asarray(locs.frame, int)[idx]
    ph = np.asarray(locs.photons, float)[idx]
    ok = (fr >= 0) & (fr < y.size)
    y[fr[ok]] = ph[ok]
    return y


def mle_spatial(locs, centre_yx, *, radius=2.0, n_frames=None, tree=None):
    """逐帧取半径内离中心最近的那个定位点的 MLE photons, 没有就是 nan。

    用于**供体**通道: 孔径中心来自受体, 供体的定位要靠空间归属。
    批量调用时把 ``tree`` 传进来, 别每个分子都重建一次 KD 树。
    """
    n = int(n_frames if n_frames is not None else (np.asarray(locs.frame).max() + 1))
    y = np.full(n, np.nan)
    if locs is None or len(locs) == 0:
        return y
    xy = np.column_stack([np.asarray(locs.y, float), np.asarray(locs.x, float)])
    if tree is None:
        tree = cKDTree(xy)
    ii = np.asarray(tree.query_ball_point(np.asarray(centre_yx, float),
                                          float(radius)), int)
    if ii.size == 0:
        return y
    d = np.hypot(xy[ii, 0] - centre_yx[0], xy[ii, 1] - centre_yx[1])
    fr = np.asarray(locs.frame, int)[ii]
    ph = np.asarray(locs.photons, float)[ii]
    o = np.lexsort((d, fr))                    # 先按帧, 帧内按距离
    f, p = fr[o], ph[o]
    first = np.ones(f.size, bool)
    first[1:] = f[1:] != f[:-1]                # 每帧只留最近的一个
    ok = (f[first] >= 0) & (f[first] < n)
    y[f[first][ok]] = p[first][ok]
    return y


def flux_ratio(trace, mle, *, n_clip=3.0):
    """孔径流量 vs MLE 流量的比值诊断。

    Returns
    -------
    dict
        ``ratio``     两者的中位比值(只在两边都有值的帧上算)
        ``n``         参与估计的帧数
        ``log_spread`` log 比值的 p16~p84 宽度; 大 = 形状/背景不一致
        ``corr``      两条曲线的相关系数
    """
    ap = np.asarray(trace.intensity, float)
    ml = np.asarray(mle, float)
    n = min(ap.size, ml.size)
    ok = np.isfinite(ap[:n]) & np.isfinite(ml[:n]) & (np.abs(ml[:n]) > 1e-9)
    if ok.sum() < 5:
        return {'ratio': np.nan, 'n': int(ok.sum()),
                'log_spread': np.nan, 'corr': np.nan}
    r = ap[:n][ok] / ml[:n][ok]
    with np.errstate(all='ignore'):
        lr = np.log(np.abs(r))
        lr = lr[np.isfinite(lr)]
    corr = np.corrcoef(ap[:n][ok], ml[:n][ok])[0, 1] if ok.sum() > 2 else np.nan
    return {'ratio': float(np.median(r)),
            'n': int(ok.sum()),
            'log_spread': float(np.percentile(lr, 84) - np.percentile(lr, 16))
            if lr.size > 3 else np.nan,
            'corr': float(corr) if np.isfinite(corr) else np.nan}


def flag(diag, *, ratio_lo=0.5, ratio_hi=2.0, spread_hi=1.5, corr_lo=0.5):
    """给一条诊断打标: 返回原因字符串列表(空 = 通过)。"""
    why = []
    r = diag.get('ratio', np.nan)
    if np.isfinite(r) and not (ratio_lo <= r <= ratio_hi):
        why.append(f'孔径/MLE 比 {r:.2f} 偏离 1')
    s = diag.get('log_spread', np.nan)
    if np.isfinite(s) and s > spread_hi:
        why.append(f'两曲线形状差 {s:.2f}')
    c = diag.get('corr', np.nan)
    if np.isfinite(c) and c < corr_lo:
        why.append(f'相关 {c:.2f} 偏低')
    return why
