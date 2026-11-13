# smlm/multiplex/fit.py
"""从配对点拟合通道间变换 (poly / affine), 以及拟合质量的报告。

稳健流程(与旧 notebook 一致):
    拟合 -> 算残差 -> 用 median + 4*MAD 剔除错配 -> 再拟合, 最多 n_iter 轮。
错配点主要来自两个分子恰好落在同一个搜索半径内, 用 4*MAD 这种对离群不敏感的
尺度估计比用标准差稳得多。
"""
from __future__ import annotations

import numpy as np

from .transform import DualViewCalibration
from .pair import PairSet

# 每参数建议的最少配对数, 低于它就自动降阶 (10 参数至少要 20 对才不至于过拟合)
MIN_PAIRS_PER_PARAM = 2.0


def poly_terms(order):
    """完整 order 阶二维多项式的项 (i = x 幂, j = y 幂), i + j <= order。

    **必须**与 DualViewCalibration._terms() 完全一致, 否则旧 npz 的系数会错位。
    """
    return [(i, d - i) for d in range(int(order) + 1) for i in range(d, -1, -1)]


def n_params(order):
    return len(poly_terms(order))


def design(xy, center, scale, order=3):
    """在归一化坐标上构造完整 order 阶二维多项式设计矩阵。"""
    n = (np.asarray(xy, float).reshape(-1, 2) - np.asarray(center, float)) / float(scale)
    x, y = n[:, 0], n[:, 1]
    return np.column_stack([x ** i * y ** j for i, j in poly_terms(order)])


def _auto_order(n, min_pairs_per_param=MIN_PAIRS_PER_PARAM, max_order=3):
    for o in range(int(max_order), 0, -1):
        if n >= min_pairs_per_param * n_params(o):
            return o
    return 1


def _fit_poly(R, G, center, scale, order):
    X = design(R, center, scale, order)
    a = np.linalg.lstsq(X, G[:, 0], rcond=None)[0]
    b = np.linalg.lstsq(X, G[:, 1], rcond=None)[0]
    return np.column_stack([X @ a, X @ b]), a, b


def _fit_affine(R, G):
    A = np.column_stack([R, np.ones(len(R))])
    aff = np.linalg.lstsq(A, G, rcond=None)[0]
    return A @ aff, aff


# ---------------------------------------------------------------------------
#  用交叉验证决定多项式阶数
# ---------------------------------------------------------------------------
def cv_select_order(R, G, center, scale, max_order=3, min_order=1,
                    n_folds=10, max_samples=4000, cv_tol=1.05,
                    cv_metric='cv_p84_px', rng=None):
    """留一/留K交叉验证选多项式阶数, 避免"点少 + 参数多"的过拟合。

    为什么必须做: 标准样品每幅只配出十几个点, 3 阶多项式有 10 个参数 ——
    此时训练残差可以做到 0.05 px, 但视场边缘的外推完全是噪声。交叉验证直接
    回答"多出来的那阶非线性到底是不是真的"。

    为什么用 p84 而不是中位数: 实测 26 个点拟合时, 3 阶的 CV **中位数**比 1 阶
    还好一点(0.085 vs 0.101 px), 但 p84 从 0.139 恶化到 0.267 px —— 多出来的
    参数只是在拟合噪声, 代价全压在尾部。而尾部误差正是会污染 FRET 的那部分。

    选阶规则: 取 cv_metric 最小的阶, 然后在"误差不超过 best*cv_tol"的阶里
    选**最低**的那个(偏好简单模型, 外推更稳)。

    Returns
    -------
    (best_order, table)
        table = [{'order', 'n_params', 'cv_median_px', 'cv_p84_px',
                  'cv_rms_px', 'cv_max_px'}, ...]
    """
    R = np.asarray(R, float)
    G = np.asarray(G, float)
    n = len(R)
    if rng is None:
        rng = np.random.default_rng(0)

    # 点太多就抽样, 控制 CV 代价 (CV 的结论不依赖于点数的绝对多少)
    if n > int(max_samples):
        idx = rng.choice(n, int(max_samples), replace=False)
        R, G, n = R[idx], G[idx], int(max_samples)

    if n <= 3:
        order = max(int(min_order), 1)
        return order, [{'order': order, 'n_params': n_params(order),
                        'cv_median_px': float('nan'), 'cv_p84_px': float('nan'),
                        'cv_rms_px': float('nan'), 'cv_max_px': float('nan'),
                        'note': 'too few pairs'}]

    # 小样本用留一, 大样本用 K 折
    folds = n if n <= 60 else int(min(n_folds, n))
    kf = np.array_split(rng.permutation(n), folds)

    table = []
    for o in range(int(min_order), int(max_order) + 1):
        X = design(R, center, scale, o)
        p = X.shape[1]
        if n - len(kf[0]) < p:                 # 训练集都装不下参数, 不可能
            table.append({'order': o, 'n_params': int(p),
                          'cv_median_px': float('inf'),
                          'cv_p84_px': float('inf'),
                          'cv_rms_px': float('inf'),
                          'cv_max_px': float('inf'),
                          'note': 'not enough pairs'})
            continue
        errs = []
        for test in kf:
            train = np.setdiff1d(np.arange(n), test)
            a = np.linalg.lstsq(X[train], G[train, 0], rcond=None)[0]
            b = np.linalg.lstsq(X[train], G[train, 1], rcond=None)[0]
            e = np.hypot(X[test] @ a - G[test, 0], X[test] @ b - G[test, 1])
            errs.append(e)
        e = np.concatenate(errs)
        table.append({'order': o, 'n_params': int(p),
                      'cv_median_px': float(np.round(np.median(e), 4)),
                      'cv_p84_px': float(np.round(np.percentile(e, 84), 4)),
                      'cv_rms_px': float(np.round(np.sqrt((e ** 2).mean()), 4)),
                      'cv_max_px': float(np.round(e.max(), 4))})

    valid = [t for t in table if np.isfinite(t.get(cv_metric, np.inf))]
    if not valid:
        return int(min_order), table
    best = min(t[cv_metric] for t in valid)
    thr = best * float(cv_tol)
    order = min(t['order'] for t in valid if t[cv_metric] <= thr)
    return int(order), table


def _rms(e):
    return float(np.sqrt((np.asarray(e, float) ** 2).sum(axis=1).mean()))


def _summary(pred, G):
    e = pred - G
    r = np.hypot(e[:, 0], e[:, 1])
    return {'rms': float(np.round(_rms(e), 4)),
            'median': float(np.round(np.median(r), 4)),
            'p84': float(np.round(np.percentile(r, 84), 4)),
            'max': float(np.round(r.max(), 4))}


def fit_transform(pairs, *, order='auto', method='poly', max_order=3,
                  robust_k=4.0, robust_iter=4, min_abs_px=0.02,
                  min_pairs_per_param=MIN_PAIRS_PER_PARAM,
                  cv_tol=1.05, cv_metric='cv_p84_px', n_folds=10,
                  roi=None, coarse_offset=(0, 0),
                  reference='g', moving='r', verbose=True):
    """从 PairSet 拟合 mov -> ref 的变换。

    Parameters
    ----------
    pairs : PairSet
        要求的位移量是 ``ref - mov`` (多项式直接拟合这个映射)。
    order : int 或 'auto'
        'auto' 用交叉验证选阶 (见 ``cv_select_order``)。标准样品常常每幅只配出
        十几个点, 硬套 3 阶多项式会得到"训练残差很小、外推全是噪声"的结果。
    method : 'poly' | 'affine'
        决定 ``metrics`` 里以哪个为主, 但两个都会算并都写进 metrics。
    robust_k, robust_iter : 迭代剔除阈值 = median + robust_k * MAD, MAD 下限 min_abs_px。
    cv_tol : 'auto' 选阶的容差, 误差在 best*cv_tol 内的最低阶优先。

    Returns
    -------
    DualViewCalibration
    """
    if len(pairs) == 0:
        raise ValueError('配对数为 0, 无法拟合')

    R_raw = np.asarray(pairs.mov, float)
    G = np.asarray(pairs.ref, float)
    # 注意: coarse_offset 不参与拟合。两通道各自从半幅起点裁切, local 坐标本来就
    # 重合(256 px 平移已被裁切吸收), 多项式只需拟合残差形变 —— 见 transform.py。
    R = R_raw
    _dy, _dx = coarse_offset

    center = R.mean(axis=0)
    scale = max(float(np.ptp(R, axis=0).max()), 1.0) / 2.0

    # ---------- 定阶 ----------
    cv_table = None
    if order == 'auto':
        ord_req, cv_table = cv_select_order(
            R, G, center, scale, max_order=max_order, cv_tol=cv_tol,
            cv_metric=cv_metric, n_folds=n_folds)
        if verbose:
            print('  交叉验证选阶 (指标 %s):' % cv_metric)
            for t in cv_table:
                print(f'    order {t["order"]}  ({t["n_params"]:2d} 参数)  '
                      f'med {t["cv_median_px"]:.3f}  p84 {t["cv_p84_px"]:.3f}  '
                      f'rms {t["cv_rms_px"]:.3f}  max {t["cv_max_px"]:.3f} px'
                      + (f'   [{t["note"]}]' if 'note' in t else ''))
            print(f'    -> 选 order {ord_req}')
        # 硬约束: 参数个数不能超过 min_pairs_per_param 允许的上限
        while ord_req > 1 and len(R) < min_pairs_per_param * n_params(ord_req):
            ord_req -= 1
    else:
        ord_req = min(int(order), int(max_order))
        while ord_req > 1 and len(R) < min_pairs_per_param * n_params(ord_req):
            ord_req -= 1

    # ---------- 稳健迭代 (只围绕多项式做, affine 用最后一轮的 keep) ----------
    keep = np.ones(len(R), bool)
    history = []
    for it in range(int(robust_iter)):
        if keep.sum() < n_params(ord_req):
            keep[:] = True
            break
        pred, _, _ = _fit_poly(R[keep], G[keep], center, scale, ord_req)
        e = np.hypot(*(pred - G[keep]).T)
        mad = 1.4826 * np.median(np.abs(e - np.median(e)))
        thr = float(np.median(e) + robust_k * max(mad, min_abs_px))
        new = e < thr
        history.append({'iter': it, 'kept': int(keep.sum()),
                        'threshold_px': float(np.round(thr, 3)),
                        'median_px': float(np.round(np.median(e), 3))})
        if verbose:
            print(f'  iter {it}: kept {keep.sum():5d}/{len(keep):5d}   '
                  f'thr={thr:.2f} px   median={np.median(e):.3f} px')
        if new.all():
            break
        idx = np.flatnonzero(keep)[new]
        keep = np.zeros_like(keep)
        keep[idx] = True
    # 剔除过狠就全部回退(否则会拿 3 个点拟合 10 个参数)
    if keep.sum() < max(n_params(ord_req),
                        int(min_pairs_per_param * n_params(ord_req))):
        keep[:] = True

    order_used = _auto_order(int(keep.sum()), min_pairs_per_param, ord_req)

    # ---------- 最终拟合 ----------
    Rk, Gk = R[keep], G[keep]
    pred_poly, poly_a, poly_b = _fit_poly(Rk, Gk, center, scale, order_used)
    pred_aff, aff = _fit_affine(Rk, Gk)

    metrics = {
        'n_pairs_raw': int(len(R)),
        'n_pairs_kept': int(keep.sum()),
        'poly_order': int(order_used),
        'poly_order_requested': int(ord_req),
        'poly': _summary(pred_poly, Gk),
        'affine': _summary(pred_aff, Gk),
        'fit_method': str(method),
        'robust': history,
        'coverage': {'x': [float(R_raw[:, 0].min()), float(R_raw[:, 0].max())],
                     'y': [float(R_raw[:, 1].min()), float(R_raw[:, 1].max())]},
        'coarse_offset_dy_dx': [float(_dy), float(_dx)],
        'files': list(pairs.files()),
    }
    if cv_table is not None:
        metrics['cv'] = cv_table
        metrics['cv_tol'] = float(cv_tol)
        metrics['cv_metric'] = str(cv_metric)
    # 分文件残差: 哪个文件在拖后腿, 一眼就能看出来
    per_file = {}
    for f in pairs.files():
        m = (pairs.file == f)                 # 全量配对上的掩码
        sel = m[keep]                         # 只保留稳健循环留下的点
        if sel.sum() < 3:
            continue
        per_file[str(f)] = {'n': int(m.sum()), 'n_kept': int(sel.sum())}
        per_file[str(f)].update(_summary(pred_poly[sel], Gk[sel]))
    if per_file:
        metrics['per_file'] = per_file

    if verbose:
        print(f'  多项式阶数 = {order_used}  ({n_params(order_used)} 参数)')
        print(f'  poly  : median {metrics["poly"]["median"]:.3f} px   '
              f'rms {metrics["poly"]["rms"]:.3f} px')
        print(f'  affine: median {metrics["affine"]["median"]:.3f} px   '
              f'rms {metrics["affine"]["rms"]:.3f} px')

    return DualViewCalibration(
        poly_a=poly_a, poly_b=poly_b, center=center, scale=float(scale),
        affine=aff, poly_order=int(order_used),
        roi=dict(roi or {}), metrics=metrics,
        coarse_offset=tuple(coarse_offset),
        reference=str(reference), moving=str(moving))


def residual_report(calib, pairs, method='poly'):
    """给定标定和配对点, 返回总体 + 分文件的残差统计。"""
    R = np.asarray(pairs.mov, float)
    G = np.asarray(pairs.ref, float)
    pred = calib.transform(R, method=method)
    out = {'all': _summary(pred, G)}
    if len(pairs) and len(pairs.files()) > 1:
        pf = {}
        for f in pairs.files():
            m = pairs.file == f
            if m.sum() >= 3:
                pf[str(f)] = _summary(pred[m], G[m])
        out['per_file'] = pf
    return out
