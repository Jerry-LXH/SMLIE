# smlm/fret/io.py
"""FRET 结果的 CSV 落盘 (长表 / 宽表)。"""
from __future__ import annotations

import csv
import os
import numpy as np

COLS = ['idx', 't_bleach_s', 'frame_j', 'lag', 'x_px', 'y_px', 'n_pre', 'n_on',
        'I_A_pre', 'I_D_pre', 'I_A_post', 'I_D_post', 'S_A',
        'E', 'E_don', 'snr_A']

# 质量控制列 (来自 PairRecord): 下游可以据此只取 "受体平稳" 的那一批。
# 单独列出来是因为它们**可空**: 没有诊断时整列写空, 但列本身总在。
QC_COLS = ['red_label', 'red_duty', 'spread', 'n_d_pre']

COLS = COLS + QC_COLS + ['gamma', 'alpha', 'beta']

# 未筛选 trace 的长表
TRACE_COLS = ['emitter_id', 'frame', 't_s', 'x_px', 'y_px',
              'I_donor', 'I_acceptor', 'bg_donor', 'bg_acceptor']


def _num(v, nd):
    v = float(v)
    return round(v, nd) if np.isfinite(v) else ''


def save_csv(results, path, cal=None, *, centres=None, qc=None, verbose=True):
    """一行一个分子; γ/α/β 也写进每一行, 单个 csv 就能完整复现。

    Parameters
    ----------
    centres : dict, optional
        ``{emitter_id: (x, y)}``, 写进 ``x_px / y_px``。
    qc : dict, optional
        ``{emitter_id: {列名: 值}}``, 填 :data:`QC_COLS` 那几列。缺的分子写空
        (列本身总在, 下游按名字读就不会因为有没有 qc 而变)。
    """
    rows = sorted(results, key=lambda r: r.emitter_id)
    g = float(cal.gamma) if cal is not None else np.nan
    a = float(cal.alpha) if cal is not None else np.nan
    b = float(cal.beta) if cal is not None else np.nan
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    qc = qc or {}

    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for r in rows:
            rid = int(r.emitter_id)
            xy = (centres or {}).get(rid, (np.nan, np.nan))
            extra = qc.get(rid, {})
            ndp = extra.get('n_d_pre', -1)
            ndp = '' if (ndp is None or int(ndp) < 0) else int(ndp)
            w.writerow([
                rid, _num(r.t_bleach, 3), int(r.frame), int(r.lag),
                _num(r.x if np.isfinite(r.x) else xy[0], 3),
                _num(r.y if np.isfinite(r.y) else xy[1], 3),
                int(r.n_pre), int(r.n_on),
                _num(r.IA_pre, 1), _num(r.ID_pre, 1),
                _num(r.IA_post, 1), _num(r.ID_post, 1), _num(r.S_A, 1),
                _num(r.E, 4), _num(r.E_don, 4), _num(r.snr_A, 2),
                extra.get('red_label', ''), _num(extra.get('red_duty', np.nan), 3),
                _num(extra.get('spread', np.nan), 2), ndp,
                _num(g, 4), _num(a, 4), _num(b, 2)])

    if verbose:
        print(f'saved {len(rows)} rows -> {path}  '
              f'({os.path.getsize(path)} bytes, {len(COLS)} cols)')
    return path


def save_traces_csv(traces_d, traces_a, t_s, path, verbose=True):
    """双通道 trace 的长表: 一行 = 一个分子的一帧。

    ``traces_d`` / ``traces_a`` 必须按 emitter_id 对齐(供体=绿, 受体=红)。
    未扣背景的积分强度可自己还原: ``raw = I + bg``。

    传全部分子就是 ``traces_all.csv``; 只传通过判据的那些就是 ``good_traces.csv``
    —— 两者格式完全相同, 下游读法一样。
    """
    n = int(len(t_s))
    nan = np.full(n, np.nan)
    fmt = ['%.0f', '%.0f', '%.3f', '%.2f', '%.2f', '%.1f', '%.1f', '%.1f', '%.1f']

    def _bg(trace, key):
        v = getattr(trace, key, None)
        return nan if v is None else np.asarray(v, float)[:n]

    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, 'w', encoding='utf-8-sig', newline='') as fh:
        fh.write(','.join(TRACE_COLS) + '\n')
        for td, ta in zip(traces_d, traces_a):
            blk = np.column_stack([
                np.full(n, int(td.emitter_id)), np.arange(n),
                np.asarray(t_s, float),
                np.full(n, float(td.x)), np.full(n, float(td.y)),
                np.asarray(td.intensity, float)[:n],
                np.asarray(ta.intensity, float)[:n],
                _bg(td, 'bg_flux'), _bg(ta, 'bg_flux')])
            np.savetxt(fh, blk, fmt=fmt, delimiter=',')
    if verbose:
        print(f'saved {len(traces_d)} traces x {n} frames -> {path}  '
              f'({os.path.getsize(path)} bytes, {len(TRACE_COLS)} cols)')
    return path


def load_traces_csv(path, verbose=False):
    """把 :func:`save_traces_csv` 的长表读回成按分子分块的数组。

    Returns
    -------
    dict
        ``{'emitter_id': (n_mol,), 'frame': (n_mol, n_frames), 't_s': ...,
        'x_px': ..., 'y_px': ..., 'I_donor': ..., 'I_acceptor': ...,
        'bg_donor': ..., 'bg_acceptor': ...}``

    要求同一文件里所有分子帧数相同 (``save_traces_csv`` 写出来就满足);
    不一致时按到达顺序无法安全 reshape, 会直接报错并提示改用原始长表。
    """
    arr = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    ids = arr[:, 0].astype(int)
    order = np.lexsort((arr[:, 1], ids))
    arr, ids = arr[order], ids[order]
    uniq, counts = np.unique(ids, return_counts=True)
    if counts.size and counts.min() != counts.max():
        raise ValueError(
            f'{path}: 各分子帧数不一致 {counts.min()}..{counts.max()}, '
            f'不能按块 reshape; 请直接用长表 (一行 = 一分子一帧)')
    n = int(counts[0]) if counts.size else 0
    out = {'emitter_id': uniq.astype(int)}
    for k, name in enumerate(TRACE_COLS[1:], start=1):
        out[name] = arr[:, k].reshape(uniq.size, n)
    if verbose:
        print(f'loaded {uniq.size} traces x {n} frames <- {path}')
    return out

