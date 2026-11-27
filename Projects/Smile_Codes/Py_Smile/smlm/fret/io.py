# smlm/fret/io.py
"""FRET 结果的 CSV 落盘 / 回读。"""
from __future__ import annotations

import csv
import os
import numpy as np

COLS = ['idx', 't_bleach_s', 'frame_j', 'lag', 'x_px', 'y_px', 'n_pre', 'n_on',
        'I_A_pre', 'I_D_pre', 'I_A_post', 'I_D_post', 'S_A',
        'E', 'E_don', 'snr_A', 'gamma', 'alpha', 'beta']

# 未筛选 trace 的长表
TRACE_COLS = ['emitter_id', 'frame', 't_s', 'x_px', 'y_px',
              'I_donor', 'I_acceptor', 'bg_donor', 'bg_acceptor']


def _num(v, nd):
    v = float(v)
    return round(v, nd) if np.isfinite(v) else ''


def save_csv(results, path, cal=None, *, centres=None, verbose=True):
    """一行一个分子; γ/α/β 也写进每一行, 单个 csv 就能完整复现。

    Parameters
    ----------
    centres : dict, optional
        ``{emitter_id: (x, y)}``, 写进 ``x_px / y_px``。
    """
    rows = sorted(results, key=lambda r: r.emitter_id)
    g = float(cal.gamma) if cal is not None else np.nan
    a = float(cal.alpha) if cal is not None else np.nan
    b = float(cal.beta) if cal is not None else np.nan
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for r in rows:
            xy = (centres or {}).get(int(r.emitter_id), (np.nan, np.nan))
            w.writerow([
                int(r.emitter_id), _num(r.t_bleach, 3), int(r.frame), int(r.lag),
                _num(r.x if np.isfinite(r.x) else xy[0], 3),
                _num(r.y if np.isfinite(r.y) else xy[1], 3),
                int(r.n_pre), int(r.n_on),
                _num(r.IA_pre, 1), _num(r.ID_pre, 1),
                _num(r.IA_post, 1), _num(r.ID_post, 1), _num(r.S_A, 1),
                _num(r.E, 4), _num(r.E_don, 4), _num(r.snr_A, 2),
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
                np.asarray(td.meta.get('bg_flux', nan), float)[:n],
                np.asarray(ta.meta.get('bg_flux', nan), float)[:n]])
            np.savetxt(fh, blk, fmt=fmt, delimiter=',')
    if verbose:
        print(f'saved {len(traces_d)} traces x {n} frames -> {path}  '
              f'({os.path.getsize(path)} bytes, {len(TRACE_COLS)} cols)')
    return path


def load_csv(path):
    """回读成 dict 列表(数值列转 float, 空串转 nan)。"""
    with open(path, newline='', encoding='utf-8-sig') as f:
        rd = list(csv.reader(f))
    hdr, body = rd[0], rd[1:]
    out = []
    for row in body:
        d = {}
        for k, v in zip(hdr, row):
            if v == '':
                d[k] = np.nan
            else:
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
        out.append(d)
    return out


def verify(path, n_expect=None):
    """落盘后自检: 列数、行数、E 的分布。"""
    rows = load_csv(path)
    E = np.array([r['E'] for r in rows if np.isfinite(r.get('E', np.nan))])
    ok = (n_expect is None or len(rows) == n_expect)
    print(f'  回读 {len(rows)} 行, E 中位 {np.median(E):.4f}, '
          f'范围 {E.min():.3f}~{E.max():.3f}, E>1 {int((E > 1).sum())} 条'
          if E.size else f'  回读 {len(rows)} 行, 无有效 E')
    return ok, rows
