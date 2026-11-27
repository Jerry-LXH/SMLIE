# smlm/fret/export.py
"""把一次 FRET 分析的全部成品落盘到一个目录。

目录结构 (与 ``smTest_dual_batch.ipynb`` 的约定一致)::

    <sif 主名>_fret/
        fret.csv             逐分子 E 与 γ/α/β —— 下游统计只读它
        E_distribution.png   E 分布 + 两个交叉检查
        good_traces_all.png  通过判据的 trace 全览
        traces_all.csv       全部 trace (一行 = 一分子一帧), 可读回
        good_traces.csv      通过判据的那些 trace (同格式)

设计: 这里只做**编排**, 写表在 :mod:`.io`, 画图在 :mod:`.viz`。
所以"换个画法"或"多写一张表"都不必动这个文件。
"""
from __future__ import annotations

import os

import numpy as np

from .io import save_csv, save_traces_csv
from .viz import plot_efficiency, plot_efficiency_vs_brightness, plot_pair_gallery

__all__ = ['save_bundle', 'FILES']

# 目录里的固定文件名 (下游脚本会按名字找, 不要改)
FILES = {
    'fret_csv': 'fret.csv',
    'e_png': 'E_distribution.png',
    'gallery_png': 'good_traces_all.png',
    'traces_all': 'traces_all.csv',
    'traces_good': 'good_traces.csv',
}

_E_MED = 'E 中位 {med:.3f} (n={n})'


def save_bundle(out_dir, *, records, results, traces_d, traces_a, cal=None,
                ex_time=1.0, title='', dpi=300, ncols=6, include_all=True,
                pre_win=100, post_win=400, verbose=True):
    """把一个视野的成品整体写进 ``out_dir`` (不存在就建)。

    Parameters
    ----------
    records : list[PairRecord]
        :func:`~smlm.fret.pairs.select_pairs` 的完整输出 (合格与否都在里面)。
    results : list[FretResult]
        算出了 E 的那些分子 (可能比合格对少: 有的对凑不出窗口)。
    traces_d, traces_a : list[Trace]
        供体(绿) / 受体(红) 的全部 trace, **按 emitter_id 对齐**。
    cal : FretCalibration, optional
        γ/α/β 会写进 ``fret.csv`` 的每一行, 单个表就能复现。
    ex_time : float
        每帧秒数 (csv 的 ``t_s`` 与图的横轴都用它)。
    title : str
        图标题前缀 (通常是 sif 文件名)。
    include_all : bool
        是否连**未筛选**的 ``traces_all.csv`` 一起写。一个 760 分子 x 2000 帧的
        视野大约 76 MB, 只关心合格分子时关掉能省一大截。

    Returns
    -------
    dict
        实际写出的 ``{名字: 路径}``; 因数据为空而跳过的项不出现。
    """
    import matplotlib.pyplot as plt

    out_dir = os.fspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    if len(traces_d) != len(traces_a):
        raise ValueError(f'两通道 trace 数不一致: {len(traces_d)} vs {len(traces_a)}')

    written = {'dir': out_dir}
    path = lambda key: os.path.join(out_dir, FILES[key])      # noqa: E731

    # ---- 1) 全帧长表: 全部 trace / 合格 trace ----
    n_frames = min((len(t.intensity) for t in traces_d), default=0)
    t_s = np.arange(n_frames) * float(ex_time)
    if traces_d and include_all:
        written['traces_all'] = save_traces_csv(
            traces_d, traces_a, t_s, path('traces_all'), verbose=False)

    keep = [r for r in records if r.keep]
    good_ids = sorted({int(r.emitter_id) for r in keep})
    if good_ids:
        written['traces_good'] = save_traces_csv(
            [traces_d[i] for i in good_ids], [traces_a[i] for i in good_ids],
            t_s, path('traces_good'), verbose=False)

    # ---- 2) 逐分子 E 表 ----
    centres = {int(t.emitter_id): (float(t.x), float(t.y)) for t in traces_d}
    written['fret_csv'] = save_csv(results or [], path('fret_csv'), cal,
                                   centres=centres, verbose=False)

    # ---- 3) 合格 trace 总览 ----
    if keep:
        fig = plot_pair_gallery(keep, traces_a, traces_d, which='kept',
                                n=len(keep), ncols=ncols, ex_time=ex_time,
                                pre_win=pre_win, post_win=post_win)
        if title:
            fig.suptitle(f'{title}   kept traces: {len(keep)}', fontsize=11,
                         y=0.998)
            # plot_pair_gallery 已经 tight_layout 过, 加 suptitle 要重新留出边
            fig.tight_layout(rect=(0, 0, 1, 1.0 - 0.35 / max(fig.get_figheight(), 1e-6)))
        fig.savefig(path('gallery_png'), dpi=dpi, bbox_inches='tight',
                    facecolor='white')
        plt.close(fig)
        written['gallery_png'] = path('gallery_png')

    # ---- 4) E 分布 + 两个交叉检查 ----
    E = np.array([r.E for r in (results or [])], float)
    if np.isfinite(E).any():
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.0))
        plot_efficiency(results, ax=axes[0])
        plot_efficiency_vs_brightness(results, ax=axes[1], key='IA_pre')
        plot_efficiency_vs_brightness(results, ax=axes[2], key='ID_pre')
        g = cal.gamma if cal is not None else np.nan
        a = cal.alpha if cal is not None else np.nan
        b = cal.beta if cal is not None else np.nan
        fig.suptitle(f'{title}   smFRET  (gamma={g:.2f}, alpha={a:.3f}, '
                     f'beta={b:.0f}; n={int(np.isfinite(E).sum())})',
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(path('e_png'), dpi=200, bbox_inches='tight',
                    facecolor='white')
        plt.close(fig)
        written['e_png'] = path('e_png')

    if verbose:
        print(f'[export] {out_dir}')
        if 'traces_all' in written:
            print(f'  {FILES["traces_all"]:<20} {len(traces_d)} 分子 x {n_frames} 帧')
        if 'traces_good' in written:
            print(f'  {FILES["traces_good"]:<20} {len(good_ids)} 分子 x {n_frames} 帧')
        print(f'  {FILES["fret_csv"]:<20} {len(results or [])} 分子 (含 γ/α/β)')
        if 'gallery_png' in written:
            print(f'  {FILES["gallery_png"]:<20} {len(keep)} 条合格 trace')
        if 'e_png' in written:
            ok = E[np.isfinite(E)]
            print(f'  {FILES["e_png"]:<20} ' +
                  _E_MED.format(med=float(np.median(ok)), n=ok.size))
        if not keep:
            print('  [!] 没有合格 trace: 只写了空表')
    return written
