# smlm/fret/viz.py
"""FRET 诊断图: trace 对画廊、γ 标定点、E 分布。

一律返回 matplotlib ``Figure``, 由调用方决定怎么显示/保存。横轴统一用秒
(``ex_time``), 因为"多久"比"第几帧"更容易读。
"""
from __future__ import annotations

import numpy as np

from .io import TracesByID


# ---------------------------------------------------------------------------
#  图像块
# ---------------------------------------------------------------------------
def crop_centered(vol, cy, cx, half):
    """(T,H,W) -> (T,2h+1,2h+1), 越界补 0。"""
    vol = np.asarray(vol)
    F, H, W = vol.shape
    n = 2 * int(half) + 1
    r0, c0 = int(round(cy)) - int(half), int(round(cx)) - int(half)
    out = np.zeros((F, n, n), vol.dtype if vol.dtype.kind == 'f' else np.float32)
    rr0, rr1 = max(0, r0), min(H, r0 + n)
    cc0, cc1 = max(0, c0), min(W, c0 + n)
    if rr1 > rr0 and cc1 > cc0:
        out[:, rr0 - r0:rr1 - r0, cc0 - c0:cc1 - c0] = vol[:, rr0:rr1, cc0:cc1]
    return out


# ---------------------------------------------------------------------------
#  trace
# ---------------------------------------------------------------------------
def plot_pair(ax, trace_a, trace_d, ev=None, *, ex_time=1.0, mle_a=None,
              mle_d=None, pre_win=100, post_win=400, title=None,
              x_max=None, axhline_zero=True):
    """在一张 axes 上画双通道 trace (红=受体/绿=供体), x 轴是**秒**。

    实线 = 孔径积分(已扣背景); 点线 = 从孔径里扣掉的背景;
    虚线 = 逐帧 MLE 光子数 (不受环形背景影响, 用来判断背景是否可信);
    灰底 = 算 E 用的 pre / post 窗口, 黑虚线 = 台阶帧。

    ``ev`` 为 None 时只画曲线 (例如台阶没找到的否决样例)。
    """
    ia = np.asarray(trace_a['intensity'], float)
    idd = np.asarray(trace_d['intensity'], float)
    n = min(ia.size, idd.size)
    ia, idd = ia[:n], idd[:n]
    t = np.arange(n) * float(ex_time)

    if ev is not None:
        j, lag = int(ev.frame), int(ev.lag)
        ax.axvspan(t[max(0, j - pre_win)], t[min(n - 1, j)], color='0.90', lw=0)
        ax.axvspan(t[min(n - 1, j + lag)],
                   t[min(n - 1, j + lag + post_win)], color='0.80', lw=0)
        ax.axvline(t[min(j, n - 1)], color='k', ls='--', lw=0.8)
    if axhline_zero:
        ax.axhline(0.0, color='0.75', lw=0.6)

    ax.plot(t, idd, color='green', lw=0.8, label='donor aperture')
    bg_d = trace_d.get('bg_flux')
    if bg_d is not None:
        ax.plot(t, np.asarray(bg_d)[:n], color='green', lw=0.9, ls=':',
                alpha=0.9, label='donor bg (subtracted)')
    if mle_d is not None:
        ax.plot(t, np.asarray(mle_d)[:n], color='green', lw=0.5, ls='--',
                alpha=0.75, label='donor MLE')

    ax.plot(t, ia, color='red', lw=0.8, alpha=0.9, label='acceptor aperture')
    bg_a = trace_a.get('bg_flux')
    if bg_a is not None:
        ax.plot(t, np.asarray(bg_a)[:n], color='red', lw=0.9, ls=':',
                alpha=0.9, label='acceptor bg (subtracted)')
    if mle_a is not None:
        ax.plot(t, np.asarray(mle_a)[:n], color='red', lw=0.5, ls='--',
                alpha=0.75, label='acceptor MLE')

    ax.set_xlabel('time (s)', fontsize=7)
    ax.set_ylabel('photons / frame', fontsize=7)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=5.5, loc='upper right', frameon=False, ncol=2)
    if x_max is not None:
        ax.set_xlim(0.0, float(x_max))
    if title is None and ev is not None:
        title = (f'#{ev.emitter_id}   A {ev.IA_pre:.0f}→{ev.IA_post:.0f}   '
                 f'D {ev.ID_pre:.0f}→{ev.ID_post:.0f}   lag={ev.lag}')
    if title:
        ax.set_title(title, fontsize=7)
    return ax


# ---------------------------------------------------------------------------
#  汇总图(静态, 给 notebook 用)
# ---------------------------------------------------------------------------
def _pair_title(rec, ex_time=None):
    """画廊标题: 编号 / keep|原因 (＋具体哪条物理门限) / E / 关键诊断。"""
    tag = 'KEEP' if rec.keep else 'REJECT: ' + ', '.join(rec.reasons)
    if not rec.keep and rec.failed_flags:
        tag += ' (' + ','.join(rec.failed_flags) + ')'
    e = f'  E={rec.E:.2f}' if np.isfinite(rec.E) else ''
    ev = rec.event
    if ev is None:
        return f'#{rec.emitter_id}  {tag}{e}\nI_A<={rec.IA_max:.0f}  (no step)'
    return (f'#{rec.emitter_id}  {tag}{e}\n'
            f'I_A {ev.IA_pre:.0f}→{ev.IA_post:.0f}  '
            f'I_D {ev.ID_pre:.0f}→{ev.ID_post:.0f}  lag={ev.lag}\n'
            f'spread={rec.spread:.1f}σ  n_d^pre={rec.n_d_pre}  '
            f'neigh={rec.neighbours}')


def plot_pair_gallery(records, traces_a, traces_d, *, which='both', n=4,
                      ex_time=1.0, pre_win=100, post_win=400, ncols=None,
                      figsize=None, dpi=110):
    """通过 / 否决样例网格, 用来目视校验阈值和找"漏杀/误杀"。

    Parameters
    ----------
    records : list[PairRecord]
        :func:`~smlm.fret.pairs.select_pairs` 的输出。想只看某个子集 (例如
        ``E < 0.85`` 的那些) 就先把列表筛好再传进来, 配合 ``which='kept'``。
    traces_a, traces_d : list or dict-like
        trace 序列或 ``{emitter_id: trace}``。内部**按 ``emitter_id`` 取**,
        所以只载入一部分分子时 (下标 ≠ 编号) 也不会取错。
    which : 'kept' | 'rejected' | 'both'
    n : int
        每组最多画几条。否决组按"首个否决原因"排序后均匀取样, 所以能看到各种
        原因, 不是随机挑。
    ex_time : float
        每帧秒数 (x 轴用秒)。
    ncols : int, optional
        每行几个面板; 超过就折行。默认每组一行 (条数少时最紧凑)。

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    groups = []
    if which in ('kept', 'both'):
        pool = [r for r in records if r.keep]
        groups.append((f'kept ({len(pool)})', pool))
    if which in ('rejected', 'both'):
        from .pairs import FUNNEL
        order = {k: i for i, (k, _) in enumerate(FUNNEL)}
        pool = [r for r in records if not r.keep]
        pool.sort(key=lambda r: (order.get(r.reasons[0], 99) if r.reasons else 99,
                                 r.emitter_id))
        groups.append((f'rejected ({len(pool)})', pool))

    picked = []
    for label, pool in groups:
        if len(pool) > n:
            sel = np.unique(np.linspace(0, len(pool) - 1, n).round().astype(int))
            pool = [pool[i] for i in sel]
        picked.append((label, pool))

    ncol = int(ncols) if ncols else max([len(p) for _, p in picked] + [1])
    blocks = []                       # (label, 该行 records, 是否打印组名)
    for label, pool in picked:
        if not pool:
            blocks.append((label, [], True))
            continue
        for k in range(0, len(pool), ncol):
            blocks.append((label, pool[k:k + ncol], k == 0))

    figsize = figsize or (3.1 * ncol, 1.9 * max(len(blocks), 1))
    fig, axes = plt.subplots(len(blocks), ncol, figsize=figsize, dpi=dpi,
                             squeeze=False)
    tr_a, tr_d = TracesByID(traces_a), TracesByID(traces_d)
    for r, (label, chunk, first) in enumerate(blocks):
        for c in range(ncol):
            ax = axes[r][c]
            if c >= len(chunk):
                ax.axis('off')
                continue
            rec = chunk[c]
            plot_pair(ax, tr_a[rec.emitter_id], tr_d[rec.emitter_id],
                      rec.event, ex_time=ex_time, pre_win=pre_win,
                      post_win=post_win, title=_pair_title(rec))
            for side in ('top', 'right'):
                ax.spines[side].set_visible(False)
            if c == 0 and first:
                ax.text(-0.28, 0.5, label, transform=ax.transAxes, rotation=90,
                        va='center', ha='center', fontsize=8, color='0.3')
    fig.tight_layout()
    return fig


def plot_efficiency(results, ax=None, *, bins=30, show_don=True, title=None):
    """E 分布直方图。

    灰阶梯是 ``E_don = 1 − I_D^pre/I_D^post``: 它完全不看受体通道, 所以不依赖
    γ/α/β —— 两条分布差很多就说明 γ 或串色修正有问题。

    ``E ≳ 1`` 的堆积要单独看: 那些分子的 ``I_D^pre`` 与 0 不可区分, E 实际上是
    **下界**而不是测量值 (分母被 S_A 主导), 图上会把比例标出来。
    """
    import matplotlib.pyplot as plt

    E = np.array([r.E for r in results], float)
    E = E[np.isfinite(E)]
    if E.size == 0:
        raise ValueError('没有有限的 E 可画')
    lo, hi = min(0.0, float(E.min())), max(1.0, float(E.max()))
    edges = np.linspace(lo, hi, int(bins) + 1)
    ax = ax or plt.gca()
    ax.hist(E, bins=edges, color='tab:purple', alpha=0.7,
            label=f'E (n={E.size})')
    if show_don:
        Ed = np.array([r.E_don for r in results], float)
        Ed = Ed[np.isfinite(Ed) & (Ed >= lo) & (Ed <= hi)]
        if Ed.size:
            ax.hist(Ed, bins=edges, histtype='step', color='0.35', lw=1.0,
                    label='E_don (γ-free)')

    med = float(np.median(E))
    p16, p84 = np.percentile(E, [16, 84])
    ax.axvline(med, color='k', ls='--', lw=1.0)
    for v in (p16, p84):
        ax.axvline(v, color='k', ls=':', lw=0.9)
    n_hi = int((E > 0.95).sum())
    ax.text(0.02, 0.97,
            f'median {med:.3f}\np16~p84 {p16:.3f}~{p84:.3f}\n'
            f'E>0.95: {n_hi} ({n_hi / E.size:.0%})  lower bound\n'
            f'E<0: {int((E < 0).sum())}',
            transform=ax.transAxes, va='top', ha='left', fontsize=7)
    ax.set_xlabel('FRET efficiency E')
    ax.set_ylabel('count')
    ax.set_title(title or 'E distribution', fontsize=9)
    ax.legend(fontsize=7, loc='center right', frameon=False)
    return ax


def plot_efficiency_vs_brightness(results, ax=None, *, key='IA_pre', nbins=5):
    """E vs 亮度: 检查暗分子是否给出系统性偏差。

    ``key`` 可以是 :class:`~smlm.fret.datatypes.FretResult` 的任一字段, 常用
    ``IA_pre`` (受体亮度) / ``ID_post`` (去淬灭后供体) / ``snr_A``。
    黑点 = 按亮度分位数分箱后的中位 ± p16~p84。
    """
    import matplotlib.pyplot as plt

    x = np.array([getattr(r, key) for r in results], float)
    y = np.array([r.E for r in results], float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() == 0:
        raise ValueError(f'{key} / E 没有有限值可画')
    x, y = x[ok], y[ok]
    ax = ax or plt.gca()
    ax.scatter(x, y, s=14, alpha=0.55, c='tab:purple', edgecolors='none')

    edges = np.unique(np.quantile(x, np.linspace(0.0, 1.0, int(nbins) + 1)))
    for a, b in zip(edges[:-1], edges[1:]):
        m = (x >= a) & (x <= b)
        if m.sum() < 3:
            continue
        c = float(np.median(x[m]))
        lo, mid, hi = np.percentile(y[m], [16, 50, 84])
        ax.errorbar([c], [mid], yerr=[[mid - lo], [hi - mid]], fmt='o-',
                    color='0.15', lw=1.0, ms=3.5, capsize=0, zorder=3)
    if (x > 0).all():
        ax.set_xscale('log')
    ax.set_xlabel(key.replace('_', ' '))
    ax.set_ylabel('E')
    ax.set_title(f'E vs {key}', fontsize=9)
    ax.grid(alpha=0.25)
    return ax


def plot_gamma(events, cal, ax=None, raw=True):
    """ΔI_A vs ΔI_D 散点 + γ 三条估计线。"""
    import matplotlib.pyplot as plt
    dD = np.array([e.dID for e in events], float)
    dA = np.array([e.dIA for e in events], float)
    ok = np.isfinite(dD) & np.isfinite(dA)
    dD, dA = dD[ok], dA[ok]
    ax = ax or plt.gca()
    ax.scatter(dD, dA, s=16, alpha=0.75, c='purple')
    if dD.size:
        xs = np.linspace(0, float(np.nanmax(dD)) * 1.05, 50)
        for g, c, ls, lab in ((cal.gamma, 'r', '-', 'γ_med'),
                              (cal.gamma_rma, 'b', '-.', 'γ_rma'),
                              (cal.gamma_ols, 'g', '--', 'γ_ols')):
            if np.isfinite(g):
                ax.plot(xs, g * xs, color=c, ls=ls, lw=1.1,
                        label=f'{lab}={g:.2f}')
        lim = [0, max(float(np.nanmax(dD)), float(np.nanmax(dA))) * 1.05]
        ax.plot(lim, lim, 'k--', lw=0.8, label='γ=1')
    ax.set_xlabel('ΔI_D (photons)')
    ax.set_ylabel('ΔI_A (photons)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    return ax
