"""two_emitter.viz —— 双发射体数据集的显示辅助。"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch

_CJK_FONTS = ['Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC',
              'Source Han Sans SC', 'Arial Unicode MS', 'PingFang SC', 'SimSun']

STATE_COLORS = {0: '#2ca02c', 1: '#c8c8c8', 2: '#d62728'}
STATE_LABELS = {0: 'on', 1: 'dark', 2: 'bleached'}


def set_cjk_font():
    """把 matplotlib 默认字体换成可用的中文字体 (避免标题出现方块)。"""
    available = {f.name for f in font_manager.fontManager.ttflist}
    zh = next((f for f in _CJK_FONTS if f in available), None)
    if zh:
        plt.rcParams['font.sans-serif'] = [zh, 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    return zh


# --------------------------------------------------------------------------- #
# 内部工具
# --------------------------------------------------------------------------- #
def _data(ds, photons=True):
    return ds.photons() if photons else ds.images.astype(np.float32)


def _frame_image(ds, frames=None, photons=True):
    data = _data(ds, photons)
    if frames is None:
        return data.mean(axis=0)
    frames = np.atleast_1d(np.asarray(frames, dtype=int))
    return data[frames].mean(axis=0)


def _extent_nm(shape, pixel_size, origin=(0, 0)):
    """返回 ``imshow(origin='upper')`` 用的 extent (nm): ``(left, right, bottom, top)``。

    与全仓库统一约定的坐标约定一致 (见 ``dataset.roi`` / ``level0`` /
    ``level1.JointHMMTwoEmitter``): **像素索引 i 覆盖 [i·px, (i+1)·px],
    中心在 (i+0.5)·px**。所以 ``origin=(r0, c0)`` 的 ROI 的 extent 是
    ``(c0·px, (c0+w)·px, (r0+h)·px, r0·px)`` —— 注意 ``bottom > top``,
    这正是 ``origin='upper'`` 想要的 (第 0 行画在上面, y 向下增大),
    与 ``imshow`` 的默认 extent ``(-0.5, w-0.5, h-0.5, -0.5)`` 同构。

    早期版本这里误用了 ``(i-0.5)·px``, 会让整幅图相对真值位置**平移半个像素**
    (106.7 nm 像素下就是 53 nm), 表现为 GT 标记偏离 PSF 中心。
    """
    h, w = shape
    r0, c0 = origin
    return (c0 * pixel_size, (c0 + w) * pixel_size,
            (r0 + h) * pixel_size, r0 * pixel_size)


def _imshow(ax, img, pixel_size, origin=(0, 0), cmap='gray',
            percentile=(0.5, 99.8), photons=True, colorbar=True):
    vmin, vmax = np.percentile(img, percentile) if percentile else (None, None)
    im = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=_extent_nm(img.shape, pixel_size, origin),
                   origin='upper', interpolation='nearest')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                           label='photons' if photons else 'ADU')
    return im


# --------------------------------------------------------------------------- #
# 图像
# --------------------------------------------------------------------------- #
def show_frame(ds, frame: int = 0, ax=None, photons: bool = True, cmap='gray',
               percentile=(0.5, 99.8), title=None, colorbar=True):
    """显示单帧 (默认转成光子数)。"""
    img = _data(ds, photons)[frame]
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    _imshow(ax, img, ds.pixel_size, cmap=cmap, percentile=percentile,
            photons=photons, colorbar=colorbar)
    ax.set_title(title if title is not None else
                 f"{ds.kind}: frame {frame}  (d={ds.separation_nm:.0f} nm)")
    return ax


def show_average(ds, frames=None, ax=None, photons: bool = True, cmap='gray',
                 percentile=(0.5, 99.8), title=None, colorbar=True):
    """显示所选帧的平均图。``frames`` 可为 int / 列表 / None(全部)。"""
    img = _frame_image(ds, frames=frames, photons=photons)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    _imshow(ax, img, ds.pixel_size, cmap=cmap, percentile=percentile,
            photons=photons, colorbar=colorbar)
    if title is None:
        n = ds.n_frames if frames is None else np.atleast_1d(np.asarray(frames)).size
        title = f"{ds.kind}: mean of {n} frames  (d={ds.separation_nm:.0f} nm)"
    ax.set_title(title)
    return ax


def show_pair_roi(ds, pair_idx: int, frames=None, ax=None, half_px: int = 6,
                  photons: bool = True, show_gt: bool = True, cmap='gray',
                  percentile=(0.5, 99.8), title=None, colorbar=True):
    """放大显示某个 pair 的 ROI(所选帧平均), 并叠加真值位置。"""
    ys, xs = ds.roi(pair_idx, half_px=half_px)
    img = _frame_image(ds, frames=frames, photons=photons)[ys, xs]
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    _imshow(ax, img, ds.pixel_size, origin=(ys.start, xs.start), cmap=cmap,
            percentile=percentile, photons=photons, colorbar=colorbar)

    if show_gt:
        idx = ds.pair_emitters(pair_idx)
        pts = ds.gt_xyz[idx, :2]
        ax.scatter(pts[:, 0], pts[:, 1], s=90, facecolors='none',
                   edgecolors='r', linewidths=1.5, label='GT')
        if pts.shape[0] == 2:
            ax.plot(pts[:, 0], pts[:, 1], 'r--', lw=1.0)
        ax.legend(loc='upper right', fontsize=8)

    if title is None:
        n = ds.n_frames if frames is None else np.atleast_1d(np.asarray(frames)).size
        title = (f"pair {pair_idx}  d={ds.pair_sep_nm[pair_idx]:.0f} nm  "
                 f"({n} frames averaged)")
    ax.set_title(title)
    return ax


# --------------------------------------------------------------------------- #
# 状态
# --------------------------------------------------------------------------- #
def plot_states(ds, pair_idx: int, ax=None, t_max: int | None = None):
    """画该 pair 两个发射体的 ON/OFF/bleached 状态时间线 (真值)。"""
    idx = ds.pair_emitters(pair_idx)
    n = ds.n_frames if t_max is None else min(ds.n_frames, int(t_max))
    S = ds.states[idx][:, :n]

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 0.9 + 0.55 * len(idx)))

    for k in range(len(idx)):
        for st, col in STATE_COLORS.items():
            m = S[k] == st
            if np.any(m):
                ax.fill_between(np.arange(n), k, k + 1, where=m, step='mid',
                                color=col, linewidth=0)
    ax.set_yticks(np.arange(len(idx)) + 0.5)
    ax.set_yticklabels([f"E{int(i)}" for i in idx])
    ax.set_xlim(0, n)
    ax.set_ylim(0, len(idx))
    ax.set_xlabel('frame')
    ax.set_title(f"pair {pair_idx}: ON/OFF 状态 (真值, d={ds.pair_sep_nm[pair_idx]:.0f} nm)")
    ax.legend(handles=[Patch(color=STATE_COLORS[s], label=STATE_LABELS[s])
                       for s in (0, 1, 2)],
              loc='upper right', fontsize=8, ncol=3)
    return ax


# --------------------------------------------------------------------------- #
# level0 汇总图
# --------------------------------------------------------------------------- #
def plot_separation_summary(rows, ax=None, d_min_nm: float | None = None,
                            yerr_key='sep_est_std_nm', logy=False):
    """把 :func:`two_emitter.level0.summarize_oracle` 的结果画成"估计间距 vs 真值。"""
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 5.5))

    sep = np.array([r["separation_nm"] for r in rows], dtype=float)
    est = np.array([r["sep_est_mean_nm"] for r in rows], dtype=float)
    std = np.array([r.get(yerr_key, np.nan) for r in rows], dtype=float)

    ax.errorbar(sep, est, yerr=std, marker='o', ms=4, lw=1.0, capsize=2.5,
                label='oracle 估计 (mean ± std)', color='#1f77b4')

    lo = float(np.nanmin([np.nanmin(sep), np.nanmin(est)])) if len(rows) else 0.0
    hi = float(np.nanmax([np.nanmax(sep), np.nanmax(est)])) if len(rows) else 1.0
    pad = 0.08 * max(hi - lo, 1.0)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k--', lw=1.0,
            label='y = x (理想)')

    if d_min_nm is not None and np.isfinite(d_min_nm):
        ax.axvline(d_min_nm, color='#d62728', ls=':', lw=1.4,
                   label=f'level0 d_min = {d_min_nm:.2f} nm')

    ax.set_xlabel('真值间距 d (nm)')
    ax.set_ylabel('oracle 估计间距 (nm)')
    ax.set_title('level0 (已知态): 间距恢复')
    if logy:
        ax.set_yscale('log')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    return ax
