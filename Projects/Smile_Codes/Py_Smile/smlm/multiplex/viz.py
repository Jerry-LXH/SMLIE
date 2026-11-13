# smlm/multiplex/viz.py
"""配准相关的可视化: 位移场、残差分布、残差空间分布、两通道叠加检查。

只用 matplotlib, 不依赖 notebook 的全局状态 —— 传 PairSet / 标定对象进来即可。
需要中文字体时自行设置 ``plt.rcParams``, 这里统一用英文标签。
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from .fit import design, poly_terms


def _residuals(calib, pairs, method='poly'):
    R = np.asarray(pairs.mov, float)
    G = np.asarray(pairs.ref, float)
    pred = calib.transform(R, method=method)
    d = pred - G
    return R, np.hypot(d[:, 0], d[:, 1])


def plot_residual_hist(calib, pairs, ax=None, method='poly', bins=40):
    """配准残差直方图 (poly vs affine)。"""
    _, rp = _residuals(calib, pairs, 'poly')
    _, ra = _residuals(calib, pairs, 'affine')
    ax = ax or plt.gca()
    hi = max(float(rp.max()), float(ra.max()), 1e-3)
    b = np.linspace(0, hi * 1.05, int(bins))
    ax.hist(ra, bins=b, alpha=0.55,
            label=f'affine  med {np.median(ra):.2f}  rms {np.sqrt((ra**2).mean()):.2f}')
    ax.hist(rp, bins=b, alpha=0.55,
            label=f'poly({calib.poly_order})  med {np.median(rp):.2f}  '
                  f'rms {np.sqrt((rp**2).mean()):.2f}')
    ax.set_xlabel('|residual| (px)')
    ax.set_ylabel('pairs')
    ax.set_title(f'registration residual (n={len(pairs)})', fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    return ax


def plot_residual_map(calib, pairs, ax=None, method='poly'):
    """残差的空间分布: 若还有系统性花纹 -> 阶数不够或视场边缘形变大。"""
    R, r = _residuals(calib, pairs, method)
    ax = ax or plt.gca()
    sc = ax.scatter(R[:, 0], R[:, 1], c=r, s=26, cmap='magma', vmin=0)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xlabel('mov x (px)')
    ax.set_ylabel('mov y (px)')
    ax.set_title(f'residual map ({method})', fontsize=9)
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label='|residual| (px)')
    return ax


def plot_displacement_field(calib, pairs=None, ax=None, grid=13, gain=25.0,
                            pad=12.0):
    """位移场: 扣掉整体平移后的部分 (即"纯平移修不掉"的那一部分)。

    注意 order 1 时这部分完全是线性的(缩放/切变), 不是"非线性形变";
    标题因此叫 displacement field 而不是 nonlinear field。
    """
    ax = ax or plt.gca()
    if pairs is not None and len(pairs):
        R = np.asarray(pairs.mov, float)
        x_lo, x_hi = R[:, 0].min() - pad, R[:, 0].max() + pad
        y_lo, y_hi = R[:, 1].min() - pad, R[:, 1].max() + pad
    else:
        cov = calib.coverage()
        x_lo, x_hi = cov.get('x', [0, 256])
        y_lo, y_hi = cov.get('y', [0, 256])
        x_lo, x_hi = x_lo - pad, x_hi + pad
        y_lo, y_hi = y_lo - pad, y_hi + pad

    gx, gy = np.meshgrid(np.linspace(x_lo, x_hi, int(grid)),
                         np.linspace(y_lo, y_hi, int(grid)))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    mapped = calib.apply_poly(pts)
    du, dv = (mapped - pts).T
    # 只留非线性部分: 减掉位移场自身的均值
    du_n, dv_n = du - du.mean(), dv - dv.mean()
    mag = np.hypot(du_n, dv_n)
    q = ax.quiver(pts[:, 0], pts[:, 1], du_n, dv_n, mag, cmap='viridis',
                  angles='xy', scale_units='xy', scale=1.0 / float(gain),
                  width=0.008)
    if pairs is not None and len(pairs):
        R = np.asarray(pairs.mov, float)
        ax.scatter(R[:, 0], R[:, 1], s=8, c='w', alpha=0.7, label='pairs')
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_hi, y_lo)
    ax.set_xlabel('mov x (px)')
    ax.set_ylabel('mov y (px)')
    ax.set_title(f'poly({calib.poly_order}) displacement field, mean removed '
                 f'(x{gain:.0f})\nmax |d - mean| = {mag.max():.2f} px', fontsize=9)
    plt.colorbar(q, ax=ax, fraction=0.046, pad=0.04, label='|d - mean| (px)')
    return ax


def plot_overview(calib, pairs=None, figsize=(14, 4.4)):
    """三联图: 位移场 / 残差直方图 / 残差空间分布。

    ``pairs=None`` 时只画位移场(标定文件里没存配对点, 常见)。
    """
    if pairs is None or len(pairs) == 0:
        fig, ax = plt.subplots(figsize=(4.6, 4.2))
        plot_displacement_field(calib, None, ax=ax)
        fig.suptitle(f'{calib.moving} -> {calib.reference}   '
                     f'order {calib.poly_order}   '
                     f'coarse={tuple(calib.coarse_offset)}', fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        return fig

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 3, wspace=0.34)
    plot_displacement_field(calib, pairs, ax=fig.add_subplot(gs[0, 0]))
    plot_residual_hist(calib, pairs, ax=fig.add_subplot(gs[0, 1]))
    plot_residual_map(calib, pairs, ax=fig.add_subplot(gs[0, 2]))
    fig.suptitle(f'{calib.moving} -> {calib.reference}   n={len(pairs)}   '
                 f'coarse={tuple(calib.coarse_offset)}', fontsize=10)
    return fig


def plot_channel_check(data_ref, data_mov, calib=None, frame=0,
                       zoom=45, figsize=(11, 5.2)):
    """两通道平均图 + 点叠加, 看是否真的重合。

    ``data_*`` 是 (T,H,W) 或 (H,W)。给了 calib 就把 mov 的点映射到 ref 再画。
    """
    img = np.asarray(data_ref, float)
    img = img.mean(axis=0) if img.ndim == 3 else img
    im_mov = np.asarray(data_mov, float)
    im_mov = im_mov.mean(axis=0) if im_mov.ndim == 3 else im_mov

    vmin, vmax = np.percentile(img, (0.2, 99.8))
    h, w = img.shape
    cy, cx = h // 2, w // 2          # 放大看视场中央

    fig, axs = plt.subplots(1, 2, figsize=figsize)
    for ax, second, title in [(axs[0], None, 'reference channel'),
                              (axs[1], im_mov, 'reference + moving (mapped)')]:
        ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, origin='upper')
        if second is not None:
            v2 = np.percentile(im_mov, (0.2, 99.8))
            ax.imshow(im_mov, cmap='magma', vmin=v2[0], vmax=v2[1],
                      origin='upper', alpha=0.45)
        ax.set_xlim(cx - zoom, cx + zoom)
        ax.set_ylim(cy + zoom, cy - zoom)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
#  取框预览: 让人眼决定两个通道切到哪里
# ---------------------------------------------------------------------------
def plot_cut_windows(frame, crops=None, *, split_x=None, seam_hint=None,
                     reference='g', moving='r', figsize=(14, 6.2),
                     profile=True):
    """在双视角平均图上画出两通道的**取框**, 供人眼确认后再切。

    为什么必须看一眼: 探测器上两条半幅交界处有一条鬼影/串扰亮带, 落在里面的
    "分子"可能是另一个通道的尾巴。但这条带有多宽、从哪一列开始, 跟光路和焦面
    都有关系 —— 不能拍一个常数。所以取框是人的判断, 本函数只负责把判断画出来。

    Parameters
    ----------
    frame : (T,H,W) 或 (H,W)
        整幅双视角图(未切通道)。
    crops : dict or None
        ``{'g': dict(y=(y0, y1), x=(x0, x1)), 'r': {...}}``。None 时只画
        ``split_x`` 分界线与鬼影带提示, 不画框(先看图再定框)。
    split_x : int
        分界线列号; None 时取图宽的一半。
    seam_hint : int
        鬼影带的**提示**宽度(px), 用灰色带标出 split_x ± seam_hint。

    Returns
    -------
    matplotlib Figure
    """
    img = np.asarray(getattr(frame, 'data', frame), float)
    img = img.mean(axis=0) if img.ndim == 3 else img
    h, w = img.shape
    if split_x is None:
        split_x = w // 2

    ncol = 2 if profile else 1
    fig, axs = plt.subplots(1, ncol, figsize=figsize,
                            gridspec_kw={'width_ratios': [1.25, 1]} if profile else None)
    ax = np.atleast_1d(axs)[0]

    vmin, vmax = np.percentile(img, (0.2, 99.8))
    ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, origin='upper')

    # 分界线 + 鬼影带
    ax.axvline(split_x - 0.5, color='yellow', ls='--', lw=1.1)
    if seam_hint:
        for x0 in (split_x - seam_hint, split_x + seam_hint):
            ax.axvline(x0 - 0.5, color='orange', ls=':', lw=1.0)
        ax.axvspan(split_x - seam_hint, split_x + seam_hint,
                   color='orange', alpha=0.16, lw=0)
        # 标在图底部, 免得和两个框的标签重叠
        ax.text(split_x, h * 0.985, f'ghost band\n±{seam_hint} px',
                color='orange', fontsize=7, ha='center', va='bottom')

    # 两通道取框
    colors = {reference: 'lime', moving: 'red'}
    if crops:
        for ch in (reference, moving):
            if ch not in crops:
                continue
            y0, y1 = (int(v) for v in crops[ch]['y'])
            x0, x1 = (int(v) for v in crops[ch]['x'])
            ax.add_patch(plt.Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0,
                                       fill=False, ec=colors[ch], lw=1.6))
            ax.text(x0 + 2, y0 + 10, f'{ch}  x{x0}..{x1}  y{y0}..{y1}\n'
                                      f'{x1 - x0} x {y1 - y0} px',
                    color=colors[ch], fontsize=7, va='top',
                    bbox=dict(fc='black', alpha=0.45, lw=0, pad=2))
            # 被让出去的部分涂暗, 一眼看出丢了什么
            ax.axvspan(x0, x1, ymin=1 - (y0 + (h - y1)) / h, ymax=1 - y0 / h,
                       color=colors[ch], alpha=0.07, lw=0)
    ax.set_title('dual-view frame — confirm the two cut boxes', fontsize=9)
    ax.set_xlabel('detector column')
    ax.set_ylabel('detector row')

    if profile:
        axp = np.atleast_1d(axs)[1]
        prof = img.mean(axis=0)
        axp.plot(np.arange(w), prof, lw=0.8, color='0.35')
        axp.axvline(split_x - 0.5, color='yellow', ls='--', lw=1.1)
        if seam_hint:
            axp.axvspan(split_x - seam_hint, split_x + seam_hint,
                        color='orange', alpha=0.16, lw=0)
        if crops:
            for ch in (reference, moving):
                if ch in crops:
                    x0, x1 = (int(v) for v in crops[ch]['x'])
                    axp.axvspan(x0, x1, color=colors[ch], alpha=0.08, lw=0)
        axp.set_xlabel('detector column')
        axp.set_ylabel('mean photons / px')
        axp.set_title('column profile — look for the bright ghost band', fontsize=9)
        axp.grid(alpha=0.3)
    fig.tight_layout()
    return fig
