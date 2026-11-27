# smlm/multiplex/viz.py
"""配准相关的可视化: 位移场、残差分布、残差空间分布、两通道叠加检查。

只用 matplotlib, 不依赖 notebook 的全局状态 —— 传 PairSet / 标定对象进来即可。
需要中文字体时自行设置 ``plt.rcParams``, 这里统一用英文标签。
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ..datatypes import MultiImageStack
from .drift import drift_side_needs
from .fit import design, poly_terms


def _multichannel_mean_montage(multi, frame=None):
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    panels = []
    for channel_index in range(multi.n_channels):
        data = np.asarray(multi.data[channel_index])
        if frame is None:
            # 累加用 float32, 避免把整个 (T,H,W) 升成 float64 再求均值
            image = data.mean(axis=0, dtype=np.float32).astype(float)
        else:
            image = np.asarray(data[int(frame)], dtype=float)
        lo, hi = np.percentile(image, (0.2, 99.8))
        scale = max(float(hi - lo), np.finfo(float).eps)
        panels.append(np.clip((image - lo) / scale, 0.0, 1.0))
    _, _, height, width = multi.shape
    return np.concatenate(panels, axis=1), height, width


def _label_montage_channels(ax, multi, width, height):
    for channel_index, channel in enumerate(multi.channels):
        x_offset = channel_index * width
        ax.text(x_offset + width / 2, -8, channel, ha='center',
                va='bottom', color='cyan', clip_on=False)
        if channel_index:
            ax.axvline(x_offset - 0.5, color='white', linewidth=0.8,
                       alpha=0.8)
    ax.set_xlim(-0.5, width * multi.n_channels - 0.5)
    ax.set_ylim(height - 0.5, -12)
    ax.set_axis_off()


def plot_multichannel_crop(multi, y, x, figsize=None, frame=None):
    """Show one shared local crop window across adjacent channel images.

    ``frame=None`` averages over time (default); pass an integer (e.g. ``0``)
    to display that single frame, which is the drift reference frame.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    y, x = tuple(y), tuple(x)
    preview = multi.crop(y=y, x=x)
    y0, y1 = (int(value) for value in y)
    x0, x1 = (int(value) for value in x)
    canvas, height, width = _multichannel_mean_montage(multi, frame=frame)
    _, _, crop_height, crop_width = preview.shape

    if figsize is None:
        figsize = (max(6.0, 5.0 * multi.n_channels), 5.0)
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(canvas, cmap='gray', origin='upper', interpolation='nearest',
              vmin=0.0, vmax=1.0)
    for channel_index, _ in enumerate(multi.channels):
        x_offset = channel_index * width
        ax.add_patch(Rectangle(
            (x_offset + x0 - 0.5, y0 - 0.5),
            x1 - x0, y1 - y0,
            fill=False, edgecolor='yellow', linewidth=1.5))
    _label_montage_channels(ax, multi, width, height)
    source = 'mean' if frame is None else f'frame {int(frame)}'
    ax.set_title(
        f'crop window ({source}): y={y0}:{y1}, x={x0}:{x1} '
        f'({crop_height}x{crop_width})')
    fig.tight_layout()
    return fig


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


def plot_multichannel_roi(
        frame,
        windows: Mapping[str, tuple[int, int, int, int]],
        figsize: tuple[float, float] = (10, 6)):
    """Show multichannel ROI windows on a source image or image stack."""
    img = np.asarray(getattr(frame, 'data', frame), float)
    if img.ndim == 3:
        img = img.mean(axis=0)
    if img.ndim != 2:
        raise ValueError(f'frame must be (T,H,W) or (H,W), got {img.shape}')

    h, w = img.shape
    fig, ax = plt.subplots(figsize=figsize)
    vmin, vmax = np.percentile(img, (0.2, 99.8))
    ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax, origin='upper')
    cmap = plt.get_cmap('tab10')
    for i, (channel, (y0, x0, height, width)) in enumerate(windows.items()):
        color = cmap(i % cmap.N)
        ax.add_patch(plt.Rectangle(
            (x0 - 0.5, y0 - 0.5), width, height,
            fill=False, ec=color, lw=1.6))
        ax.text(x0 + 2, y0 + 10, f'{channel}  x{x0}..{x0 + width}  '
                                f'y{y0}..{y0 + height}',
                color=color, fontsize=8, va='top',
                bbox=dict(fc='black', alpha=0.45, lw=0, pad=2))
    ax.set_xlim(-0.5, w - 0.5)
    ax.set_ylim(h - 0.5, -0.5)
    ax.set_title('multichannel ROI windows')
    ax.set_xlabel('detector column')
    ax.set_ylabel('detector row')
    fig.tight_layout()
    return fig


def _ordered_channel_locs(multi, locs):
    """Validate a ``{channel: Localizations}`` mapping and order it."""
    if not isinstance(locs, Mapping):
        raise TypeError('locs must be a {channel: Localizations} mapping, '
                        f'got {type(locs).__name__}')
    missing = [channel for channel in multi.channels if channel not in locs]
    if missing:
        raise ValueError(f'missing localization sets for channels {missing}')
    ordered = [locs[channel] for channel in multi.channels]
    if any(current is None for current in ordered):
        raise ValueError('localization sets must not be None')
    return ordered


def plot_multichannel_localizations(multi, locs, figsize=None):
    """Plot channel localizations on adjacent channel images in local coordinates.

    ``locs`` maps each channel name to its localizations.
    """
    locs = _ordered_channel_locs(multi, locs)

    if figsize is None:
        figsize = (max(6.0, 5.0 * multi.n_channels), 5.0)
    canvas, height, width = _multichannel_mean_montage(multi)
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(canvas, cmap='gray', origin='upper', interpolation='nearest',
              vmin=0.0, vmax=1.0)
    for index, (channel, channel_locs) in enumerate(
            zip(multi.channels, locs)):
        ax.scatter(channel_locs.x + index * width, channel_locs.y,
                   s=12, marker='x', c='lime', linewidths=0.8)
        ax.text(index * width + width / 2, -8,
                f'{channel}: {len(channel_locs)} locs', ha='center',
                va='bottom', color='cyan', clip_on=False)
        if index:
            ax.axvline(index * width - 0.5, color='white',
                       linewidth=0.8, alpha=0.8)
    ax.set_xlim(-0.5, width * multi.n_channels - 0.5)
    ax.set_ylim(height - 0.5, -12)
    ax.set_axis_off()
    fig.tight_layout()
    return fig


def plot_multi_drift(drift, ex_time=None, unit='pixel', pixel_size=None,
                     figsize=None):
    """Plot every channel's drift on one shared time axis.

    Two stacked panels (dx and dy, both in each channel's own pixels) share the
    horizontal time axis, so channel-to-channel differences are directly
    comparable. ``ex_time`` switches the axis from frames to seconds.
    """
    if not isinstance(drift, Mapping):
        raise TypeError('drift must be a {channel: DriftTrace} mapping')
    if not drift:
        raise ValueError('drift mapping is empty')
    frames = {len(trace.dx) for trace in drift.values()}
    if len(frames) != 1:
        raise ValueError(f'channels have different frame counts {sorted(frames)};'
                         ' a shared time axis needs one frame count')
    n_frames = frames.pop()

    from ..viz import _unit_scale
    scale, axis_unit = _unit_scale(unit, pixel_size)
    if ex_time is not None:
        time = np.arange(n_frames, dtype=float) * float(ex_time)
        xlabel = 'time (s)'
    else:
        time = np.arange(n_frames, dtype=float)
        xlabel = 'frame'

    if figsize is None:
        figsize = (7.0, 5.0)
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
    for channel, trace in drift.items():
        axes[0].plot(time, np.asarray(trace.dx, float) * scale,
                     lw=1.0, label=channel)
        axes[1].plot(time, np.asarray(trace.dy, float) * scale,
                     lw=1.0, label=channel)
    axes[0].set_ylabel(f'dx ({axis_unit})')
    axes[1].set_ylabel(f'dy ({axis_unit})')
    axes[1].set_xlabel(xlabel)
    axes[0].set_title('drift in each channel\'s own pixels', fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8, loc='best')
    fig.tight_layout()
    return fig, axes


def plot_drift_fov(multi, drift, *, channels=None, frame=0, figsize=None,
                   transpose=False):
    """Mark the drift-limited common area on each channel and label its margins.

    No per-frame boxes: only the region shared by every frame (the largest
    usable area) is outlined with a thin border, and the minimum pixels lost on
    each side are printed. The background is one channel frame, which is the
    drift reference frame by default.

    ``transpose=True`` rotates the images into landscape orientation; margin
    labels stay attached to their original side.
    """
    if not isinstance(multi, MultiImageStack):
        raise TypeError('multi must be a MultiImageStack')
    if not isinstance(drift, Mapping):
        raise TypeError('drift must be a {channel: DriftTrace} mapping')
    names = list(multi.channels) if channels is None else [
        channel for channel in channels if channel in drift]
    if not names:
        raise ValueError('no matching channels to plot')

    _, height, width = multi.shape[-3:]
    disp_width, disp_height = (height, width) if transpose else (width, height)
    to_display = ({'top': 'left', 'bottom': 'right', 'left': 'top',
                   'right': 'bottom'} if transpose
                  else {side: side for side in ('top', 'bottom', 'left',
                                                'right')})
    positions = {
        'top': ((disp_width / 2.0, -9.0), 0),
        'bottom': ((disp_width / 2.0, disp_height + 13.0), 0),
        'left': ((-10.0, disp_height / 2.0), 90),
        'right': ((disp_width + 10.0, disp_height / 2.0), 90),
    }

    if figsize is None:
        figsize = (5.6 * len(names), 4.4)
    fig, axes = plt.subplots(1, len(names), figsize=figsize, squeeze=False)
    axes = axes.ravel()
    for ax, channel in zip(axes, names):
        data = np.asarray(multi.data[multi.index(channel)])
        if frame is None:
            image = data.mean(axis=0, dtype=np.float32).astype(float)
        else:
            image = np.asarray(data[int(frame)], dtype=float)
        if transpose:
            image = image.T
        lo, hi = np.percentile(image, (0.2, 99.8))
        scale = max(float(hi - lo), np.finfo(float).eps)
        ax.imshow(np.clip((image - lo) / scale, 0.0, 1.0), cmap='gray',
                  origin='upper', interpolation='nearest', vmin=0.0, vmax=1.0)

        trace = drift[channel]
        dx = np.asarray(trace.dx, float)
        dy = np.asarray(trace.dy, float)
        need = drift_side_needs(trace)
        x_lo, x_hi = -float(dx.min()), width - float(dx.max())
        y_lo, y_hi = -float(dy.min()), height - float(dy.max())
        rect = ((y_lo, x_lo, y_hi - y_lo, x_hi - x_lo) if transpose
                else (x_lo, y_lo, x_hi - x_lo, y_hi - y_lo))

        ax.add_patch(Rectangle(rect[:2], rect[2], rect[3], fill=False,
                               edgecolor='lime', linewidth=0.8,
                               linestyle='--'))

        for side in ('top', 'bottom', 'left', 'right'):
            (x, y), rotation = positions[to_display[side]]
            ax.text(x, y, f'{side} {need[side]} px', color='white',
                    fontsize=8, ha='center', va='center', rotation=rotation,
                    clip_on=False,
                    bbox=dict(boxstyle='round,pad=0.2', fc='black',
                              ec='none', alpha=0.55))

        ax.set_xlim(-22, disp_width + 22)
        ax.set_ylim(disp_height + 26, -22)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        source = 'mean' if frame is None else f'frame {int(frame)}'
        ax.set_title(f'{channel} ({source}): total drift '
                     f'{trace.total_drift:.2f} px  |  common area '
                     f'{int(x_hi - x_lo)}x{int(y_hi - y_lo)}', fontsize=9)
    fig.tight_layout()
    return fig, axes


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


# ---------------------------------------------------------------------------
#  多通道 trace: 每个子图叠加所选通道, 颜色 = 通道, 亮度 = 字段
# ---------------------------------------------------------------------------
_CHANNEL_COLORS = {'g': 'tab:green', 'green': 'tab:green', 'donor': 'tab:green',
                   'r': 'tab:red', 'red': 'tab:red', 'acceptor': 'tab:red'}
_SPARE_COLORS = ('tab:blue', 'tab:orange', 'tab:purple', 'tab:brown', '0.35')
_FIELD_ALPHAS = (1.0, 0.45, 0.30, 0.22, 0.18)


def _channel_colors(labels, colors=None):
    """通道 -> 颜色: 显式给 > 按名字猜(g 绿 / r 红) > 备用色。"""
    given = dict(colors or {})
    unknown = [k for k in given if k not in labels]
    if unknown:
        raise ValueError(f'colors 里有未选中的通道 {unknown}; 已选: {list(labels)}')
    spare = iter([c for c in _SPARE_COLORS if c not in given.values()])
    out = {}
    for label in labels:
        if label in given:
            out[label] = given[label]
            continue
        out[label] = _CHANNEL_COLORS.get(str(label).lower(),
                                         next(spare, '0.3'))
    return out


def _scalars_per_channel(scalars, labels):
    """``scalars`` 可以是一个序列(所有通道共用) 或 ``{通道: 序列}``。"""
    if scalars is None:
        return {label: () for label in labels}
    if isinstance(scalars, Mapping):
        unknown = [k for k in scalars if k not in labels]
        if unknown:
            raise ValueError(f'scalars 里有未选中的通道 {unknown}; 已选: {list(labels)}')
        out = {}
        for label in labels:
            value = scalars.get(label, ())
            out[label] = (value,) if isinstance(value, str) else tuple(value)
        return out
    names = (scalars,) if isinstance(scalars, str) else tuple(scalars)
    return {label: names for label in labels}


def plot_channel_traces(traces, *, channels=None, colors=None, groups=None,
                        indices=None, n_max=6, x=None, xlabel='Frame',
                        alphas=None, linewidth=0.8, scalars=None, bands=None,
                        figsize=None, title=None, yscale=None, legend=True,
                        ylim_bottom=0.0):
    """同一批分子的多通道 trace: **每个子图叠加所选通道**。

    与 :func:`smlm.trace.viz.plot_trace_fields` 的分工:

    * ``plot_trace_fields`` —— 一条 trace 的多个字段分列 (列 = 字段组);
    * 这里 —— 列 = ``groups``, 每个子图里叠加**所有选中通道**。
      颜色专门用来区分通道, 组内不同字段改用**亮度**(alpha) 与线宽区分,
      所以 ``fields=('raw', 'bg')`` 画出来是"同一颜色的亮/暗两条线"。

    Parameters
    ----------
    traces : dict
        ``{通道: [Trace]}``, 通常就是 ``extract_multichannel_traces`` 的返回值。
        各通道列表长度必须一致且按 ``emitter_id`` 对齐。
    channels : sequence, optional
        要叠加的通道 (默认全部); 顺序决定叠放顺序 (后画的在上层)。
    colors : dict, optional
        ``{通道: 颜色}``; 默认 ``g``/``donor`` 绿、``r``/``acceptor`` 红。
    groups : dict | sequence, optional
        与 :func:`plot_trace_fields` 完全同构 —— 每组一个子图列::

            groups={'signal': dict(fields=('intensity',)),
                    'raw + bg': dict(fields=('raw', 'bg'))}

        组内第一个字段最亮, 后面的依次变暗 (见 ``alphas``)。
    indices : sequence of int, optional
        画哪些分子; 默认前 ``n_max`` 个。
    alphas : sequence of float, optional
        组内字段的亮度梯度, 默认 ``(1.0, 0.45, 0.30, ...)``。
    scalars : str | sequence | dict, optional
        单个数字字段, 写进图例 (每个通道各一份, 前缀通道名);
        也可以给 ``{通道: 序列}`` 只给部分通道。
    bands : str, optional
        状态字段名; 只在**第一个**通道上画 (两条色带叠起来会糊)。
    x : array-like, optional
        横轴; 默认用第一个通道的 ``frames``。
    ylim_bottom : float | None
        纵轴下限 (默认 ``0``): 上限自适应, 下限统一压到 0, 各行零点对齐才比得出
        差异。``None`` = 完全自适应。

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> mp.plot_channel_traces(
    ...     multiTraces, channels=('g', 'r'),
    ...     groups={'signal': dict(fields=('intensity',)),
    ...             'raw + bg': dict(fields=('raw', 'bg'))},
    ...     scalars={'g': ('radius',), 'r': ('radius',)})
    """
    from matplotlib.lines import Line2D

    from ..trace.viz import (_clamp_bottom, _draw_bands, _field_label,
                             _format_scalar, _normalize_groups, _one_trace)

    if not isinstance(traces, Mapping) or not traces:
        raise TypeError('traces 必须是 {通道: [Trace]} 的映射')
    labels = list(traces) if channels is None else list(channels)
    if not labels:
        raise ValueError('channels 是空的')
    unknown = [c for c in labels if c not in traces]
    if unknown:
        raise ValueError(f'channels 里有未知通道 {unknown}; '
                         f'可用: {list(traces)}')

    lists = {label: _one_trace(traces[label]) for label in labels}
    n = len(lists[labels[0]])
    for label in labels[1:]:
        if len(lists[label]) != n:
            raise ValueError(f'通道 {label!r} 有 {len(lists[label])} 条 trace, '
                             f'{labels[0]!r} 有 {n} 条 —— 索引必须对齐')
    if n == 0:
        raise ValueError('traces 里没有 trace')

    indices = (list(range(min(n, int(n_max)))) if indices is None
               else list(indices))
    if not indices:
        raise ValueError('indices 是空的')
    for idx in indices:
        if not 0 <= idx < n:
            raise IndexError(f'index {idx} 超出 0..{n - 1}')

    specs = _normalize_groups(groups if groups is not None else ('intensity',))
    alphas = tuple(alphas) if alphas is not None else _FIELD_ALPHAS
    if not alphas:
        raise ValueError('alphas 不能为空')
    color_of = _channel_colors(labels, colors)
    scalar_names = _scalars_per_channel(scalars, labels)

    nrows, ncols = len(indices), len(specs)
    if figsize is None:
        figsize = (3.4 * ncols + 0.7, 1.6 * nrows + 0.9)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False,
                             sharex='col')

    for r, idx in enumerate(indices):
        xv = (np.asarray(lists[labels[0]][idx]['frames']) if x is None
              else np.asarray(x))
        for c, group in enumerate(specs):
            ax = axes[r][c]
            handles, texts = [], []
            if bands is not None:
                primary = lists[labels[0]][idx]
                if bands not in primary:
                    raise KeyError(f'通道 {labels[0]!r} 的 trace #{idx} '
                                   f'没有字段 {bands!r}')
                _draw_bands(ax, xv, np.asarray(primary[bands]),
                            bleach_frame=int(primary.get('bleach_frame', -1)))
            for label in labels:
                trace = lists[label][idx]
                for k, name in enumerate(group['fields']):
                    if name not in trace:
                        raise KeyError(f'通道 {label!r} 的 trace #{idx} 没有字段 '
                                       f'{name!r}; 可用: {sorted(trace)}')
                    alpha = float(alphas[k % len(alphas)])
                    lw = linewidth * (1.0 - 0.15 * min(k, 3))
                    ax.plot(xv, np.asarray(trace[name], float),
                            color=color_of[label], alpha=alpha, lw=lw,
                            linestyle=group['linestyle'] or '-')
                    handles.append(Line2D([], [], color=color_of[label],
                                          alpha=alpha, lw=lw))
                    texts.append(f'{label} · {_field_label(name)}')
            for label in labels:
                trace = lists[label][idx]
                for name in scalar_names[label]:
                    if name in trace:
                        texts.append(f'{label} {_format_scalar(name, trace[name])}')
                        handles.append(Line2D([], [], linestyle='none', marker=''))
            if legend and texts and (len(labels) > 1 or len(group['fields']) > 1
                                     or scalar_names[labels[0]]):
                ax.legend(handles, texts, fontsize=6.0, frameon=False,
                          loc='best')
            ax.set_title(f"#{lists[labels[0]][idx].get('emitter_id', idx)}  "
                         f"{group['name']}", fontsize=9)
            ax.tick_params(labelsize=7)
            if len(group['fields']) == 1:
                ax.set_ylabel(_field_label(group['fields'][0]), fontsize=8)
            if yscale:
                ax.set_yscale(yscale)
            if ylim_bottom is not None and not yscale:
                _clamp_bottom(ax, ylim_bottom)
            if r == nrows - 1:
                ax.set_xlabel(xlabel, fontsize=8)

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
