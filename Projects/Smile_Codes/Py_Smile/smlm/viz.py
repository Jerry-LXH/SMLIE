# smlm/viz.py
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

def _resolve_frames(stack, frames, stride=1):
    """把 frames 参数规整成整数数组(相对 data 的帧下标)。None -> 全部帧。

    frames 可为 int / 列表 / range / slice / 布尔掩码; stride>1 时再等间隔抽帧。
    """
    n = int(stack.n_frames)
    if frames is None:
        idx = np.arange(n)
    elif isinstance(frames, slice):
        idx = np.arange(n)[frames]
    else:
        sel = np.asarray(frames)
        if sel.dtype == bool:
            if sel.size != n:
                raise ValueError(f"布尔掩码长度 {sel.size} 与帧数 {n} 不一致。")
            idx = np.flatnonzero(sel)
        else:
            idx = np.atleast_1d(sel).astype(int)
    if idx.size and (idx.min() < 0 or idx.max() >= n):
        raise IndexError(f"帧索引越界, 合法范围 [0, {n - 1}]")
    if int(stride) > 1:
        idx = idx[::int(stride)]
    if idx.size == 0:
        raise ValueError("frames / stride 选完没有任何帧。")
    return idx


def average_frames(stack, frames=None):
    """对所选帧取平均, 返回 (H, W) 图像。frames 可为 int / 列表 / range / None。"""
    frames = _resolve_frames(stack, frames)
    return stack.data[frames,:, :].mean(axis=0)


def _unit_scale(unit, pixel_size):
    """按 unit 返回 (scale, axis_unit)。scale 是像素坐标要乘的因子。

    unit: 'pixel' | 'nm' | 'um'。后两者需要 pixel_size(nm/pixel)。
    """
    if unit == 'pixel':
        return 1.0, 'px'
    if unit in ('nm', 'um'):
        if pixel_size is None:
            raise ValueError(f"unit={unit!r} 需要 pixel_size (nm/pixel)。")
        ps = float(pixel_size)
        if unit == 'nm':
            return ps, 'nm'
        return ps / 1000.0, 'µm'      # nm/pixel -> µm/pixel
    raise ValueError(f"unknown unit={unit!r}, 只支持 'pixel' / 'nm' / 'um'。")


def show_frames(stack, frames=None, ax=None, cmap='gray',
                percentile=(0.1, 99.9), title=None, colorbar=True,
                unit='pixel', pixel_size=None):
    """显示所选帧的平均图。

    percentile: 用分位数做对比度拉伸, SMLM 单分子图动态范围大, 这样比 min/max 好看。
    unit:       'pixel'(默认) 或 'nm'。'nm' 时需传 pixel_size(nm/pixel);
                仅改变坐标轴刻度(通过 extent), 不缩放图像数据。
    pixel_size: nm/pixel, unit='nm' 时必填。
    """
    img = average_frames(stack, frames)
    scale, axis_unit = _unit_scale(unit, pixel_size)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    vmin, vmax = np.percentile(img, percentile) if percentile else (None, None)

    # 像素中心 i 对应 [i-0.5, i+0.5]; 乘 scale 换到目标单位
    nrow, ncol = img.shape
    extent = (-0.5 * scale, (ncol - 0.5) * scale,
              (nrow - 0.5) * scale, -0.5 * scale)   # 上下翻转, 配 origin='upper'

    im = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax,
                   origin='upper', extent=extent)

    if title is None:
        fr = _resolve_frames(stack, frames)
        title = f"avg of {fr.size} frame(s)" if fr.size > 1 else f"frame {fr[0]}"
    ax.set_title(title)
    ax.set_xlabel(f"x ({axis_unit})"); ax.set_ylabel(f"y ({axis_unit})")
    if colorbar:
        ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="photons")
    return ax


def overlay_points(stack, points, frames=None, ax=None,
                   marker='o', color='red', size=40, facecolor='none',
                   linewidth=1.2, unit='pixel', pixel_size=None, **imshow_kw):
    """在平均图上叠加 Candidates 或 Localizations 的坐标。

    points: Candidates(整数像素) 或 Localizations(亚像素) 均可, 只要有 .frame/.x/.y。
    frames: 选定要显示的帧; 只画落在这些帧内的点。None = 全部。
    unit:   'pixel'(默认) 或 'nm'; 与底图保持一致, 点坐标随之换算。
    pixel_size: nm/pixel, unit='nm' 时必填。
    """
    frames = _resolve_frames(stack, frames)
    scale, _ = _unit_scale(unit, pixel_size)

    # unit / pixel_size 透传给 show_frames, 保证底图与点用同一坐标系
    ax = show_frames(stack, frames=frames, ax=ax,
                     unit=unit, pixel_size=pixel_size, **imshow_kw)

    # 只保留落在所选帧里的点, 并换算到目标单位
    sel = np.isin(points.frame, frames)
    x = np.asarray(points.x)[sel] * scale
    y = np.asarray(points.y)[sel] * scale

    label = f"{type(points).__name__} (n={x.size})"
    unfilled = marker in ('+', 'x', '|', '_', '.', ',', '1', '2', '3', '4')
    if unfilled:
        ax.scatter(x, y, s=size, marker=marker, c=color,
                   linewidths=linewidth, label=label)
    else:
        ax.scatter(x, y, s=size, marker=marker, facecolors=facecolor,
                   edgecolors=color, linewidths=linewidth, label=label)

    ax.legend(loc='upper right', fontsize=8, framealpha=0.6)
    # xlim/ylim 也按单位缩放, 否则 nm 模式下点会跑出可视范围
    ncol, nrow = stack.data.shape[2], stack.data.shape[1]
    ax.set_xlim(-0.5 * scale, (ncol - 0.5) * scale)
    ax.set_ylim((nrow - 0.5) * scale, -0.5 * scale)   # y 翻转, 与 imshow 一致
    return ax



def plot_drift(trace, ex_time=None, ax=None, title='Tracking',
               cmap='jet_r', linewidth=2, unit='pixel', show=True,
               pixel_size=None):
    """可视化 XY 漂移轨迹, 按时间着色。

    对应 MATLAB plotTracking。
      x 轴 = dx (列/col 方向), y 轴 = dy (行/row 方向)。
      ex_time    : 每帧曝光(秒)。给了则 colorbar 显示时间(s), 否则显示帧号。
      unit       : 'pixel'(默认) 或 'nm'。
      pixel_size : nm/pixel。unit='nm' 时必填; 若为 None 则退回 trace.pixel_size。
      show       : True 时函数结束自动 plt.show()。
    """
    x = np.asarray(trace.dx, dtype=float)
    y = np.asarray(trace.dy, dtype=float)
    F = len(x)

    # 显式 pixel_size 优先, 否则用 trace 自带的
    ps = pixel_size if pixel_size is not None else trace.pixel_size
    scale, axis_unit = _unit_scale(unit, ps)
    x = x * scale
    y = y * scale

    # --- 着色量: 有曝光时间就用秒, 否则用帧号 ---
    if ex_time is not None:
        t = np.arange(F) * ex_time
        cbar_label = 'Time (s)'
    else:
        t = np.arange(F, dtype=float)
        cbar_label = 'Frame'

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    # --- 构造分段, 每段颜色 = 该段起点的时间 (interp 着色) ---
    pts = np.column_stack([x, y]).reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)   # (F-1, 2, 2)

    lc = LineCollection(segs, cmap=cmap, linewidth=linewidth)
    lc.set_array(t[:-1])                 # 每段用起点时间上色
    lc.set_clim(t.min(), t.max())        # 对应 caxis([min max])
    line = ax.add_collection(lc)

    # LineCollection 不会自动撑开坐标轴, 手动设范围 + 留边
    pad_x = 0.05 * (x.max() - x.min() + 1e-9)
    pad_y = 0.05 * (y.max() - y.min() + 1e-9)
    ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
    ax.set_ylim(y.min() - pad_y, y.max() + pad_y)
    ax.set_aspect('equal')               # axis equal

    cb = fig.colorbar(line, ax=ax)
    cb.set_label(cbar_label)

    ax.set_xlabel(f'x/col ({axis_unit})')
    ax.set_ylabel(f'y/row ({axis_unit})')
    ax.set_title(title)
    ax.tick_params(labelsize=11)

    if show:
        plt.show()

    return fig, ax


def _plot_field_hist(ax, data, label, unit, *, bins='auto',
                     clip_pct=(0.5, 99.5), color='steelblue'):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        ax.text(0.5, 0.5, f"{label}: no data", ha='center', va='center',
                transform=ax.transAxes)
        ax.axis('off')
        return

    if clip_pct is not None:
        lo, hi = np.percentile(data, clip_pct)
    else:
        lo, hi = data.min(), data.max()
    if hi <= lo:
        hi = lo + 1.0
    shown = data[(data >= lo) & (data <= hi)]

    edges = _fd_edges(shown, lo, hi) if bins == 'auto' \
        else np.linspace(lo, hi, bins + 1)

    n_clipped = data.size - shown.size
    ax.hist(shown, bins=edges, color=color, edgecolor='none')
    ax.set_xlim(lo, hi)
    ax.set_xlabel(f"{label} ({unit})" if unit else label)
    ax.set_ylabel("count")
    ax.set_title(label)

    txt = f"median={np.median(data):.3g}\nN={data.size}"
    if n_clipped:
        txt += f"\n+{n_clipped} clipped"
    ax.text(0.97, 0.97, txt, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.8))

def _fd_edges(data, lo, hi, max_bins=120, min_bins=10):
    """Freedman–Diaconis 规则算 bin 数, 限定在 [lo, hi] 内。"""
    if data.size < 2:
        return np.linspace(lo, hi, min_bins + 1)
    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25
    if iqr <= 0:
        nb = int(np.sqrt(data.size))     # 退化时用 sqrt 规则
    else:
        width = 2 * iqr / np.cbrt(data.size)
        nb = int(np.ceil((hi - lo) / width)) if width > 0 else min_bins
    nb = int(np.clip(nb, min_bins, max_bins))
    return np.linspace(lo, hi, nb + 1)

_FIELDS = {
    'photons': ('Photons', 'photons/molecule'),
    'bg':      ('Background', 'photons/pixel'),
    'sigma':   ('PSF width', 'pixel'),
    'lp':      ('Localization precision', 'nm'),
}

def plot_distributions(locs, fields=('photons', 'bg', 'sigma', 'lp'),
                       bins='auto', figsize=(10, 8),
                       clip_pct=(0.5, 99.5)):
    """绘制各物理量分布直方图, 每张图独立自适应范围与 bins。"""
    fields = [f for f in fields if getattr(locs, f, None) is not None]
    n = len(fields)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()

    for ax, f in zip(axes, fields):
        label, unit = _FIELDS[f]
        _plot_field_hist(ax, getattr(locs, f), label, unit,
                         bins=bins, clip_pct=clip_pct)

    for ax in axes[n:]:
        ax.axis('off')

    fig.tight_layout()
    return fig, axes

def _track_xy(em: Emitter):
    """从单个 Emitter 取出逐帧 (x, y) 坐标。字段名不同就改这里。"""
    return np.asarray(em.x, dtype=float), np.asarray(em.y, dtype=float)

_PROJECTIONS = {'mean': np.mean, 'max': np.max,
               'std': np.std, 'median': np.median}


def _make_background(data, projection='mean', frame_idx=None):
    """把 data 归约成一张 2D 底图。
    data       : 3D 影片 (T,H,W) / 2D 图 / None
    frame_idx  : 选帧方式, 三者之一:
                 None             -> 用 projection 汇总全部帧;
                 int              -> 直接取该帧原图, 忽略 projection;
                 slice / 索引数组  -> 先取出这些帧得到子影片, 再按 projection 汇总。
                                    例如 slice(0, 200) / np.arange(0, 200) /
                                    [0, 50, 99]; 只选中 1 帧时等价于取原图。
    projection : 'mean' | 'max' | 'std' | 'median'  (需要汇总多帧时生效)
    """
    if data is None:
        return None
    arr = np.asarray(data)

    if arr.ndim == 2:                      # 已经是单张图
        return arr.astype(float)

    if arr.ndim == 3:                      # 影片 → 压成一张
        if frame_idx is not None:
            if (isinstance(frame_idx, (int, np.integer))
                    and not isinstance(frame_idx, bool)):
                return arr[int(frame_idx)].astype(float)      # 单帧原图
            # 注意: tuple 会被 numpy 当成多维下标 (arr[(1, 2)] == arr[1, 2]),
            # 所以除 slice 外一律转成 ndarray 再做花式索引。
            sel = frame_idx if isinstance(frame_idx, slice) \
                else np.asarray(frame_idx)
            sub = arr[sel]
            if sub.ndim == 2:                                 # 只选中 1 帧
                return sub.astype(float)
            if sub.shape[0] == 0:
                raise ValueError(
                    f"frame_idx={frame_idx!r} 没有选中任何帧 "
                    f"(影片共 {arr.shape[0]} 帧)。")
            return _PROJECTIONS[projection](sub, axis=0).astype(float)
        return _PROJECTIONS[projection](arr, axis=0).astype(float)

    raise ValueError(f"data 维度不支持: ndim={arr.ndim}")

def plot_trajectories(data, emitters, *, ax=None,
                      projection='mean', frame_idx=None, bg_cmap='gray',
                      color_by='time', cmap='viridis', linewidth=0.8,
                      mark_start=False, min_len=1, figsize=(8, 8),
                      interpolation='nearest', dpi=150, px_per_pixel=None,
                      zoom=None, show_grid=False):
    """把所有 emitter 的运动轨迹叠画在原始影像上。

    data        : 原始影片 (T,H,W) 或单张图 (H,W) 或 None。
    projection  : data 为影片时的投影方式 'mean'/'max'/'std'/'median'。
    frame_idx   : 底图选帧方式(优先于整段投影):
                  int            -> 直接用该帧原图作底;
                  slice / 数组    -> 先取出这些帧, 再按 projection 汇总, 例如
                                    slice(0, 200) / np.arange(0, 200) /
                                    [0, 50, 99];
                  None           -> 全部帧按 projection 汇总。
                  注意 (100, 200) 意为"第 100 和第 200 两帧", 不是区间;
                  要区间请用 slice(100, 201) 或 np.arange(100, 201)。
    color_by    : 'time'(每条轨迹沿时间渐变上色) / 'index' / 'none'。
    """
    # --- 先做底图, 才能据其尺寸决定 figsize (原代码在此处误用了未定义的 bg) ---
    bg = _make_background(data, projection=projection, frame_idx=frame_idx)

    if ax is None:
        if px_per_pixel is not None and bg is not None:
            h, w = bg.shape
            figsize = (w * px_per_pixel, h * px_per_pixel)
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure

    # --- 画底图 ---
    if bg is not None:
        ax.imshow(bg, cmap=bg_cmap, origin='upper',
                  interpolation=interpolation,
                  vmin=np.percentile(bg, 1), vmax=np.percentile(bg, 99.5))

    # --- 收集线段 ---
    # 对 'time': 每条轨迹按帧号逐段拆开, 每段颜色 = 该段起点帧
    # 对 'index': 每条轨迹整体一个颜色 (emitter 序号)
    # 对 'none': 单色
    segments, cvals, starts = [], [], []
    all_pts = []                     # 用于无底图时确定坐标范围

    for i, em in enumerate(emitters):
        x, y = _track_xy(em)
        if x.size < max(min_len, 1):
            continue

        pts = np.column_stack([x, y])
        all_pts.append(pts)
        starts.append(pts[0])

        if pts.shape[0] < 2:
            continue

        if color_by == 'time':
            # 逐段拆分, 每段用起点帧号上色
            frames = np.asarray(em.frame, dtype=float)
            segs = np.stack([pts[:-1], pts[1:]], axis=1)   # (n-1, 2, 2)
            segments.extend(segs)
            cvals.extend(frames[:-1])                       # 每段一个颜色值
        elif color_by == 'index':
            # 整条轨迹作为一条 polyline, 用 emitter 序号上色
            segments.append(pts)
            cvals.append(float(i))
        else:  # 'none'
            segments.append(pts)

    # --- 画线 ---
    if segments:
        if color_by in ('time', 'index'):
            cvals = np.asarray(cvals, dtype=float)
            lc = LineCollection(segments, cmap=cmap, linewidths=linewidth)
            lc.set_array(cvals)
            lc.set_clim(cvals.min(), cvals.max())   # 所有轨迹共用同一时间标尺
            ax.add_collection(lc)
            cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label('frame' if color_by == 'time' else 'emitter index')
        else:
            lc = LineCollection(segments, colors='tab:orange',
                                linewidths=linewidth)
            ax.add_collection(lc)

    if mark_start and starts:
        s = np.array(starts)
        ax.scatter(s[:, 0], s[:, 1], s=4, c='red', zorder=3)

    # --- 坐标范围 ---
    if bg is not None:                     # 有底图: 锁定到图像像素范围
        h, w = bg.shape
        ax.set_xlim(0, w)
        ax.set_ylim(h, 0)                  # y 向下, 与 imshow 一致
    else:                                  # 无底图: 按轨迹自动
        allpts = np.vstack(all_pts) if all_pts else np.empty((0, 2))
        if allpts.size:
            ax.set_xlim(allpts[:, 0].min(), allpts[:, 0].max())
            ax.set_ylim(allpts[:, 1].max(), allpts[:, 1].min())

    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x (pixel)')
    ax.set_ylabel('y (pixel)')
    ax.set_title(f"Trajectories (N={len(emitters)})")
    fig.tight_layout()
    return fig, ax

_CENTER_FUNCS = {'median': np.median, 'mean': np.mean,
                 'first': lambda a: a[0], 'last': lambda a: a[-1]}


def emitter_roi(stack, emitter, half=10, frames=None, stride=1, center='median'):
    """取以 emitter 为中心的"单分子窗口"影片, 返回 ImageStack。

    结果可直接喂给 ``smlm.io.save_movie`` 导出 mp4。返回的 metadata 里带了
    ``crop_offset``(y0, x0) 和 ``source_frames``(原始帧号), 所以 save_movie 会
    自动把 overlay 的定位点平移到窗口坐标, 并按原始帧号对齐。

    Parameters
    ----------
    half : int
        半窗宽, 窗口尺寸是 (2*half+1)^2 像素。
    frames : None | int | slice | 序列 | 布尔掩码
        帧选择, 语义同 save_movie; None = 全部帧。
    stride : int
        在 frames 之上再等间隔抽帧。
    center : 'median' | 'mean' | 'first' | 'last'
        窗口中心取 emitter 轨迹的哪个位置, 默认 'median'(抗离群点)。
    """
    from dataclasses import replace
    from .imgprocess import Crop

    if center not in _CENTER_FUNCS:
        raise ValueError(f"center 只支持 {list(_CENTER_FUNCS)}, 收到 {center!r}。")
    nf, h, w = stack.data.shape
    side = 2 * int(half) + 1
    if side > min(h, w):
        raise ValueError(f"窗口 {side}x{side} 比图像 {h}x{w} 还大, 请调小 half。")

    f = _CENTER_FUNCS[center]
    cx = int(round(float(f(np.asarray(emitter.x, dtype=float)))))
    cy = int(round(float(f(np.asarray(emitter.y, dtype=float)))))
    x0 = max(0, min(cx - int(half), w - side))     # 贴边时整体平移, 不缩窗
    y0 = max(0, min(cy - int(half), h - side))

    roi = Crop(y=(y0, y0 + side), x=(x0, x0 + side))(stack)
    idx = _resolve_frames(stack, frames, stride)
    meta = dict(roi.metadata)
    meta['source_frames'] = idx                    # 原始帧号, 给 save_movie 对齐用
    return replace(roi, data=roi.data[idx], metadata=meta)


def _as_index_list(idx, n):
    """把 idx 规整成整数列表。None -> 全部; 支持 int / slice / 序列 / 布尔掩码。"""
    if idx is None:
        return list(range(n))
    if isinstance(idx, (int, np.integer)) and not isinstance(idx, bool):
        idx = [int(idx)]
    elif isinstance(idx, slice):
        idx = list(range(n))[idx]
    else:
        sel = np.atleast_1d(np.asarray(idx))
        if sel.dtype == bool:
            if sel.size != n:
                raise ValueError(f"布尔掩码长度 {sel.size} 与 emitter 数 {n} 不一致。")
            sel = np.flatnonzero(sel)
        idx = [int(v) for v in sel]
    bad = [i for i in idx if i < 0 or i >= n]
    if bad:
        raise IndexError(f"emitter 序号越界 {bad}, 合法范围 [0, {n - 1}]。")
    return idx


def show_emitter_movies(stack, emitters, idx=None, half=10, frames=None,
                        stride=None, fps=20, overlay=True, ncols=6, size=1.6,
                        cmap='gray', percentile=(0.1, 99.9), clim=None,
                        max_frames=100, dpi=80, marker_kw=None):
    """内联播放若干 emitter 的"单分子窗口"小视频, 用来人工挑选。

    每个 emitter 一个窗口子图, 全部同步播放。返回 IPython.display.Image
    (动图), 把它放在 cell 最后一行即可在 notebook 里显示。
    挑好之后用 ``emitter_roi`` + ``smlm.io.save_movie`` 导出成 mp4。

    Parameters
    ----------
    idx : None | int | slice | 序列
        要看哪些 emitter。None = 全部(emitters 可能很多, 建议显式给一段)。
    half, frames : 见 ``emitter_roi``。
    stride : int, optional
        在 frames 上再抽帧。默认 None -> 自动抽到不超过 max_frames 帧,
        免得内联动图过大。
    fps : float
        notebook 里的预览播放帧率(只影响预览, 与导出无关)。
    overlay : bool
        是否叠加该 emitter 的定位点。
    ncols : int
        每行排几个窗口。
    max_frames : int
        stride=None 时自动抽帧的上限。
    """
    import io as _io
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from PIL import Image as _PILImage
    from IPython.display import Image as _IPyImage

    idxs = _as_index_list(idx, len(emitters))
    if not idxs:
        raise ValueError("idx 没有选中任何 emitter。")

    if stride is None:                       # 自动抽帧, 控制动图体积
        n_all = _resolve_frames(stack, frames).size
        stride = max(1, int(np.ceil(n_all / max(1, int(max_frames)))))
    rois = [emitter_roi(stack, emitters[i], half=half, frames=frames,
                        stride=stride) for i in idxs]

    n, nf = len(rois), rois[0].data.shape[0]
    side = rois[0].data.shape[1]
    if clim is None:
        pick = slice(None, None, max(1, nf // 5))      # 抽几帧估灰度范围即可
        sample = np.concatenate([r.data[pick].ravel() for r in rois])
        vmin, vmax = (float(v) for v in np.percentile(sample, percentile))
        if vmax <= vmin:
            vmax = vmin + 1.0
    else:
        vmin, vmax = float(clim[0]), float(clim[1])

    kw = dict(c='r', marker='x', s=14, linewidths=0.5)
    kw.update(marker_kw or {})

    nrows = int(np.ceil(n / ncols))
    fig = Figure(figsize=(size * min(n, ncols), size * nrows), dpi=dpi)
    canvas = FigureCanvasAgg(fig)

    ims, scs, pts_by_frame = [], [], []
    for k, roi in enumerate(rois):
        ax = fig.add_subplot(nrows, ncols, k + 1)
        ax.set_title(f"#{idxs[k]}", fontsize=7, pad=2)
        ax.set_xticks([]); ax.set_yticks([])
        ims.append(ax.imshow(np.zeros((side, side)), cmap=cmap, vmin=vmin,
                             vmax=vmax, origin='upper',
                             interpolation='nearest'))
        scs.append(ax.scatter([], [], **kw) if overlay else None)
        ax.set_xlim(-0.5, side - 0.5)          # 建好 scatter 后再锁范围
        ax.set_ylim(side - 0.5, -0.5)

        bucket = {}
        if overlay:
            em = emitters[idxs[k]]
            oy, ox = roi.metadata['crop_offset']
            xs = np.asarray(em.x, dtype=float) - ox
            ys = np.asarray(em.y, dtype=float) - oy
            pos = {int(f): j for j, f in enumerate(roi.metadata['source_frames'])}
            for t, f in enumerate(np.asarray(em.frame).astype(int)):
                j = pos.get(int(f))
                if j is not None:
                    bucket.setdefault(j, []).append((xs[t], ys[t]))
        pts_by_frame.append(bucket)

    def _grab(j):
        for k, roi in enumerate(rois):
            ims[k].set_data(roi.data[j])
            if scs[k] is not None:
                pts = pts_by_frame[k].get(j)
                scs[k].set_offsets(np.empty((0, 2)) if pts is None
                                   else np.asarray(pts, dtype=float))
        canvas.draw()
        return np.asarray(canvas.buffer_rgba())[..., :3].copy()

    pil = [_PILImage.fromarray(_grab(j)) for j in range(nf)]
    buf = _io.BytesIO()
    pil[0].save(buf, format='GIF', save_all=True, append_images=pil[1:],
                duration=max(20, int(round(1000.0 / max(float(fps), 1e-6)))),
                loop=0, optimize=True)
    print(f"[viz] 预览 {n} 个 emitter x {nf} 帧 (stride={stride}), "
          f"窗口 {side}x{side}, GIF {len(buf.getvalue()) / 1e6:.2f} MB")
    return _IPyImage(data=buf.getvalue(), format='gif')


_STAT_FIELDS = {
    'survival_sec':    ('Survival time',      's'),
    'n_locs':          ('Locs per emitter',   'count'),
    'duty_cycle':      ('Duty cycle',         ''),
    'eccentricity':    ('Eccentricity',       ''),
    'photons_per_sec': ('Photon rate',        'photons/s'),
    'photons_sum':     ('Total photons',      'photons'),
    'sigma_mean':      ('Mean PSF width',     'pixel'),
    'bg_mean':         ('Mean background',    'photons/pixel'),
}


def plot_emitter_stats(stats, fields=('survival_sec', 'n_locs', 'duty_cycle',
                                      'photons_sum', 'eccentricity', 'sigma_mean'),
                       bins='auto', clip_pct=(0.5, 99.5), figsize=(11, 7)):
    """把 EmitterStats 里的各标量字段画成分布直方图。"""
    fields = [f for f in fields
              if getattr(stats, f, None) is not None
              and np.asarray(getattr(stats, f)).ndim == 1]
    n = len(fields)
    ncol = 3 if n > 4 else 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()

    for ax, f in zip(axes, fields):
        label, unit = _STAT_FIELDS.get(f, (f, ''))
        _plot_field_hist(ax, getattr(stats, f), label, unit,
                         bins=bins, clip_pct=clip_pct)

    for ax in axes[n:]:
        ax.axis('off')
    fig.tight_layout()
    return fig, axes

def plot_bleaching_summary(stats, figsize=(4.5, 4)):
    """漂白 vs 未漂白(截断)的 emitter 数量对比。"""
    fig, ax = plt.subplots(figsize=figsize)
    b = np.asarray(stats.bleached, dtype=bool)
    n_bl, n_tr = int(b.sum()), int((~b).sum())
    ax.bar(['bleached', 'truncated'], [n_bl, n_tr],
           color=['tab:blue', 'tab:gray'])
    total = n_bl + n_tr
    for i, v in enumerate((n_bl, n_tr)):
        ax.text(i, v, f"{v}\n({v/total:.0%})", ha='center', va='bottom')
    ax.set_ylabel('emitter count')
    ax.set_title(f"Bleaching (N={total})")
    fig.tight_layout()
    return fig, ax


_PANEL_ORDER = ('intensity', 'bg', 'snr', 'blinking')


def _resolve_panels(panels, has_bg, has_snr):
    """把 panels 参数规整成有序的面板名元组(并按数据可用性过滤)。

    'all'/None 属于"自动模式": 缺数据的面板静默跳过, 不打扰用户;
    显式列出面板时, 若某个面板缺数据则给出提示。
    """
    if panels is None:
        names, explicit = list(_PANEL_ORDER), False
    elif isinstance(panels, str):
        key = panels.strip().lower()
        if key == 'all':
            names, explicit = list(_PANEL_ORDER), False
        elif key == 'trace':
            names, explicit = ['intensity'], True
        else:
            names, explicit = [key], True
    else:
        names, explicit = list(panels), True

    unknown = [n for n in names if n not in _PANEL_ORDER]
    if unknown:
        raise ValueError(
            f"未知面板 {unknown}; 可用: 'all' / 'trace' / "
            f"{list(_PANEL_ORDER)} 里的任意组合。")

    keep = [n for n in _PANEL_ORDER if n in names]
    dropped = [n for n in keep
               if (n == 'bg' and not has_bg) or (n == 'snr' and not has_snr)]
    if dropped and explicit:
        print(f"  [warn] trace 里没有 {'/'.join(dropped)} 数据, 已跳过该面板。")
    keep = [n for n in keep if n not in dropped]
    if not keep:
        raise ValueError("panels 解析后没有任何可绘制的面板。")
    return tuple(keep)


def check_trace(trace, emitter=None, ex_time=1.0, mle_brightness=None,
                use_time=True, show_track=True, frame_range=None,
                panels='all', show_raw=True, title=True, save_path=None):
    """可视化单个 Emitter:强度/背景/SNR/闪烁 + 空间轨迹,合并为一张图。

    Parameters
    ----------
    trace : Trace
        主数据对象。frames 为稠密帧轴(含 OFF 帧)。
    emitter : Emitter, optional
        用于 (1) 推断 ON/OFF 闪烁模式  (2) 绘制空间轨迹  (3) 映射 MLE 亮度。
    ex_time : float
        单帧曝光时间 (秒)。use_time=False 时忽略。
    mle_brightness : np.ndarray, optional
        MLE 估计的逐帧亮度(如 photons)。可为:
          - 稀疏(长度 == emitter.frame,仅 ON 帧):自动映射到稠密轴,OFF 填 NaN;
          - 稠密(长度 == trace.frames):直接使用。
        叠加在 raw / intensity 上作为第三条线。
    use_time : bool
        True -> 横轴为时间(s);False -> 横轴为帧号。
    show_track : bool
        是否绘制空间轨迹(需要 emitter)。
    frame_range : tuple[int | None, int | None], optional
        只展示指定帧号窗口(含两端),形如 (start, end)。
        start / end 是 trace.frames 中的**帧号**(不是数组下标);任一端可为
        None 表示该端不设限,例如 (1000, None) / (None, 2000) / (500, 800)。
        默认 None 表示展示全部帧。
        横轴沿用绝对帧号/时间(不重新归零);右侧轨迹图仍绘制完整轨迹。
    panels : 'all' | 'trace' | sequence of str
        要画哪些面板, 默认 'all' (= intensity + bg* + snr* + blinking,
        * 表示 trace 里没有该数据就自动跳过)。
        'trace' 等价于 ('intensity',), 只画强度曲线, 即"干净的 trace 图"。
        也可显式列出, 例如 ('intensity', 'blinking');
        绘制顺序固定为 intensity -> bg -> snr -> blinking。
    show_raw : bool
        是否在 intensity 面板里画 Raw 曲线, 默认 True。
        想得到最简洁的信号图可设 False, 例如:
        panels='trace', show_raw=False, show_track=False, title=False
        此时只剩一条背景扣减后的 intensity 曲线。
    title : bool | str
        True(默认) -> 自动标题 "Emitter #id (frames a-b, n=N)";
        False      -> 不画标题(配合 panels='trace' 得到最干净的图);
        str        -> 用该字符串作为标题。
    save_path : str, optional
        若给定,保存整张图为该路径(.png)。

    Returns
    -------
    fig : matplotlib Figure 句柄。
    """
    frames = np.asarray(trace.frames)
    raw = np.asarray(trace.raw, dtype=float)
    intensity = np.asarray(trace.intensity, dtype=float)

    # ---- 横轴 ----
    if use_time:
        x_axis = (frames - frames[0]) * ex_time
        x_label = "Time (s)"
    else:
        x_axis = frames.astype(float)
        x_label = "Frame"

    # ---- ON/OFF 闪烁模式 ----
    if emitter is not None:
        is_on = np.isin(frames, np.asarray(emitter.frame))
    else:
        is_on = ~np.isnan(intensity)

    # ---- MLE 亮度:稀疏 -> 稠密映射 ----
    mle = None
    if mle_brightness is not None:
        mle_arr = np.asarray(mle_brightness, dtype=float)
        if emitter is not None and mle_arr.size == np.asarray(emitter.frame).size:
            mle = np.full(frames.size, np.nan)
            mask = np.isin(frames, np.asarray(emitter.frame))
            mle[mask] = mle_arr
        elif mle_arr.size == frames.size:
            mle = mle_arr
        else:
            print(f"  [warn] mle_brightness 长度 {mle_arr.size} 既不匹配 "
                  f"emitter.frame 也不匹配 trace.frames,已忽略。")

    # ---- 可选通道(bg / snr)----
    # 先取出序列, 便于下面统一按帧窗口裁剪。
    bg = getattr(trace, "bg", None)
    bg = None if bg is None else np.asarray(bg, dtype=float)
    snr = getattr(trace, "snr", None)
    snr = None if snr is None else np.asarray(snr, dtype=float)

    # ---- 帧窗口裁剪(frame_range)----
    # 只裁"显示用"的序列; 横轴沿用绝对帧号/时间, 不重新归零, 便于与全量对照。
    if frame_range is not None:
        try:
            start, end = frame_range
        except (TypeError, ValueError):
            raise ValueError(
                "frame_range 需要能解包成两项的 (start, end), "
                "例如 (1000, 2000); 任一端可用 None 表示不设限。") from None
        win = np.ones(frames.size, dtype=bool)
        if start is not None:
            win &= frames >= start
        if end is not None:
            win &= frames <= end
        if not win.any():
            raise ValueError(
                f"frame_range={frame_range!r} 在 trace.frames "
                f"[{frames[0]}, {frames[-1]}] 内没有选中任何帧。")
        frames, raw, intensity = frames[win], raw[win], intensity[win]
        x_axis, is_on = x_axis[win], is_on[win]
        if mle is not None:
            mle = mle[win]
        if bg is not None:
            bg = bg[win]
        if snr is not None:
            snr = snr[win]

    # ---- 面板选择与布局 ----
    has_bg = bg is not None
    has_snr = snr is not None
    panel_names = _resolve_panels(panels, has_bg, has_snr)
    n_rows = len(panel_names)

    # 布局: 左列时间序列(n_rows 行), 右列轨迹(跨全高)
    draw_track = show_track and emitter is not None
    height = max(2.2 * n_rows, 2.8)
    if draw_track:
        fig = plt.figure(figsize=(12, height), constrained_layout=True)
        gs = GridSpec(n_rows, 2, figure=fig, width_ratios=[2, 1])
        axes = []
        for r in range(n_rows):
            ax = fig.add_subplot(gs[r, 0], sharex=axes[0] if axes else None)
            axes.append(ax)
        axt = fig.add_subplot(gs[:, 1])
    else:
        fig, axes = plt.subplots(n_rows, 1, figsize=(8, height),
                                 sharex=True, constrained_layout=True)
        axes = [axes] if n_rows == 1 else list(axes)
        axt = None

    if title is True:
        fig.suptitle(f"Emitter #{trace.emitter_id}  "
                     f"(frames {frames[0]}–{frames[-1]}, n={len(frames)})",
                     fontsize=13, fontweight="bold")
    elif title:
        fig.suptitle(str(title), fontsize=13, fontweight="bold")

    # ---- 逐面板绘制 ----
    for row, panel in enumerate(panel_names):
        ax = axes[row]

        if panel == 'intensity':            # raw / intensity / MLE
            if show_raw:
                ax.plot(x_axis, raw, "-", color="0.35", lw=1.0, label="Raw")
            ax.plot(x_axis, intensity, "-", color=(0, 0.447, 0.741), lw=1.5,
                    label="Intensity (bg-sub)")
            if mle is not None:
                ax.plot(x_axis, mle, "-", color=(0.85, 0.16, 0.55), lw=1.5,
                        marker=".", ms=4, label="MLE brightness")
            ax.axhline(0, ls="--", color="0.6", lw=1.0)
            ax.set_ylabel("Intensity")
            ax.legend(loc="best", fontsize=9)

        elif panel == 'bg':
            ax.plot(x_axis, bg, "-", color=(0.85, 0.33, 0.10), lw=1.2)
            ax.set_ylabel("Background\n(per px)")

        elif panel == 'snr':
            ax.plot(x_axis, snr, "-", color=(0.47, 0.67, 0.19), lw=1.2)
            ax.axhline(0, ls="--", color="0.6", lw=0.8)
            ax.set_ylabel("SNR")

        else:                               # 'blinking': ON=1 的竖线
            on_x = x_axis[is_on]
            ax.vlines(on_x, 0, 1, color=(0, 0.447, 0.741), lw=1.0)
            ax.plot(on_x, np.ones_like(on_x), ".", color=(0, 0.447, 0.741), ms=4)
            ax.set_ylim(0, 1.2)
            ax.set_yticks([0, 1])
            ax.set_ylabel("Emitter ON")
            ax.set_title(f"Blinking  (duty cycle = {is_on.mean():.2f})",
                         fontsize=10)

    axes[-1].set_xlabel(x_label)

    # ---- 轨迹图(右列)----
    if draw_track:
        x = np.asarray(emitter.x, dtype=float)
        y = np.asarray(emitter.y, dtype=float)
        valid = ~np.isnan(x) & ~np.isnan(y)

        axt.plot(x[valid], y[valid], "-o", lw=1.2, ms=4, label="Track")
        if valid.any():
            xi, yi = x[valid], y[valid]
            axt.plot(xi[0], yi[0], "o", mfc="g", mec="g", ms=9, label="Start")
            axt.plot(xi[-1], yi[-1], "o", mfc="r", mec="r", ms=9, label="End")
        axt.invert_yaxis()          # image 坐标:row 向下增长
        axt.set_aspect("equal")
        axt.grid(True, alpha=0.3)
        axt.set_xlabel("X (pixel)")
        axt.set_ylabel("Y (pixel)")
        axt.set_title(f"Track  (ecc={emitter.eccentricity():.2f})", fontsize=10)
        axt.legend(loc="best", fontsize=9)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  ✓ saved -> {save_path}")

    plt.show()
    return fig


def _emitter_movie_kwargs(emitter, traces, locs, index, extra):
    """组装单个分子的 ``save_movie`` 关键字: trace / overlay + 用户额外参数。"""
    kw = {}
    if traces is not None:
        kw['trace'] = traces[index]
    if locs is not None:
        li = getattr(emitter, 'loc_idx', None)
        kw['overlay'] = locs[li] if li is not None else locs
    kw.update(extra)                               # 用户显式给的优先
    return kw


def preview_emitters(stack, traces, emitters, locs=None, indices=None, *,
                     half=12, frames=None, stride=1, center='median',
                     width=420, ncols=2, **kwargs):
    """**一次预览多个**单分子窗口, 每个都和 ``save_movie`` 同款渲染。

    相当于对 ``indices`` 里的每个序号调用 :func:`smlm.io.preview_movie`, 再拼成
    一个可滚动的图墙; ``traces`` / ``locs`` 与 ``emitters`` 同索引时自动取对应
    的 trace 与定位点, 不用手写循环。

    Parameters
    ----------
    stack : ImageStack
        全 FOV 影片 (通常是漂移校正后的 ``corr_image``)。
    traces : sequence | None
        与 ``emitters`` 同索引的 trace 序列; None = 不画右侧曲线。
    emitters : sequence
        Emitter 列表, ``indices`` 是它的下标。
    locs : Localizations | None
        用于画 ``×`` 的定位点; None = 不叠加。有 ``emitter.loc_idx`` 时只取该分子的点。
    indices : None | int | slice | 序列 | 布尔掩码
        要预览哪些 emitter, 语义同 :func:`_as_index_list`; None = 全部。
    half, frames, stride, center
        窗口参数, 透传 :func:`emitter_roi`。
    width : int
        每个播放器的显示宽度 (像素)。
    ncols : int
        每行放几个; 只影响外框宽度, 不会裁剪视频。
    **kwargs
        其余透传 :func:`smlm.io.preview_movie` —— ``exposure_time`` / ``speed`` /
        ``clim`` ...

    Returns
    -------
    IPython.display.HTML
        图墙; notebook 里直接显示。
    """
    from IPython.display import HTML

    from .io import preview_movie

    idxs = _as_index_list(indices, len(emitters))
    extra = dict(kwargs)
    extra.setdefault('width', int(width))

    cells = []
    for i in idxs:
        roi = emitter_roi(stack, emitters[i], half=half, frames=frames,
                          stride=stride, center=center)
        kw = _emitter_movie_kwargs(emitters[i], traces, locs, i, extra)
        video = preview_movie(roi, **kw)._repr_html_()
        cells.append(
            '<div style="display:inline-block;margin:0 8px 8px 0;'
            'vertical-align:top">'
            f'<div style="font:12px/1.4 monospace;color:#666">#{i}</div>'
            f'{video}</div>')

    box = f'<div style="max-width:{ncols * (int(width) + 16)}px">'
    return HTML(box + ''.join(cells) + '</div>')


def save_emitters(stack, traces, emitters, locs=None, indices=None, *,
                  out_dir='.', pattern='em{idx:03d}.mp4', tag=None,
                  half=12, frames=None, stride=1, center='median', **kwargs):
    """把多个单分子窗口**批量**存成 mp4, 返回写出的路径列表。

    与 :func:`preview_emitters` 用同一套渲染, 只是落盘; 文件名由 ``pattern`` 决定,
    可用 ``{tag}`` 和 ``{idx}`` 两个占位符 (``tag`` 用 ``tag=`` 传入)。

    Examples
    --------
    >>> save_emitters(corr_image, traces, emitters, corr_locs, indices=KEEP,
    ...               out_dir=os.path.dirname(src), pattern='{tag}_em{idx:03d}.mp4',
    ...               tag='sif001', exposure_time=0.05, speed=5)
    """
    import os

    from .io import save_movie

    idxs = _as_index_list(indices, len(emitters))
    os.makedirs(out_dir, exist_ok=True)
    extra = dict(kwargs)

    paths = []
    for i in idxs:
        roi = emitter_roi(stack, emitters[i], half=half, frames=frames,
                          stride=stride, center=center)
        name = pattern.format(tag=tag if tag is not None else '', idx=i)
        path = os.path.join(out_dir, name)
        save_movie(roi, path, **_emitter_movie_kwargs(emitters[i], traces, locs,
                                                      i, extra))
        paths.append(path)
    return paths