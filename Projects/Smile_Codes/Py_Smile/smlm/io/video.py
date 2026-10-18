# smlm/io/video.py
"""把图像栈导出成视频(默认 MP4 / H.264, 拖进 PowerPoint 即可播放)。

帧率默认按真实采集速率:  ``fps = speed / exposure_time``

    exposure_time=0.1 s, speed=1  ->  10 fps  (原速, 默认)
    exposure_time=0.1 s, speed=5  ->  50 fps  (5 倍速)

对比度拉伸在全部帧上只算一次(默认 0.1/99.9 分位), 这样整片亮度不会闪。

依赖 imageio + imageio-ffmpeg(后者带一个静态 ffmpeg 可执行文件)::

    pip install imageio-ffmpeg
"""
from __future__ import annotations
import os
import numpy as np

__all__ = ['save_movie']

_GRAY_CMAPS = (None, 'gray', 'grey')


def _as_stack(data):
    """ImageStack / RawDataset / ndarray -> (T, H, W) 数组。"""
    arr = np.asarray(getattr(data, 'data', data))
    if arr.ndim == 2:
        arr = arr[None]
    if arr.ndim != 3:
        raise ValueError(
            f"data 需要 (T,H,W) 影片 / (H,W) 单图 / ImageStack, 现在 ndim={arr.ndim}。")
    return arr


def _resolve_frame_idx(n_frames, frames, stride):
    """把 frames 规整成整数索引数组。None -> 全部帧。"""
    if frames is None:
        idx = np.arange(n_frames)
    elif isinstance(frames, slice):
        idx = np.arange(n_frames)[frames]
    else:
        sel = np.asarray(frames)
        if sel.dtype == bool:
            if sel.size != n_frames:
                raise ValueError(
                    f"布尔掩码长度 {sel.size} 与帧数 {n_frames} 不一致。")
            idx = np.flatnonzero(sel)
        else:
            idx = np.atleast_1d(sel).astype(int)
    if idx.size and (idx.min() < 0 or idx.max() >= n_frames):
        raise IndexError(f"帧索引越界, 合法范围 [0, {n_frames - 1}]。")
    if int(stride) > 1:
        idx = idx[::int(stride)]
    if idx.size == 0:
        raise ValueError("frames / stride 选完没有任何帧。")
    return idx


def _even_pad(rgb):
    """H.264 的 yuv420p 要求宽高都是偶数; 奇数时往右/下补 1 像素。"""
    h, w = rgb.shape[:2]
    if h % 2 == 0 and w % 2 == 0:
        return rgb
    return np.pad(rgb, ((0, h % 2), (0, w % 2), (0, 0)), mode='edge')


def _cmap_lut(cmap):
    """取 256 级 RGB 查色表(比逐帧调 matplotlib 快很多)。"""
    try:
        from matplotlib import colormaps
        cmap_obj = colormaps[cmap]
    except Exception:                          # matplotlib 老版本
        from matplotlib import cm
        cmap_obj = cm.get_cmap(cmap)
    lut = np.asarray(cmap_obj(np.linspace(0, 1, 256)))[:, :3]
    return (lut * 255.0 + 0.5).astype(np.uint8)


def _to_rgb(frame, lut, vmin, vmax):
    q = np.clip((frame - vmin) / (vmax - vmin), 0.0, 1.0)
    q = (q * 255.0 + 0.5).astype(np.uint8)
    if lut is None:                            # 灰度: 省一份 LUT 拷贝
        return np.repeat(q[:, :, None], 3, axis=2)
    return lut[q]


def _stack_meta(data, n_frames):
    """从 ImageStack.metadata 里取窗口偏移与原始帧号。

    Crop / emitter_roi 会在 metadata 里写 crop_offset(y0,x0) 和
    source_frames(原始帧号)。没有这两项时退化为"不偏移 + 帧号 0..n-1"。
    """
    meta = getattr(data, 'metadata', None) or {}
    off = meta.get('crop_offset') or (0, 0)
    src = meta.get('source_frames')
    if src is None or len(src) != n_frames:
        src = np.arange(n_frames)
    return (int(off[0]), int(off[1])), np.asarray(src, dtype=int)


def _bucket_overlay(overlay, frame_map, offset=(0, 0)):
    """把 overlay 的 (frame, x, y) 按"输出帧序号"分桶。

    frame_map : {原始帧号 -> 输出帧序号}, 由 frames 选择决定。
    offset    : (y0, x0) 窗口偏移, 从 overlay 坐标里减掉, 换到窗口坐标系。
    """
    oy, ox = offset
    fr = np.asarray(overlay.frame).astype(int)
    x = np.asarray(overlay.x, dtype=float) - ox
    y = np.asarray(overlay.y, dtype=float) - oy
    buckets = {}
    for f, xi, yi in zip(fr, x, y):
        k = frame_map.get(int(f))
        if k is not None:
            buckets.setdefault(k, []).append((xi, yi))
    return {k: np.asarray(v, dtype=float) for k, v in buckets.items()}


def _auto_clim(sample, percentile=(0.1, 99.9), bg_sigma=5.0,
               signal_percentile=99.5):
    """从取样子集估显示范围 -> (vmin, vmax, stats)。

    为什么不能直接对整个体素用高分位: 单分子亮斑只占极小比例像素(典型
    <1%), 第 99.9 百分位基本=背景/噪声水平, 分子会被削成饱和白盘
    (信号像素可被削掉一半以上)。

    做法:
      1. 中位数 + MAD 估背景 med 与噪声 sigma(对亮点鲁棒);
      2. "信号像素" = 大于 med + bg_sigma*sigma 的像素;
      3. 信号像素足够 -> vmax = 它们在 signal_percentile 的分位(分子峰值
         的典型高度, 只削掉最亮的一小撮);
         几乎没有亮点 -> 退回全局 percentile[1];
      4. vmin = 全局 percentile[0](压背景)。
    """
    sample = np.asarray(sample, dtype=float)
    lo_pct, hi_pct = percentile
    vmin = float(np.percentile(sample, lo_pct))
    med = float(np.median(sample))
    sigma = 1.4826 * float(np.median(np.abs(sample - med)))

    if sigma > 0:
        sig = sample[sample > med + float(bg_sigma) * sigma]
    else:
        sig = sample[:0]

    if sig.size >= 100:
        vmax = float(np.percentile(sig, signal_percentile))
        clip = float((sig >= vmax).mean()) * 100.0
    else:                                   # 没有明显的稀疏亮点
        vmax = float(np.percentile(sample, hi_pct))
        clip = float((sample >= vmax).mean()) * 100.0

    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = float(sample.max())
        clip = 0.0
        if vmax <= vmin:
            vmax = vmin + 1.0
    return vmin, vmax, dict(median=med, sigma=sigma, n_sig=int(sig.size),
                            clip_pct=clip)


def auto_clim(data, *, frames=None, stride=1, percentile=(0.1, 99.9),
              bg_sigma=5.0, signal_percentile=99.5, max_frames=20,
              verbose=True):
    """给影片估一个"不会削掉分子"的显示范围 (vmin, vmax)。

    想手动微调时: 先调用它看看数值, 再把结果传给 save_movie(clim=...)。
    参数含义见 save_movie 的 clim / bg_sigma / signal_percentile。
    """
    arr = _as_stack(data)
    idx = _resolve_frame_idx(arr.shape[0], frames, stride)
    sample = arr[idx[:: max(1, idx.size // int(max_frames))]]
    vmin, vmax, st = _auto_clim(sample, percentile, bg_sigma,
                                signal_percentile)
    if verbose:
        print(f"[io] auto_clim: vmin={vmin:.4g}, vmax={vmax:.4g}  "
              f"(背景中位={st['median']:.4g}, 噪声σ={st['sigma']:.4g}, "
              f"信号像素={st['n_sig']}, 信号被削={st['clip_pct']:.2f}%)")
    return vmin, vmax


def _split_trace(trace):
    """把 trace 参数规整成 (帧号数组, 强度数组)。

    支持三种写法:
      - Trace 对象(.frames / .intensity, 即扣背景后的强度);
      - (frames, values) 二元组;
      - 1 维强度数组(帧号默认 0..n-1)。
    """
    if hasattr(trace, 'intensity'):
        v = np.asarray(trace.intensity, dtype=float)
        f = getattr(trace, 'frames', None)
        f = np.arange(v.size) if f is None else np.asarray(f, dtype=int)
        return f, v
    if isinstance(trace, (tuple, list)) and len(trace) == 2:
        return np.asarray(trace[0], dtype=int), np.asarray(trace[1], dtype=float)
    v = np.asarray(trace, dtype=float)
    if v.ndim != 1:
        raise ValueError(
            "trace 需要是 Trace 对象 / (frames, values) / 一维强度数组, "
            f"现在 shape={v.shape}。")
    return np.arange(v.size), v


class _FigureRenderer:
    """逐帧渲染: 图像 [+ 定位点] [+ 右侧强度 trace 面板]。

    只有需要画附加内容(定位点 / trace)时才走 matplotlib; 否则 save_movie
    用纯 numpy 的快速路径。
    """

    def __init__(self, h, w, cmap, vmin, vmax, overlay_kw=None,
                 trace_x=None, trace_y=None, x_label='Frame',
                 y_label='Intensity (bg-sub)', trace_panel=3.0,
                 trace_size=2.6, dpi=100, trace_kw=None):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        self.dpi = dpi
        with_trace = trace_x is not None
        # 有 trace 时图像面板用固定物理尺寸: 小窗口 ROI 会被放大显示
        img_w_in = trace_size if with_trace else w / dpi
        img_h_in = (trace_size * h / w) if with_trace else h / dpi
        fig_w = img_w_in + (trace_panel if with_trace else 0.0)
        self.fig = Figure(figsize=(fig_w, max(img_h_in, 1e-3)), dpi=dpi)
        self.canvas = FigureCanvasAgg(self.fig)

        frac = img_w_in / fig_w
        ax = self.fig.add_axes([0, 0, frac, 1.0])      # 满幅, 不留边距
        ax.set_axis_off()
        self.im = ax.imshow(np.zeros((h, w)), cmap=cmap, vmin=vmin, vmax=vmax,
                            origin='upper', interpolation='nearest')
        ax.set_aspect('auto')                          # 面板宽高比 = 图像宽高比
        kw = dict(c='r', marker='x', s=30, linewidths=0.5)
        kw.update(overlay_kw or {})
        self.sc = ax.scatter([], [], **kw)
        self.ax = ax

        self.trace_x = None
        self.ax_tr = None
        if with_trace:
            tkw = dict(color=(0, 0.447, 0.741), lw=1.2)
            tkw.update(trace_kw or {})
            x0 = frac + 0.05 * (1 - frac)
            self.ax_tr = self.fig.add_axes([x0, 0.17, 0.96 - x0, 0.74])
            self.trace_x = np.asarray(trace_x, dtype=float)
            self.trace_y = np.asarray(trace_y, dtype=float)
            self.ax_tr.plot(self.trace_x, self.trace_y, **tkw)
            self.ax_tr.axhline(0, ls='--', color='0.7', lw=0.8)
            self.cursor = self.ax_tr.axvline(float(self.trace_x[0]),
                                             color='0.35', lw=1.0, ls='--')
            self.dot, = self.ax_tr.plot([], [], 'o', ms=5,
                                        color=tkw['color'], zorder=4)
            self.ax_tr.set_xlabel(x_label, fontsize=9)
            self.ax_tr.set_ylabel(y_label, fontsize=9)
            self.ax_tr.tick_params(labelsize=8)
            self.ax_tr.grid(True, alpha=0.25)

    def __call__(self, frame, pts, j):
        self.im.set_data(frame)
        self.sc.set_offsets(np.empty((0, 2)) if pts is None else pts)
        if self.ax_tr is not None:
            xv = float(self.trace_x[j])
            self.cursor.set_xdata([xv, xv])
            yv = self.trace_y[j]
            self.dot.set_data([xv], [yv] if np.isfinite(yv) else [])
        self.canvas.draw()
        return np.asarray(self.canvas.buffer_rgba())[..., :3].copy()


def save_movie(data, path, *, exposure_time=None, fps=None, speed=1.0,
               frames=None, stride=1, clim=None, percentile=(0.1, 99.9),
               bg_sigma=5.0, signal_percentile=99.5,
               cmap='gray', overlay=None, overlay_kw=None,
               trace=None, trace_kw=None, trace_panel=3.0, trace_size=2.6,
               trace_y_label='Intensity (bg-sub)',
               codec='libx264', quality=8, verbose=True):
    """把图像栈写成视频(默认 MP4/H.264/yuv420p, 可直接插进 PowerPoint)。

    Parameters
    ----------
    data : ImageStack | RawDataset | np.ndarray
        影片 (T,H,W) 或单图 (H,W)。ImageStack/RawDataset 会自动取 .data。
    path : str
        输出路径, 建议 .mp4。
    exposure_time : float, optional
        单帧曝光时间(秒)。给了它就按真实速率播: fps = speed / exposure_time。
    fps : float, optional
        直接指定帧率, 优先于 exposure_time。两者都缺会报错。
    speed : float
        倍速。默认 1.0 = 原速; 2.0 = 两倍快, 0.5 = 放慢一半。
    frames : None | slice | int | 索引数组 | 布尔掩码
        要写进视频的帧(相对于 data 的帧序号)。默认 None = 全部帧。
        例如 slice(0, 300) / np.arange(0, 2000, 2) / (np.arange(T) < 300)。
    stride : int
        在 frames 基础上再等间隔抽帧, 默认 1(不抽)。
    clim : (vmin, vmax), optional
        手动灰度范围; 默认 None -> 自动估(见下方说明)。
    percentile : (lo, hi)
        自动估里用到的全局分位点, 默认 (0.1, 99.9): lo 用来压背景,
        hi 只在"几乎没有亮点"时当退回值。
    bg_sigma : float
        背景噪声倍数, 默认 5.0。像素 > 背景中位 + bg_sigma*σ 才算"信号
        像素"(σ 由 MAD 鲁棒估计)。
    signal_percentile : float
        vmax = 信号像素的该分位, 默认 99.5。也就是说只削掉最亮的一小撮
        分子峰顶, 而不是把整片分子压成白盘。

        为什么要这样: 单分子亮斑只占极小比例像素(典型 <1%), 若直接对
        全体素取 99.9 分位, 拿到的基本是背景水平 —— 分子会严重过曝
        (实测稀疏场景下可剔掉 60% 的信号像素)。想完全不过曝就显式传
        clim=(vmin, vmax) 或 clim=(vmin, data.max())。
    cmap : str | None
        伪彩, 默认 'gray'。可传 matplotlib 的名字(如 'hot' / 'viridis')。
    overlay : Localizations | Candidates, optional
        给了就在每帧上叠加该帧的定位点(需要 .frame / .x / .y)。
        frame 号与坐标会自动对齐: 若 data 是被 Crop / emitter_roi 裁过的
        ImageStack, 会读它的 metadata['source_frames'] / ['crop_offset']
        把定位点按原始帧号匹配、并平移到窗口坐标。
    overlay_kw : dict, optional
        叠加点的散点样式, 默认 {'c':'r', 'marker':'x', 's':30, 'linewidths':0.5}。
    trace : Trace | (frames, values) | np.ndarray, optional
        给了就在视频右侧加一条强度曲线(通常直接传 Trace, 用它的
        ``intensity``, 即扣背景后的信号), 并画一条随帧移动的竖线 + 圆点
        指示当前帧。三种写法:
          - Trace 对象: 取 .frames / .intensity;
          - (frames, values) 二元组;
          - 一维强度数组(帧号默认 0..n-1)。
        帧号按 ImageStack 的 metadata['source_frames'] 对齐, 所以配合
        emitter_roi 裁过的窗口也能正确对上。
    trace_kw : dict, optional
        曲线样式, 默认 {'color': (0, 0.447, 0.741), 'lw': 1.2}。
    trace_panel : float
        右侧曲线面板宽度(英寸), 默认 3.0。
    trace_size : float
        有 trace 时图像面板的宽度(英寸), 默认 2.6。小窗口 ROI 会被放大
        显示, 免得 15x15 的图挤在角落里。
    trace_y_label : str
        曲线 y 轴标签, 默认 'Intensity (bg-sub)'。
    codec : str
        ffmpeg 编码器, 默认 libx264。
    quality : int
        imageio 的 0-10 画质(越大越好), 默认 8。
    verbose : bool
        打印进度。

    Returns
    -------
    path : str
        写出的文件路径。
    """
    import imageio

    arr = _as_stack(data)
    n_frames, h, w = arr.shape
    idx = _resolve_frame_idx(n_frames, frames, stride)
    # ImageStack 可能已被 Crop / emitter_roi 裁过: 取窗口偏移和原始帧号
    (oy, ox), src_frames = _stack_meta(data, n_frames)

    # ---- 帧率: 默认真实速率 ----
    if fps is None:
        if exposure_time is None:
            raise ValueError(
                "需要 exposure_time(单帧曝光秒数) 或直接给 fps, 二者至少一个。")
        fps = float(speed) / float(exposure_time)
    fps = float(fps)
    if fps <= 0:
        raise ValueError(f"fps 必须为正, 现在是 {fps}。")

    # ---- 对比度: 只在取样帧上算一次, 避免逐帧闪烁 ----
    sample = arr[idx[:: max(1, idx.size // 20)]]
    clip_stats = None
    if clim is None:
        vmin, vmax, clip_stats = _auto_clim(sample, percentile, bg_sigma,
                                            signal_percentile)
    else:
        vmin, vmax = float(clim[0]), float(clim[1])

    lut = None if cmap in _GRAY_CMAPS else _cmap_lut(cmap)

    # ---- overlay 按输出帧序号分桶 ----
    shown_src = np.asarray(src_frames)[idx]        # 每个输出帧对应的原始帧号
    buckets, render = {}, None
    n_pts = 0
    if overlay is not None:
        frame_map = {int(f): i for i, f in enumerate(shown_src)}
        buckets = _bucket_overlay(overlay, frame_map, offset=(oy, ox))
        n_pts = int(len(np.asarray(overlay.frame)))

    # ---- 右侧强度曲线: 按原始帧号与视频帧对齐 ----
    tr_x = tr_y = None
    x_label, trace_pts = "Frame", 0
    if trace is not None:
        t_frames, t_vals = _split_trace(trace)
        pos = {int(f): k for k, f in enumerate(t_frames)}
        tr_y = np.array([t_vals[pos[int(f)]] if int(f) in pos else np.nan
                         for f in shown_src], dtype=float)
        trace_pts = int(np.isfinite(tr_y).sum())
        rel = shown_src.astype(float) - float(shown_src[0])
        if exposure_time is not None:
            tr_x, x_label = rel * float(exposure_time), "Time (s)"
        else:
            tr_x = rel

    if overlay is not None or trace is not None:
        render = _FigureRenderer(
            h, w, cmap if cmap is not None else 'gray', vmin, vmax,
            overlay_kw=overlay_kw, trace_x=tr_x, trace_y=tr_y,
            x_label=x_label, y_label=trace_y_label, trace_panel=trace_panel,
            trace_size=trace_size, trace_kw=trace_kw)

    path = os.fspath(path)
    parent = os.path.dirname(os.path.abspath(path))
    if parent and not os.path.isdir(parent):
        try:                                   # 目录缺失时先试着自己建
            os.makedirs(parent, exist_ok=True)
        except OSError as exc:
            raise FileNotFoundError(
                f"输出目录不存在, 也无法创建: {parent}\n"
                f"  ({exc})\n"
                f"  路径若在外置盘上, 请先确认盘已挂载。") from None

    writer = imageio.get_writer(path, fps=fps, codec=codec, quality=quality,
                                pixelformat='yuv420p', macro_block_size=None)
    if verbose:
        extra = (f", overlay={n_pts} 点" if overlay is not None else "") \
            + (f", trace={trace_pts} 点" if trace is not None else "")
        print(f"[io] {os.path.basename(path)}: {idx.size} 帧, {fps:g} fps, "
              f"{h}x{w}, frames=[{idx[0]}, {idx[-1]}]{extra}")
        if clip_stats is None:
            print(f"     clim=({vmin:.4g}, {vmax:.4g})  [手动]")
        else:
            print(f"     clim=({vmin:.4g}, {vmax:.4g})  "
                  f"[自动: 背景中位={clip_stats['median']:.4g}, "
                  f"噪声σ={clip_stats['sigma']:.4g}, "
                  f"信号像素={clip_stats['n_sig']}, "
                  f"其中被削={clip_stats['clip_pct']:.2f}%]"
                  f"   → 觉得还亮/还暗就传 clim= 手动定")

    try:
        total = idx.size
        step = max(1, total // 10)
        for j, i in enumerate(idx):
            frame = np.asarray(arr[i], dtype=float)
            if render is None:
                rgb = _to_rgb(frame, lut, vmin, vmax)
            else:
                rgb = render(frame, buckets.get(j), j)
            writer.append_data(_even_pad(rgb))
            if verbose and ((j + 1) % step == 0 or j + 1 == total):
                print(f"    {j + 1}/{total}")
    finally:
        writer.close()

    if verbose:
        print(f"[io] 完成 -> {path}  ({os.path.getsize(path) / 1e6:.1f} MB)")
    return path
