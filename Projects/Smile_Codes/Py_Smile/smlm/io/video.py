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
from collections.abc import Mapping
from functools import lru_cache
import numpy as np

__all__ = ['save_movie']


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


def _cmap_lut(cmap, dpi=100):
    """实测 256 级查色表, 与 ``canvas.draw()`` 的输出逐字节一致。

    matplotlib 的路径是 ``lut[floor(g*256)]`` 再由 Agg 把浮点色值截断成 byte
    (不是四舍五入), 所以直接调 ``cmap(np.linspace(0,1,256))`` 会差 ±1。
    这里渲染一张 256 像素宽的 ramp 图把结果抄下来: 一次校准, 之后逐帧只做
    数组索引。
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    ramp = np.repeat(np.linspace(0.0, 1.0, 256)[None, :], 2, axis=0)
    fig = Figure(figsize=(256 / dpi, 2 / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(ramp, cmap='gray' if cmap is None else cmap, vmin=0.0, vmax=1.0,
              origin='upper', interpolation='nearest', aspect='auto')
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[0, :256, :3]
    return np.ascontiguousarray(buf)


@lru_cache(maxsize=16)
def _cached_lut(cmap, dpi=100):
    return _cmap_lut(cmap, dpi)


def _quantize(frame, vmin, vmax):
    """帧 -> 0..255 色表索引 (与 matplotlib 的 Normalize + cmap 索引一致)。"""
    arr = np.asarray(frame, dtype=float)
    span = float(vmax) - float(vmin)
    if span <= 0:                              # 常数图: 全部落到最暗
        return np.zeros(arr.shape, np.uint8)
    q = np.clip((arr - vmin) / span, 0.0, 1.0)
    q = np.nan_to_num(q, nan=0.0) * 256.0
    return np.minimum(np.floor(q), 255.0).astype(np.uint8)


def _to_rgb(frame, lut, vmin, vmax):
    return lut[_quantize(frame, vmin, vmax)]


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

    支持四种写法:
      - Trace 对象(.frames / .intensity, 即扣背景后的强度);
      - 字段字典 ``{'frames': ..., 'intensity': ...}`` (Trace 本身就是这种);
      - (frames, values) 二元组;
      - 1 维强度数组(帧号默认 0..n-1)。
    """
    if hasattr(trace, 'intensity'):
        v = np.asarray(trace.intensity, dtype=float)
        f = getattr(trace, 'frames', None)
        f = np.arange(v.size) if f is None else np.asarray(f, dtype=int)
        return f, v
    if isinstance(trace, Mapping) and 'intensity' in trace:
        v = np.asarray(trace['intensity'], dtype=float)
        f = trace.get('frames')
        f = np.arange(v.size) if f is None else np.asarray(f, dtype=int)
        return f, v
    if isinstance(trace, (tuple, list)) and len(trace) == 2:
        return np.asarray(trace[0], dtype=int), np.asarray(trace[1], dtype=float)
    v = np.asarray(trace, dtype=float)
    if v.ndim != 1:
        raise ValueError(
            "trace 需要是 Trace 对象 / 字段字典 / (frames, values) / 一维强度数组, "
            f"现在 shape={v.shape}。")
    return np.arange(v.size), v


class _FigureRenderer:
    """逐帧渲染: 图像 [+ 定位点] [+ 右侧强度 trace 面板]。

    做法: matplotlib 只在 ``__init__`` 里画**一次**静态层 (坐标轴 / 网格 / 刻度 /
    整条曲线 / 标签), 之后每帧用纯 numpy 合成变化的部分::

        out = 静态层.copy()
        图像面板  <- 逐帧像素 -> 灰度 -> 最近邻放大
        定位点    <- 预渲染的 marker 精灵, 按 alpha 贴图
        光标/圆点 <- 一列像素 + 预渲染的圆点精灵

    为什么不用 ``canvas.draw()``: 每帧重画整张图要 ~13 ms (1500 帧 = 19 s),
    而这里每帧 ~0.4 ms (1500 帧 < 1 s)。定位点/曲线一变就整图重绘是纯浪费
    —— 它们只占几个像素。
    """

    def __init__(self, h, w, cmap, vmin, vmax, overlay_kw=None,
                 trace_x=None, trace_y=None, x_label='Frame',
                 y_label='Intensity (bg-sub)', trace_panel=3.0,
                 trace_size=2.6, dpi=100, trace_kw=None):
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        self.dpi = dpi
        self.h, self.w = int(h), int(w)
        self.vmin, self.vmax = float(vmin), float(vmax)
        cmap = 'gray' if cmap is None else cmap
        self.lut = _cached_lut(cmap, float(dpi))
        self.fast = True                    # 几何退化时退回逐帧 canvas.draw()

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

        self._fit_trace_layout()
        self._build_static(kw, tkw if with_trace else None)

    # ---------- 布局 ----------
    def _fit_trace_layout(self):
        """把 trace 面板右移, 保证 y 轴标签不压到图像面板上。

        原来的 ``x0 = frac + 0.05*(1-frac)`` 只留了 ~15 px, 而 y 轴刻度标签 +
        ylabel 要 ~40 px, 标签会画到图像上。这里按实测边界自适应。
        """
        if self.ax_tr is None:
            return
        self.canvas.draw()
        fig_w_px = self.fig.get_figwidth() * self.dpi
        need = fig_w_px
        artists = [self.ax_tr.yaxis.label] + list(self.ax_tr.get_yticklabels())
        for art in artists:
            bb = art.get_window_extent()
            if bb.width and art.get_text():
                need = min(need, bb.x0)
        box = self.ax_tr.get_position()
        # 标签在轴左侧占的宽度, 平移后基本不变
        inset = max(0.0, box.x0 * fig_w_px - need)
        img_right = self.ax.get_position().x1 * fig_w_px
        want = (img_right + 6.0 + inset) / fig_w_px    # 6 px 间隙
        if want > box.x0 and box.x1 - want > 0.15:
            self.ax_tr.set_position([want, box.y0, box.x1 - want, box.height])

    # ---------- 一次性的静态层 + 几何量 ----------
    def _build_static(self, overlay_kw, trace_kw):
        for artist in (self.im, self.sc):
            artist.set_visible(False)
        if self.ax_tr is not None:
            self.cursor.set_visible(False)
            self.dot.set_visible(False)
        self.canvas.draw()
        rgba = np.asarray(self.canvas.buffer_rgba()).copy()
        self.base = rgba[..., :3].copy()
        self._H, self._W = self.base.shape[:2]

        # 图像面板: 面板像素 -> 图像像素索引
        box = self.ax.get_window_extent()
        r0, r1 = self._rows(box.y0, box.y1)
        c0, c1 = int(round(box.x0)), int(round(box.x1))
        r0, r1 = max(0, r0), min(self._H, r1)
        c0, c1 = max(0, c0), min(self._W, c1)
        # 图像面板整块覆盖: 图像轴 set_axis_off(), 且 _fit_trace_layout 保证
        # 没有别的 artist 画进这块区域。
        self.panel = (r0, r1, c0, c1)

        # 单行/单列图像没有 "数据单位->像素" 的可逆几何, 退回逐帧重绘
        if self.w < 2 or self.h < 2:
            self.fast = False
            return
        # 面板太小时 matplotlib 自己的重采样不再是纯最近邻 (实测 <8 px),
        # 这种情况退回逐帧重绘, 保证两条路径输出一致。
        if r1 - r0 < 8 or c1 - c0 < 8:
            self.fast = False
            return
        axx, axb, ayy, ayb = self._affine(self.ax, (0, 0), (self.w - 1,
                                                            self.h - 1))
        if not (np.isfinite(axx) and np.isfinite(ayy)) or axx == 0 or ayy == 0:
            self.fast = False
            return

        # 画布像素中心 -> 数据坐标 -> 最近邻索引 (与 matplotlib 一致)
        # cols 是画布列号(原点在左), rows 是画布行号(原点在上);
        # transData 的 y 原点在下, 故需 H-r-0.5 翻转。
        # 边界恰为半像素时 matplotlib 取较小索引, 故用 ceil(x-0.5) 而非
        # floor(x+0.5)。
        cx = np.arange(c0, c1) + 0.5
        cy = self._H - np.arange(r0, r1) - 0.5
        self._j_of = np.clip(np.ceil((cx - axb) / axx - 0.5), 0,
                             self.w - 1).astype(int)
        self._i_of = np.clip(np.ceil((cy - ayb) / ayy - 0.5), 0,
                             self.h - 1).astype(int)

        # 定位点: 图像数据坐标 (像素索引) -> 画布像素
        self._ov = self._affine(self.ax, (0, 0), (1, 1))
        kw = dict(c='r', marker='x', s=30, linewidths=0.5)
        kw.update(overlay_kw or {})
        self.ov_sprite = _marker_sprite(kw, self.dpi, scatter=True)

        # trace 面板: 数据坐标 -> 画布像素, 光标可用的行范围
        if self.ax_tr is not None:
            self._tr = self._affine(self.ax_tr, (0.0, 0.0), (1.0, 1.0))
            tb = self.ax_tr.get_window_extent()
            self._cy0, self._cy1 = self._rows(tb.y0, tb.y1)
            self.dot_sprite = _marker_sprite(
                dict(marker='o', ms=5, color=trace_kw['color']), self.dpi)
            self._cursor_color = (89, 89, 89)

    def _rows(self, y0, y1):
        """显示坐标 (原点在下) -> 画布行号 (原点在上)。"""
        return (self._H - int(round(y1)), self._H - int(round(y0)))

    @staticmethod
    def _affine(ax, p0, p1):
        """``transData`` 的线性化: 返回 (ax, b, ay, b) 使 px=ax*x+b, py=ay*y+b。"""
        q0 = ax.transData.transform(np.asarray(p0, float))
        q1 = ax.transData.transform(np.asarray(p1, float))
        ax_c = (q1[0] - q0[0]) / (float(p1[0]) - float(p0[0]))
        ay_c = (q1[1] - q0[1]) / (float(p1[1]) - float(p0[1]))
        return float(ax_c), float(q0[0] - ax_c * p0[0]), \
            float(ay_c), float(q0[1] - ay_c * p0[1])

    # ---------- 每帧 ----------
    def _frame_rgb(self, frame):
        arr = np.asarray(frame, dtype=float)
        return (self.lut[_quantize(arr, self.vmin, self.vmax)],
                ~np.isfinite(arr))

    def __call__(self, frame, pts, j):
        if not self.fast:
            return self._render_reference(frame, pts, j)

        out = self.base.copy()
        r0, r1, c0, c1 = self.panel
        rgb, bad = self._frame_rgb(frame)
        view = out[r0:r1, c0:c1]
        view[...] = rgb[self._i_of][:, self._j_of]
        if bad.any():
            # NaN 在 matplotlib 里是"透明", 也就是露出静态层, 不能画成黑色
            up = bad[self._i_of][:, self._j_of]
            if up.any():
                view[up] = self.base[r0:r1, c0:c1][up]

        if pts is not None and len(pts):
            pts = np.asarray(pts, float).reshape(-1, 2)
            ax_c, ax_b, ay_c, ay_b = self._ov
            xs = ax_c * pts[:, 0] + ax_b
            ys = ay_c * pts[:, 1] + ay_b
            for x, y in zip(xs, ys):
                _blit(out, self.ov_sprite, _snap(self._H - y), _snap(x),
                      self.panel)

        if self.ax_tr is not None:
            tx, tb, ty, td = self._tr
            xv = tx * float(self.trace_x[j]) + tb
            cx = _snap(xv)
            if 0 <= cx < self._W:
                # 虚线光标 (4 亮 4 暗), 与原来的 ls='--' 视觉一致
                rows = np.arange(self._cy0, self._cy1)
                on = (rows - self._cy0) % 8 < 4
                out[rows[on], cx] = self._cursor_color
            yv = self.trace_y[j]
            if np.isfinite(yv):
                yd = ty * yv + td
                _blit(out, self.dot_sprite, _snap(self._H - yd), cx,
                      (self._cy0, self._cy1, 0, self._W))
        return out

    def _render_reference(self, frame, pts, j):
        """逐帧重画整张图 (旧实现): 只在几何退化时用于回退, 测试也用它对照。"""
        self.im.set_visible(True)
        self.sc.set_visible(True)
        self.im.set_data(frame)
        self.sc.set_offsets(np.empty((0, 2)) if pts is None else pts)
        if self.ax_tr is not None:
            self.cursor.set_visible(True)
            self.dot.set_visible(True)
            xv = float(self.trace_x[j])
            self.cursor.set_xdata([xv, xv])
            yv = self.trace_y[j]
            self.dot.set_data([xv], [yv] if np.isfinite(yv) else [])
        self.canvas.draw()
        return np.asarray(self.canvas.buffer_rgba())[..., :3].copy()


def _marker_sprite(kw, dpi, scatter=False, px=48):
    """把一个 marker 渲染成 (rgb, alpha) 精灵, 供逐帧 alpha 贴图使用。

    ``scatter=True`` 时按 ``ax.scatter`` 的语义解释 (``s`` = 面积, 单位 pt²),
    否则按 ``ax.plot`` 的语义 (``ms`` = 直径, 单位 pt)。
    裁剪框对 ``px/2`` (像素边界) **对称**, 这样 ``_blit`` 按 ``sh//2`` 对位
    就是精确的 —— matplotlib 会把 marker 吸附到最近的像素中心, 没有亚像素。
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = Figure(figsize=(px / dpi, px / dpi), dpi=dpi, facecolor='none')
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.patch.set_alpha(0.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if scatter:
        ax.scatter([0.5], [0.5], **kw)
    else:
        ax.plot([0.5], [0.5], linestyle='none', **kw)
    canvas.draw()
    rgba = np.asarray(canvas.buffer_rgba())
    alpha = rgba[..., 3]
    ys, xs = np.nonzero(alpha)
    if ys.size == 0:
        return np.zeros((0, 0, 3), np.uint8), np.zeros((0, 0), np.uint8)
    mid = px // 2
    half = int(max(mid - ys.min(), ys.max() + 1 - mid,
                   mid - xs.min(), xs.max() + 1 - mid))
    y0, y1 = mid - half, mid + half
    x0, x1 = mid - half, mid + half
    return (rgba[y0:y1, x0:x1, :3].copy(),
            alpha[y0:y1, x0:x1].copy())


def _snap(v):
    """画布坐标 -> 像素索引: matplotlib 把点吸附到最近像素中心, 平局取上。"""
    return int(np.floor(v + 0.5))


def _blit(out, sprite, cy, cx, clip):
    """把精灵以 (cy, cx) 为中心 alpha 混合到 ``out``, 只写在 ``clip`` 内。"""
    rgb, alpha = sprite
    sh, sw = alpha.shape
    if sh == 0:
        return
    r0, r1, c0, c1 = clip
    y0, x0 = cy - sh // 2, cx - sw // 2
    ys0, xs0 = max(y0, r0), max(x0, c0)
    ys1, xs1 = min(y0 + sh, r1), min(x0 + sw, c1)
    if ys0 >= ys1 or xs0 >= xs1:
        return
    a = alpha[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0]
    if not a.any():
        return
    src = rgb[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0].astype(np.int16)
    a16 = a.astype(np.int16)[..., None]
    dst = out[ys0:ys1, xs0:xs1].astype(np.int16)
    blended = (dst * (255 - a16) + src * a16 + 127) // 255
    out[ys0:ys1, xs0:xs1] = blended.astype(np.uint8)



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

    lut = _cached_lut(cmap, 100.0)

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


def preview_movie(data, path=None, *, width=480, **kwargs):
    """在 notebook 里**先看**视频, 满意再用 ``save_movie`` 存。

    渲染完全走 :func:`save_movie` (同样的对比度、定位点、右侧同步曲线), 只是把结果
    内联成 ``<video>`` 直接显示。不给 ``path`` 时写一个临时 mp4, 读完字节就删掉,
    不会在磁盘上留文件。

    Parameters
    ----------
    data : ImageStack | RawDataset | np.ndarray
        影片; 单个分子通常传 ``smlm.viz.emitter_roi(corr_image, emitter)``。
    path : str, optional
        给了就同时把这个视频存到这个路径 (等同于再调一次 ``save_movie``)。
    width : int
        播放器显示宽度 (像素)。
    **kwargs
        其余全部透传 :func:`save_movie` —— ``exposure_time`` / ``speed`` /
        ``frames`` / ``overlay`` / ``trace`` / ``clim`` ...

    Returns
    -------
    IPython.display.HTML
        内联视频, notebook 里直接播放。

    Examples
    --------
    >>> roi = emitter_roi(corr_image, emitters[i], half=12)
    >>> preview_movie(roi, exposure_time=exposure_s, speed=5,
    ...               overlay=locs_i, trace=fill_traces[i])
    """
    import base64
    import os
    import tempfile

    from IPython.display import HTML

    kwargs.setdefault('verbose', False)        # 预览默认安静

    owner = None
    if path is None:
        fd, target = tempfile.mkstemp(prefix='smlm_preview_', suffix='.mp4')
        os.close(fd)
        owner = target
    else:
        target = os.fspath(path)

    try:
        save_movie(data, target, **kwargs)
        with open(target, 'rb') as fh:
            payload = fh.read()
    finally:
        if owner is not None:
            try:
                os.unlink(owner)
            except OSError:
                pass

    b64 = base64.b64encode(payload).decode('ascii')
    return HTML(
        f'<video width="{int(width)}" controls loop autoplay muted '
        f'src="data:video/mp4;base64,{b64}"></video>')
