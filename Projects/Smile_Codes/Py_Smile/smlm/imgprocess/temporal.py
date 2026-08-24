# smlm/imgprocess/temporal.py
"""时间维度处理: 光漂白校正 / 漂移校正 (估计侧在 smlm/drift, 校正侧在本模块)。

SOFI / RIED 的前提是"宽平稳 + 无运动"的闪烁信号, 漂移与漂白都会破坏它。
这些函数输入输出都是 ImageStack (T,H,W) 或 numpy 数组, 供 pipeline 调用。

包含:
  - apply_drift : 按已知逐帧漂移做亚像素傅里叶配准 (轻量纯函数);
  - bleach_correct / BleachCorrect : 光漂白校正;
  - TemporalMean : 时间平均 (Transform 子类, 栈 -> 单帧均值投影);
  - DriftCorrector : 完整漂移校正器 (Transform 子类, 图像 + 定位, 支持裁切)。
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.fft import fft2, ifft2
from scipy.ndimage import fourier_shift
from ..datatypes import ImageStack, Localizations, DriftTrace
from .base import Transform


def apply_drift(stack, shifts):
    """按 shifts 把每帧配准回参考帧 (亚像素, 傅里叶移位)。

    stack : (T, H, W) float
    shifts: (T, 2) 每帧 (dy, dx) 累积漂移, 单位 pixel —— 与 smlm/drift 的
            DriftTrace.dy/dx 一致 (即"坐标要减去的量")。
    """
    stack = np.asarray(stack, dtype=np.float32)
    shifts = np.asarray(shifts, dtype=float)
    out = np.empty_like(stack)
    for t in range(stack.shape[0]):
        out[t] = np.real(ifft2(fourier_shift(fft2(stack[t]), -shifts[t])))
    return out


def bleach_correct(stack, mode="mono_exp"):
    """光漂白校正 —— SOFI/RIED 假设信号严格平稳, 漂白会破坏这一假设。

    mode="mono_exp": 对逐帧均值拟合单指数衰减, 再逐帧归一化 (最常用、最稳)。
    mode="none":     不做, 只靠 substack 分块。
    """
    if mode == "none":
        return stack
    stack = np.asarray(stack, dtype=np.float32)
    trace = stack.reshape(stack.shape[0], -1).mean(1)
    t = np.arange(len(trace), dtype=float)
    base = trace.min() * 0.9
    y = np.log(np.maximum(trace - base, 1e-6))
    a, b = np.polyfit(t, y, 1)                  # log y = a t + b
    fit = np.exp(b) * np.exp(a * t) + base
    return (stack * (fit[0] / fit)[:, None, None]).astype(np.float32)


class BleachCorrect(Transform):
    """光漂白校正 Transform (ImageStack -> ImageStack)。"""
    def __init__(self, mode="mono_exp"):
        self.mode = mode

    def __call__(self, stack):
        data = bleach_correct(np.asarray(stack.data, dtype=np.float32), self.mode)
        return replace(stack, data=data)

    def __repr__(self):
        return f"BleachCorrect(mode={self.mode!r})"


class TemporalMean(Transform):
    """时间平均 Transform：把整个栈按时间平均成单帧 (T,H,W) -> (1,H,W)。

    等价于 widefield mean 投影：SOFI/RIED 中常用它作为"宽场对照"。
    输出仍是 ImageStack（T=1），可直接取 .data[0] 或继续接 Pipeline。
    """
    def __init__(self, dtype=np.float32):
        self.dtype = dtype

    def __call__(self, stack):
        data = np.asarray(stack.data, dtype=self.dtype).mean(axis=0, keepdims=True)
        return replace(stack, data=data)

    def __repr__(self):
        return f"TemporalMean(dtype={self.dtype.__name__})"


# ======================================================================
# 漂移校正器 (由 smlm/drift/base.py 移入, 校正侧与估计侧分离)
# ======================================================================

def _index_locs(locs: Localizations, mask: np.ndarray) -> Localizations:
    """按布尔/整数索引筛选一个 Localizations 的所有列 (含 extra 中的数组)。"""
    def take(a):
        return None if a is None else np.asarray(a)[mask]

    extra = {}
    for k, v in (locs.extra or {}).items():
        try:
            if hasattr(v, '__len__') and len(v) == len(locs):
                extra[k] = np.asarray(v)[mask]
            else:
                extra[k] = v
        except TypeError:
            extra[k] = v

    return Localizations(
        frame=take(locs.frame), x=take(locs.x), y=take(locs.y),
        photons=take(locs.photons), bg=take(locs.bg),
        sigma=take(locs.sigma), lp=take(locs.lp),
        z=take(locs.z), z_lp=take(locs.z_lp), extra=extra,
    )


class DriftCorrector(Transform):
    """根据 DriftTrace 校正图像栈和定位坐标 (Transform 子类)。

    两种用法:
      - 作为 Transform 进 Pipeline: `DriftCorrector(trace=tr, crop=...)(stack)`
        返回校正后的 ImageStack (ImageStack 进, ImageStack 出);
      - 低层 API (兼容旧调用): `DriftCorrector(crop=...).apply(tr, image=stack,
        locs=locs)` 返回 (corr_image, corr_locs, keep_idx)。

    参数
    ----
    trace    : 可选的 DriftTrace。Transform 用法时必须在构造时绑定;
               apply()/correct_* 里也可显式传入以覆盖。
    crop     : (H_out, W_out) 或 None。给定时对平移后的图像做居中裁切,
               多出的部分作为平移缓冲区(避免边缘黑边混入)。定位坐标同步
               转换到裁切窗口的局部坐标系。None 表示不裁切、保持原尺寸。
    edge     : 裁切后边缘保护(pixel), 落在 edge 内的定位被丢弃。
    fill     : 图像平移后空缺像素的填充值。
    order    : 插值阶数, 1 = 双线性(与 MATLAB imtranslate 'linear' 对应)。
    """

    def __init__(self, trace=None, crop=None, edge=0, fill=0.0, order=1):
        super().__init__()
        self.trace = trace          # 构造时可绑定 DriftTrace (Transform 用法)
        self.crop = crop
        self.edge = edge
        self.fill = fill
        self.order = order

    def __call__(self, stack):
        """Transform 接口: ImageStack -> ImageStack (需构造时绑定 trace)。"""
        if self.trace is None:
            raise ValueError(
                "作为 Transform 使用需在构造时传入 trace=DriftTrace(...); "
                "或改用 .apply(trace, image=..., locs=...)。")
        return self.correct_image(stack, self.trace)[0]

    # ---- 图像 ----
    def correct_image(self, image: ImageStack, trace: DriftTrace):
        from scipy.ndimage import shift as nd_shift

        data = np.asarray(image.data)                  # (T, H, W)
        T, H, W = data.shape
        if trace.n_frames != T:
            raise ValueError(f"帧数不一致: image {T} vs trace {trace.n_frames}")

        # 裁切窗口
        if self.crop is None:
            H_out, W_out, r0, c0 = H, W, 0, 0
        else:
            H_out, W_out = self.crop
            if H_out > H or W_out > W:
                raise ValueError(
                    f"crop {self.crop} 超出图像尺寸 ({H},{W}), 无可用缓冲区。")
            r0 = (H - H_out) // 2
            c0 = (W - W_out) // 2

        out = np.empty((T, H_out, W_out), dtype=data.dtype)
        for t in range(T):
            # 抵消漂移: 图像反向平移 (-dy, -dx)
            shifted = nd_shift(data[t], shift=(-trace.dy[t], -trace.dx[t]),
                               order=self.order, mode='constant',
                               cval=self.fill)
            out[t] = shifted[r0:r0 + H_out, c0:c0 + W_out]

        return (ImageStack(data=out, pixel_size=image.pixel_size,
                           metadata={**image.metadata,
                                     'drift_corrected': True}),
                (r0, c0, H_out, W_out))

    # ---- 定位 ----
    def correct_locs(self, locs: Localizations, trace: DriftTrace,
                     window=None):
        frame = np.asarray(locs.frame).astype(int)
        x_new = np.asarray(locs.x, float) - trace.dx[frame]
        y_new = np.asarray(locs.y, float) - trace.dy[frame]

        if window is None:                              # 不裁切
            out = _index_locs(locs, np.ones(len(locs), bool))
            out.x, out.y = x_new, y_new
            keep_idx = np.arange(len(locs))
            return out, keep_idx

        r0, c0, H_out, W_out = window
        x_loc = x_new - c0                              # 转到局部坐标 (0-based)
        y_loc = y_new - r0

        e = self.edge
        keep = ((x_loc >= e) & (x_loc <= W_out - 1 - e) &
                (y_loc >= e) & (y_loc <= H_out - 1 - e))

        out = _index_locs(locs, keep)
        out.x = x_loc[keep]
        out.y = y_loc[keep]
        keep_idx = np.nonzero(keep)[0]
        return out, keep_idx

    # ---- 一步到位 ----
    def apply(self, trace: DriftTrace = None, image: ImageStack = None,
              locs: Localizations = None):
        """返回 (corr_image, corr_locs, keep_idx)。未提供的输入对应位置为 None。

        trace 未给时使用构造时绑定的 self.trace。
        """
        if trace is None:
            trace = self.trace
        if trace is None:
            raise ValueError(
                "缺少 DriftTrace: 请传入 apply(trace, ...) 或在构造时绑定 trace=...")

        window = None
        corr_image = None
        if image is not None:
            corr_image, window = self.correct_image(image, trace)

        corr_locs = keep_idx = None
        if locs is not None:
            corr_locs, keep_idx = self.correct_locs(locs, trace, window)

        return corr_image, corr_locs, keep_idx
