# smlm/imgprocess/interp.py
"""插值 / 上采样 / 下采样 Transform。

傅里叶插值是无损的: 原始像素已满足 Nyquist 时, OTF 之外频率本来就是 0,
补零不引入新信息。这样 m 阶 SOFI / RIED 的输出像素才配得上它的分辨率。
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.fft import fftn, ifftn, fftshift, ifftshift
from .base import Transform


def fourier_zoom(img, factor):
    """无损傅里叶插值: FFT -> 中心零填充 -> IFFT。支持 (..., H, W)。

    两个关键点 (缺一不可):
      - IFFT 自带 1/M 归一化, 零填充后需乘回 factor^2 (2D) 才能保持幅度;
      - 奇尺寸的 DC 对齐偏移要 +1, 否则 DC 落在半像素处 -> 相位误差/负值。
    """
    factor = int(factor)
    if factor == 1:
        return np.asarray(img, dtype=np.float32)
    s = np.asarray(img, dtype=np.float32)
    f = fftshift(fftn(s, axes=(-2, -1)), axes=(-2, -1))
    ny, nx = s.shape[-2:]
    Ny, Nx = ny * factor, nx * factor
    padded = np.zeros(s.shape[:-2] + (Ny, Nx), dtype=np.complex128)
    y0 = (Ny - ny + (ny % 2)) // 2
    x0 = (Nx - nx + (nx % 2)) // 2
    padded[..., y0:y0 + ny, x0:x0 + nx] = f
    out = ifftn(ifftshift(padded, axes=(-2, -1)), axes=(-2, -1))
    out = np.real(out) * (factor * factor)          # 幅度归一化 (2D)
    return np.ascontiguousarray(out, dtype=np.float32)


def fourier_upsample(stack, factor=2, taper=True):
    """fourier_zoom 的封装: 加余弦边缘窗, 抑制 FFT 周期性造成的边缘振铃。

    taper=True 会把图像四周压暗 (高阶 SOFI 按 w^m 加剧) —— 高上采样/高阶时
    建议 taper=False。
    """
    if factor == 1:
        return np.asarray(stack, dtype=np.float32)
    s = np.asarray(stack, dtype=np.float32)
    if taper:
        ny, nx = s.shape[-2], s.shape[-1]
        wy = np.hanning(max(ny // 8, 1) * 2 + 2)[1:-1]
        wx = np.hanning(max(nx // 8, 1) * 2 + 2)[1:-1]
        my = np.ones(ny)
        my[:len(wy) // 2] = wy[:len(wy) // 2]
        my[-len(wy) // 2:] = wy[-len(wy) // 2:]
        mx = np.ones(nx)
        mx[:len(wx) // 2] = wx[:len(wx) // 2]
        mx[-len(wx) // 2:] = wx[-len(wx) // 2:]
        s = s * (my[:, None] * mx[None, :])      # (H,W) 可广播到 (T,H,W)
    return fourier_zoom(s, factor)


class FourierUpsample(Transform):
    """傅里叶上采样: ImageStack(T,H,W) -> (T, H*factor, W*factor)。

    同时按 factor 缩小 pixel_size (nm/pixel), 并累乘 metadata['upsample']。
    加余弦边缘窗默认关闭 (taper=False), 避免压暗边缘。
    """
    def __init__(self, factor=2, taper=False):
        self.factor = int(factor)
        self.taper = bool(taper)

    def __call__(self, stack):
        data = fourier_upsample(stack.data, self.factor, taper=self.taper)
        meta = dict(stack.metadata)
        meta['upsample'] = int(meta.get('upsample', 1)) * self.factor
        ps = stack.pixel_size
        new_ps = (ps / self.factor) if ps is not None else None
        return replace(stack, data=data, pixel_size=new_ps, metadata=meta)

    def __repr__(self):
        return f"FourierUpsample(factor={self.factor}, taper={self.taper})"
