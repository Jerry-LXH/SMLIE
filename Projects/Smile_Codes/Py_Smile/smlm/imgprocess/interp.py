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
    大栈沿批(帧)轴自动分块 (chunk), 只改执行顺序, 数值与一次性处理完全一致;
    内存占用从"整栈"降到"一个分块"。
    """
    factor = int(factor)
    if factor == 1:
        return np.asarray(img, dtype=np.float32)
    s = np.asarray(img, dtype=np.float32)
    single = (s.ndim == 2)                 # 鸭子: 2D 单帧 -> 返回 2D
    if single:
        s = s[None, ...]
    ny, nx = s.shape[-2:]
    Ny, Nx = ny * factor, nx * factor
    y0 = (Ny - ny + (ny % 2)) // 2
    x0 = (Nx - nx + (nx % 2)) // 2

    batch = s.shape[:-2] or (1,)
    n_batch = int(np.prod(batch))
    flat = s.reshape(n_batch, ny, nx)
    out = np.empty((n_batch, Ny, Nx), dtype=np.float32)
    # 每帧估算: fftn/fftshift (ny*nx) + padded/ifftn (Ny*Nx) 的 complex128 + real
    per_frame = 32 * ny * nx + 40 * Ny * Nx
    chunk = max(1, min(n_batch, int(5.0e8 / max(per_frame, 1))))
    for s0 in range(0, n_batch, chunk):
        sub = flat[s0:s0 + chunk]
        f = fftshift(fftn(sub, axes=(-2, -1)), axes=(-2, -1))
        padded = np.zeros(sub.shape[:-2] + (Ny, Nx), dtype=np.complex128)
        padded[..., y0:y0 + ny, x0:x0 + nx] = f
        zo = np.real(ifftn(ifftshift(padded, axes=(-2, -1)), axes=(-2, -1)))
        out[s0:s0 + chunk] = zo * (factor * factor)     # 幅度归一化 (2D)
    res = np.ascontiguousarray(out.reshape(batch + (Ny, Nx)), dtype=np.float32)
    return res[0] if single else res


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


def fourier_interp_2d(img, newsz):
    """matlab `fInterp_2D` 等价: 任意目标尺寸的傅里叶插值 (四角频谱放置, 偶数轴 Nyquist 拆半)。

    支持 (..., H, W); newsz=(Ny, Nx) 目标尺寸 (>= 原尺寸)。幅度按 (Ny/ny)(Nx/nx)
    缩放以守恒原像素计数。
    """
    s = np.asarray(img, dtype=np.float64)
    ny, nx = s.shape[-2], s.shape[-1]
    Ny, Nx = int(newsz[0]), int(newsz[1])
    F = np.fft.fft2(s, axes=(-2, -1)) * ((Ny / ny) * (Nx / nx))
    Fp = np.zeros(s.shape[:-2] + (Ny, Nx), dtype=np.complex128)
    nyqy = int(np.ceil((ny + 1) / 2))           # 1-based Nyquist 索引 (matlab 约定)
    nyqx = int(np.ceil((nx + 1) / 2))
    Fp[..., 0:nyqy, 0:nyqx] = F[..., 0:nyqy, 0:nyqx]                       # xl, yl
    Fp[..., Ny - (ny - nyqy):, 0:nyqx] = F[..., nyqy:, 0:nyqx]             # xh, yl
    Fp[..., 0:nyqy, Nx - (nx - nyqx):] = F[..., 0:nyqy, nyqx:]             # xl, yh
    Fp[..., Ny - (ny - nyqy):, Nx - (nx - nyqx):] = F[..., nyqy:, nyqx:]   # xh, yh
    if ny % 2 == 0 and Ny != ny:                # 偶数轴: Nyquist 频率拆成两半
        Fp[..., nyqy - 1, :] *= 0.5
        Fp[..., nyqy - 1 + (Ny - ny), :] = Fp[..., nyqy - 1, :]
    if nx % 2 == 0 and Nx != nx:
        Fp[..., :, nyqx - 1] *= 0.5
        Fp[..., :, nyqx - 1 + (Nx - nx)] = Fp[..., :, nyqx - 1]
    return np.real(np.fft.ifft2(Fp, axes=(-2, -1)))


def fourier_lateral(img, factor=2):
    """matlab `fourierInterpolation(img, [factor,factor,...], 'lateral')` 等价。

    与默认零填充的 `fourier_zoom` 不同: 先对每帧做周期对称(镜像)延拓再傅里叶插值到
    奇数尺寸, 最后裁回 factor×原尺寸 —— 抑制非周期边界振铃。**RIED 原版 Fourier
    Interpolation 一步用的就是 'lateral'**, 逐模块对照显示: 不加它 python 与 matlab
    在 s2 起就出现 ~0.98 的相关差异并向下游放大 (s6 约 0.95); 加它则各阶段 ≈1.0。
    """
    factor = int(factor)
    if factor == 1:
        return np.asarray(img, dtype=np.float32)
    s = np.asarray(img, dtype=np.float32)
    single = (s.ndim == 2)                 # 鸭子: 2D 单帧 -> 返回 2D
    if single:
        s = s[None, ...]
    h, w = s.shape[-2], s.shape[-1]
    pre = (int(np.ceil(h / 2)), int(np.ceil(w / 2)))
    post = (int(np.floor(h / 2)), int(np.floor(w / 2)))
    # matlab: 每个横向轴 pad (ceil/2 前置, floor/2 后置) —— 对奇数尺寸 pre≠post!
    padw = list(zip(pre, post))                    # [(ceil(h/2),floor(h/2)), (ceil(w/2),floor(w/2))]
    ph, pw = h + pre[0] + post[0], w + pre[1] + post[1]
    newsz = (round(factor * ph - (factor - 1)), round(factor * pw - (factor - 1)))
    # matlab get_valid_part: sz = input_sz - ~mod(input_sz,2); idx = ceil(sz/2)+1+(f-1)floor(sz/2)
    sh = h - (1 - h % 2)
    sw = w - (1 - w % 2)
    iy = int(np.ceil(sh / 2) + 1 + (factor - 1) * np.floor(sh / 2)) - 1    # -> 0-based
    ix = int(np.ceil(sw / 2) + 1 + (factor - 1) * np.floor(sw / 2)) - 1

    batch = s.shape[:-2] or (1,)
    n_batch = int(np.prod(batch))
    flat = s.reshape(n_batch, h, w)
    out = np.empty((n_batch, factor * h, factor * w), dtype=np.float32)
    # 每帧估算: padded(float64 8B/px) + F/Fp/ifft(complex128 16B/px) + real(8B/px)
    # 注意: fourier_interp_2d 的 FFT 缓冲是在 newsz(=2*ph-1 量级) 尺寸上分配的,
    #       比 padded (ph*pw) 大 ~4 倍 —— 必须用 newsz 估算, 否则块太大仍爆内存。
    per_frame = 8 * ph * pw + 56 * newsz[0] * newsz[1]
    chunk = max(1, min(n_batch, int(5.0e8 / max(per_frame, 1))))
    for s0 in range(0, n_batch, chunk):
        sub = flat[s0:s0 + chunk].astype(np.float64)
        padded = np.pad(sub, [(0, 0), (pre[0], post[0]), (pre[1], post[1])],
                        mode="symmetric")
        ip = fourier_interp_2d(padded, newsz)
        out[s0:s0 + chunk] = ip[..., iy:iy + factor * h, ix:ix + factor * w]
    res = np.ascontiguousarray(out.reshape(batch + (factor * h, factor * w)),
                               dtype=np.float32)
    return res[0] if single else res


class FourierUpsample(Transform):
    """傅里叶上采样: ImageStack(T,H,W) -> (T, H*factor, W*factor)。

    同时按 factor 缩小 pixel_size (nm/pixel), 并累乘 metadata['upsample']。
    lateral=True 时用 matlab 等价的对称镜像 padding (fourier_lateral), 复现
    RIED 原版的 'lateral' 插值; 否则用默认零填充 (fourier_zoom)。
    加余弦边缘窗默认关闭 (taper=False), 避免压暗边缘。

    take_abs=True 时对插值结果取绝对值 (对齐 RIEDm.m 的
    `data_FI = abs(fourierInterpolation(...))`——lateral 插值的 sinc 过冲会
    产生 ~44% 负值, 缺失 abs() 会让负值进 RL3D → preRL 出现 2×2 棋盘伪影)。

    鸭子类型: 直接喂 ndarray 时返回 ndarray (不带 pixel_size 语义)。
    """
    def __init__(self, factor=2, taper=False, lateral=False, take_abs=False):
        self.factor = int(factor)
        self.taper = bool(taper)
        self.lateral = bool(lateral)
        self.take_abs = bool(take_abs)

    def _apply(self, data):
        if self.lateral:
            data = fourier_lateral(data, self.factor)
        else:
            data = fourier_upsample(data, self.factor, taper=self.taper)
        if self.take_abs:
            data = np.abs(data)
        return data

    def _wrap(self, stack, data):
        meta = dict(stack.metadata)
        meta['upsample'] = int(meta.get('upsample', 1)) * self.factor
        ps = stack.pixel_size
        new_ps = (ps / self.factor) if ps is not None else None
        return replace(stack, data=data, pixel_size=new_ps, metadata=meta)

    def __repr__(self):
        return (f"FourierUpsample(factor={self.factor}, taper={self.taper}, "
                f"lateral={self.lateral}, take_abs={self.take_abs})")
