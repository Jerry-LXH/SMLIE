from __future__ import annotations
import numpy as np
from .base import PSFModel


class GaussianPSF(PSFModel):
    """各向同性二维高斯 PSF, 局部窗口渲染 + 核缓存。

    相比"每个分子铺整张全视场核"的旧实现:
      1. 每个分子只在以它为中心的截断窗口 (±trunc·σ) 内铺核,
         复杂度从 O(N·H·W) 降到 O(N·Ww²) (本例 ~800 倍);
      2. 分子位置静态时, 各分子核只算一次并缓存, 之后每帧只做加权累加。

    截断损失: 半径 trunc·σ 外的光子占比 ≈ exp(-trunc²/2);
    trunc=4 时约 0.03%, 远低于光子 shot noise, 可忽略。
    注意: 若在外部原地修改 emitters.xyz (如引入漂移), 缓存不会自动失效,
    请通过 clear_cache() 重置。
    """
    def __init__(self, sigma_nm: float, n_channels: int = 1, trunc: float = 4.0):
        self.sigma = sigma_nm
        self.n_channels = n_channels
        self.trunc = float(trunc)
        self._cache = {}             # id(emitters) -> (xyz, kernels, slices)

    def clear_cache(self):
        self._cache.clear()

    def _build_kernels(self, emitters, fov):
        """为每个 emitter 预计算局部高斯核及其在整图中的像素切片。

        返回 (kernels, slices):
          kernels[i] : (Ww, Ww) 相对窗口内的高斯核 (未归一化之外已含 norm);
          slices[i]  : (r0, r1, c0, c1) 窗口在整图中的 [r0:r1, c0:c1] 像素索引。
        """
        px = fov.pixel_size
        sigma_px = self.sigma / px
        r = max(int(np.ceil(self.trunc * sigma_px)), 2)
        s2 = 2 * self.sigma ** 2
        norm = 1.0 / (2 * np.pi * self.sigma ** 2) * px ** 2

        kernels, slices = [], []
        for i in range(emitters.n):
            x0, y0 = emitters.xyz[i, 0], emitters.xyz[i, 1]
            c = int(np.round(x0 / px - 0.5))     # 峰值所在像素列索引 (0-based 中心)
            rw = int(np.round(y0 / px - 0.5))    # 峰值所在像素行索引
            c0, c1 = c - r, c + r + 1
            r0, r1 = rw - r, rw + r + 1
            # 窗口内像素中心相对 emitter 的位移 (nm)
            Xc = (np.arange(c0, c1) + 0.5) * px - x0   # (Ww,)
            Yc = (np.arange(r0, r1) + 0.5) * px - y0   # (Ww,)
            k = norm * np.exp(-(Xc[None, :] ** 2 + Yc[:, None] ** 2) / s2)
            kernels.append(k)
            slices.append((r0, r1, c0, c1))
        return kernels, slices

    def render(self, emitters, photon_frame, fov, rng=None):
        """把该帧光子按高斯 PSF 多项分布随机分配到像面, 返回 (C, H, W) 计数图。

        对每个 emitter/通道的 N 个光子, 以截断窗口内归一化高斯核为类别概率做
        Multinomial(N, p) 分配(逐光子落点), 得到的光子直方图即 shot noise 后的
        像面; 取代旧的"期望值加权累加"。rng=None 时内部新建 default_rng()。
        """
        if rng is None:
            rng = np.random.default_rng()
        img = np.zeros((self.n_channels, fov.height, fov.width))
        if emitters.n == 0:
            return img

        cached = self._cache.get(id(emitters))
        if cached is None or cached[0] is not emitters.xyz:
            cached = (emitters.xyz, *self._build_kernels(emitters, fov))
            self._cache[id(emitters)] = cached
        _, kernels, slices = cached

        H, W = fov.height, fov.width
        for i in range(emitters.n):
            r0, r1, c0, c1 = slices[i]
            R0, C0 = max(r0, 0), max(c0, 0)        # 裁剪到图像范围内
            R1, C1 = min(r1, H), min(c1, W)
            if R1 <= R0 or C1 <= C0:
                continue
            k = kernels[i][R0 - r0:R1 - r0, C0 - c0:C1 - c0]
            s = k.sum()
            if s <= 0:
                continue
            p = k / s                              # 可见窗口内再归一化(能量守恒)
            pflat = p.ravel()
            for c in range(self.n_channels):
                n = int(round(photon_frame[i, c]))
                if n <= 0:
                    continue
                counts = rng.multinomial(n, pflat).reshape(k.shape)
                img[c, R0:R1, C0:C1] += counts
        return img
