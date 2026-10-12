from __future__ import annotations
import numpy as np
from .base import Detector


class SimpleCMOS(Detector):
    """极简 sCMOS：泊松背景 + 高斯读出噪声 + 增益量化。

    输入光子图已是 optics 按 PSF 多项分配后的像面计数(信号 shot noise 已含),
    因此这里只做: 叠加泊松背景 -> 高斯读出噪声 -> 增益 + 偏置 -> 量化。
    注: qe 已去除——可等效合并到发射速率里, 留到 Phase 2 再做。

    参数:
      gain                 : 增益 (ADU/e-);
      bg_noise_e           : 背景均值 (e-/px);
      read_noise_e         : 读出噪声标准差 (e-);
      offset               : 偏置 (ADU);
      bg_inhomogeneous     : True=生成"空间不均匀"背景 (随机缓慢变化, 均值=bg_noise_e,
                             模式跨帧固定, 每帧只叠加泊松 shot noise); False=均匀背景(原行为);
      bg_inhomogeneity_frac: 不均匀幅度 = frac × bg_noise_e (仅 bg_inhomogeneous=True 时生效)。
    """
    def __init__(self, gain: float = 1.0, bg_noise_e: float = 5.0,
                 read_noise_e: float = 1.5, offset: float = 100.0,
                 bg_inhomogeneous: bool = False,
                 bg_inhomogeneity_frac: float = 0.3):
        self.gain = gain
        self.bg_noise = bg_noise_e
        self.read_noise = read_noise_e
        self.offset = offset
        self.bg_inhomogeneous = bool(bg_inhomogeneous)
        self.bg_inhomogeneity_frac = float(bg_inhomogeneity_frac)
        self._bg_map = None                 # 缓存不均匀背景模式 (跨帧固定)

    @staticmethod
    def _resize_bilinear(coarse, h, w):
        """把粗网格双线性上采样到 (h, w) —— 纯 numpy, 零 scipy 依赖。"""
        ch, cw = coarse.shape
        yy = (np.arange(h) + 0.5) * ch / h - 0.5
        xx = (np.arange(w) + 0.5) * cw / w - 0.5
        y0 = np.clip(np.floor(yy).astype(int), 0, ch - 1)
        x0 = np.clip(np.floor(xx).astype(int), 0, cw - 1)
        y1 = np.minimum(y0 + 1, ch - 1)
        x1 = np.minimum(x0 + 1, cw - 1)
        fy = (yy - y0)[:, None]
        fx = (xx - x0)[None, :]
        c00 = coarse[np.ix_(y0, x0)]
        c01 = coarse[np.ix_(y0, x1)]
        c10 = coarse[np.ix_(y1, x0)]
        c11 = coarse[np.ix_(y1, x1)]
        top = c00 * (1 - fx) + c01 * fx
        bot = c10 * (1 - fx) + c11 * fx
        return top * (1 - fy) + bot * fy

    def _make_bg_map(self, shape, rng):
        """随机缓慢变化背景均值场: 均值=bg_noise, 幅度=frac×bg_noise, 各通道共享。

        粗网格(约 1/16 分辨率)白噪声 -> 双线性上采样到全分辨率,
        得低空间频率的平滑背景 (空间上缓慢变化, 无高频跳变)。
        """
        h, w = shape[-2:]
        ch = max(int(round(h / 16)), 2)
        cw = max(int(round(w / 16)), 2)
        coarse = rng.normal(size=(ch, cw))
        field = self._resize_bilinear(coarse, h, w)
        field = field - field.mean()
        amp = max(self.bg_noise * self.bg_inhomogeneity_frac, 1e-9)
        bg = self.bg_noise + field / (field.std() + 1e-12) * amp
        bg = np.maximum(bg, 0.0)
        if len(shape) == 3:                       # (C,H,W): 各通道共享同一背景模式
            bg = np.broadcast_to(bg[None], shape)
        return np.ascontiguousarray(bg)

    def capture(self, photon_image, rng):
        img = photon_image.astype(float)
        if self.bg_inhomogeneous:
            # 不均匀背景: 固定空间模式(跨帧复用), 每帧只叠加泊松 shot noise
            if self._bg_map is None or self._bg_map.shape != img.shape:
                self._bg_map = self._make_bg_map(img.shape, rng)
            img = img + rng.poisson(self._bg_map)          # 逐像素均值场 → 泊松
        else:
            img = img + rng.poisson(self.bg_noise, size=img.shape)   # 均匀背景 (原行为)
        img = img + rng.normal(0.0, self.read_noise, size=img.shape)   # 读出噪声
        adu = img * self.gain + self.offset
        return np.round(adu).astype(np.uint16)
