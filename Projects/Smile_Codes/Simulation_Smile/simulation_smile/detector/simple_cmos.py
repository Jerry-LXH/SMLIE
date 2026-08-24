from __future__ import annotations
import numpy as np
from .base import Detector


class SimpleCMOS(Detector):
    """极简 sCMOS：背景 + 泊松 shot noise + 高斯读出噪声 + 增益量化。
    （EMCCD 的乘性增益噪声、逐像素噪声留到 Phase 2。）"""
    def __init__(self, gain: float = 1.0, bg_noise_e: float = 5.0, read_noise_e: float = 1.5,
                 offset: float = 100.0, quantum_eff: float = 0.9):
        self.gain = gain
        self.bg_noise = bg_noise_e
        self.read_noise = read_noise_e
        self.offset = offset
        self.qe = quantum_eff

    def capture(self, photon_image, rng):
        electrons = rng.poisson((photon_image + self.bg_noise) * self.qe)   # bg + 光电转换 + shot
        electrons = electrons + rng.normal(0, self.read_noise,
                                           size=electrons.shape)
        adu = electrons * self.gain + self.offset
        return np.round(adu).astype(np.uint16)
