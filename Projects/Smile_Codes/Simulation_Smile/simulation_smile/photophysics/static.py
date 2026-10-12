from __future__ import annotations
import numpy as np
from .base import PhotophysicsModel
from ..datatypes import PhotonStream


class Static(PhotophysicsModel):
    """最简单模型：不闪烁、不漂白，每帧发射量 = 亮度 × 激发 × dt，
    再叠加泊松 shot noise。"""
    def simulate(self, emitters, illumination, n_frames, dt, rng):
        exc = illumination.intensity(emitters.xyz)        # (N,)
        expected = emitters.brightness * exc * dt         # (N,) 每帧期望光子
        # 广播到 (N, T)，再泊松采样
        expected_frames = np.repeat(expected[:, None], n_frames, axis=1)
        counts = rng.poisson(expected_frames)             # (N, T) shot noise
        counts = counts[:, :, None]                       # -> (N, T, 1) 单通道
        return PhotonStream(counts=counts.astype(float))
