from __future__ import annotations
import numpy as np
from .base import PhotophysicsModel
from ..datatypes import PhotonStream


class Ecl(PhotophysicsModel):
    """ECL 复合泊松发射模型（burst 打包的光子发射）。

    每帧、每个 emitter：
      1) 先从泊松分布采样"发射事件数" M；
      2) 每个事件再从 burst 大小分布（1..10 个光子）采样 B_i；
      3) 本帧发射的光子总数
             N = sum_{i=1}^{M} B_i     （M = 0 时 N = 0，即 B0 = 0）

    参数
    ----
    burst_probs : (10,) array_like
        burst 大小 1..10 的概率：第 0 项对应 1 个光子、…、第 9 项对应 10 个光子。
        自动归一化到和为 1。

    泊松速率约定：为保持 ``EmitterSet.brightness`` 的"光子/s"含义，
        lambda_M = brightness * excitation * dt / E[B]
    其中 ``E[B] = sum(b * p_b)`` 是平均 burst 大小；因此每帧平均总光子数
        E[M] * E[B] = brightness * excitation * dt。

    该模型只改变发射的离散性（burst 打包），不含 on/off 闪烁或漂白。
    """

    MAX_BURST = 10

    def __init__(self, burst_probs):
        p = np.asarray(burst_probs, dtype=float).ravel()
        assert p.size == self.MAX_BURST, \
            f"burst_probs 长度必须为 {self.MAX_BURST}，实际 {p.size}"
        p = np.clip(p, 0.0, None)
        s = p.sum()
        assert s > 0, "burst_probs 不能全为 0"
        self.burst_probs = p / s
        self.burst_sizes = np.arange(1, self.MAX_BURST + 1, dtype=float)
        self.mean_burst = float((self.burst_sizes * self.burst_probs).sum())

    def simulate(self, emitters, illumination, n_frames, dt, rng):
        exc = illumination.intensity(emitters.xyz)          # (N,) 激发强度
        mean_photons = emitters.brightness * exc * dt       # (N,) 每帧期望总光子
        rate = mean_photons / self.mean_burst               # (N,) 每帧期望 burst 数
        rate_frames = np.repeat(rate[:, None], n_frames, axis=1)   # (N, T)

        M = rng.poisson(rate_frames)                        # (N, T) burst 事件数
        n_b = rng.multinomial(M, self.burst_probs)          # (N, T, 10) 各大小 burst 数
        counts = (n_b * self.burst_sizes).sum(axis=-1)      # (N, T) 每帧总光子数
        counts = counts[:, :, None]                         # (N, T, 1) 单通道
        return PhotonStream(counts=counts.astype(float))
