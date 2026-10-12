from __future__ import annotations
import numpy as np
from .base import PhotophysicsModel
from ..datatypes import PhotonStream


class Blinking(PhotophysicsModel):
    """由态转移速率矩阵 + 态亮度矩阵定义的闪烁/漂白模型。

    约定：
    - rate_matrix K: (S, S)，K[i, j] 是从态 i 到态 j 的转移速率，单位 1/s。
      对角线随意（会被忽略并重建为生成元 Q 的对角）。
    - 最后一个态（索引 S-1）定义为漂白态：强制为吸收态（无出射），亮度归零。
    - state_brightness: (S,) 或 (S, C)，每个态（每通道）的相对亮度因子。
      实际发射率 = emitter.brightness × 激发强度 × state_brightness[s, c]。

    帧内用 Gillespie 采样连续指数等待时间，积分每个态的停留时长，
    再对期望光子数施加泊松 shot noise。
    """

    def __init__(self, rate_matrix, state_brightness, initial_state=0):
        K = np.asarray(rate_matrix, dtype=float)
        assert K.shape[0] == K.shape[1], "速率矩阵必须为方阵"
        self.n_states = K.shape[0]
        self.bleached = self.n_states - 1

        # 构建生成元 Q：off-diagonal 为速率，最后一行清零（吸收态）
        offdiag = K.copy()
        np.fill_diagonal(offdiag, 0.0)
        offdiag[self.bleached, :] = 0.0          # 漂白态：无出射
        self.offdiag = offdiag                   # (S, S) 出射速率
        self.exit_rates = offdiag.sum(axis=1)    # (S,) 每个态的总出射速率

        # 态亮度矩阵 (S, C)，漂白态强制为 0
        b = np.asarray(state_brightness, dtype=float)
        if b.ndim == 1:
            b = b[:, None]
        b = b.copy()
        b[self.bleached] = 0.0
        self.state_brightness = b                # (S, C)
        self.n_channels = b.shape[1]

        self.initial_state = initial_state       # int 或 (S,) 概率分布

    def _draw_initial(self, rng) -> int:
        if np.isscalar(self.initial_state):
            return int(self.initial_state)
        return int(rng.choice(self.n_states, p=self.initial_state))

    def simulate(self, emitters, illumination, n_frames, dt, rng):
        N = emitters.n
        C = self.n_channels
        exc = illumination.intensity(emitters.xyz)      # (N,)

        counts = np.zeros((N, n_frames, C))
        states = np.zeros((N, n_frames), dtype=int)     # 记录每帧末状态

        for i in range(N):
            s = self._draw_initial(rng)
            base = emitters.brightness[i] * exc[i]       # 标量：emitter 的绝对亮度权重

            for t in range(n_frames):
                dwell = np.zeros(self.n_states)          # 本帧各态停留时长
                elapsed = 0.0

                # ---- 帧内 Gillespie ----
                while elapsed < dt:
                    rate = self.exit_rates[s]
                    if rate <= 0.0:                      # 吸收态/无出射：停到帧末
                        dwell[s] += dt - elapsed
                        break
                    tau = rng.exponential(1.0 / rate)    # 连续指数等待时间
                    if elapsed + tau >= dt:              # 下次跳变超出本帧
                        dwell[s] += dt - elapsed
                        break
                    dwell[s] += tau
                    elapsed += tau
                    # 按出射速率比例选择下一个态（cumsum + 采样，避免归一化）
                    r = rng.random() * rate
                    s = int(np.searchsorted(np.cumsum(self.offdiag[s]), r))

                # ---- 帧内积分 → 期望光子 → 泊松噪声 ----
                expected = base * (dwell @ self.state_brightness)   # (C,)
                counts[i, t, :] = rng.poisson(expected)
                states[i, t] = s

        return PhotonStream(counts=counts, states=states)
