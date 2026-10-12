from __future__ import annotations
import numpy as np
from .simple_cmos import SimpleCMOS


class SimpleEMCCD(SimpleCMOS):
    """极简 EMCCD：泊松背景 + EM 倍增(随机 gamma) + 读出噪声 + ADU 量化。

    与 SimpleCMOS 的结构对齐, 差异只在"电子倍增": 输入电子 N 先经 EM 寄存器
    随机倍增, 输出放大器读出噪声在倍增**之后**叠加 (不被 EM 放大), 再按
    sensitivity (e-/ADU) 转数字值。背景/不均匀背景逻辑直接继承 SimpleCMOS。

    物理模型 (每像素每帧):
      1) 输入电子数 N = img + 泊松背景 (bg_noise_e, 可空间不均匀)
         —— 背景电子同样经过 EM 寄存器, 所以并入 N 一起倍增;
      2) EM 倍增: M ~ Gamma(shape=N, scale=em_gain)
         —— EMCCD 标准模型: 每个输入电子独立得到指数分布增益(均值 em_gain),
            N 个输入之和 = Gamma(N, em_gain), 内禀含 excess noise factor √2;
            N=0 (无输入电子) 时 M=0;
      3) 读出噪声: M += N(0, read_noise_e)
         —— 输出放大器噪声在 EM 之后叠加, 不被倍增; 若放倍增前会被放大
            em_gain² 倍而主导总噪声 (物理错误);
      4) ADU = round(M / sensitivity + offset)   (sensitivity 单位 e-/ADU)。

    参数:
      em_gain              : EM 增益 (平均倍增倍数, 默认 300);
      sensitivity          : 转换系数 (e-/ADU, 默认 5.75);
      offset               : 偏置 (ADU, 默认 486);
      bg_noise_e           : 背景均值 (e-/px, 默认 5.0);
      read_noise_e         : 读出噪声标准差 (e-, 倍增后叠加, 默认 1.0);
      bg_inhomogeneous     : 同 SimpleCMOS (空间不均匀背景);
      bg_inhomogeneity_frac: 同 SimpleCMOS。
    """
    def __init__(self, em_gain: float = 300.0, sensitivity: float = 5.75,
                 offset: float = 486.0, bg_noise_e: float = 5.0,
                 read_noise_e: float = 1.0,
                 bg_inhomogeneous: bool = False,
                 bg_inhomogeneity_frac: float = 0.3):
        # 复用 SimpleCMOS 字段; 其 gain (ADU/e-) 在 EMCCD 中不用, 占位保持结构一致
        super().__init__(gain=1.0, bg_noise_e=bg_noise_e,
                         read_noise_e=read_noise_e, offset=offset,
                         bg_inhomogeneous=bg_inhomogeneous,
                         bg_inhomogeneity_frac=bg_inhomogeneity_frac)
        self.em_gain = em_gain
        self.sensitivity = sensitivity

    def capture(self, photon_image, rng):
        # ---- 1) 输入电子数 N = 信号 + 泊松背景 (同 SimpleCMOS) ----
        img = photon_image.astype(float)
        if self.bg_inhomogeneous:
            # 不均匀背景: 固定空间模式(跨帧复用), 每帧只叠加泊松 shot noise
            if self._bg_map is None or self._bg_map.shape != img.shape:
                self._bg_map = self._make_bg_map(img.shape, rng)
            img = img + rng.poisson(self._bg_map)
        else:
            img = img + rng.poisson(self.bg_noise, size=img.shape)   # 均匀背景

        # ---- 2) EM 倍增: M ~ Gamma(shape=N, scale=em_gain) ----
        N = np.maximum(img, 0.0)                 # 输入电子数 (>=0)
        M = rng.gamma(N, self.em_gain)           # 均值 = N * em_gain

        # ---- 3) 读出噪声: EM 之后叠加 (不被倍增) ----
        M = M + rng.normal(0.0, self.read_noise, size=M.shape)

        # ---- 4) ADU 转换 ----
        adu = M / self.sensitivity + self.offset
        return np.round(np.maximum(adu, 0.0)).astype(np.uint16)
