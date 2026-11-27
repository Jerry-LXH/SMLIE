"""由 ``raw`` 与 ``bg`` 推出 ``intensity`` / ``snr``。

这两个阶段只做算术, 不碰图像; 因此把它们和信号/背景阶段分开, 换信号或换背景都
不用动这里。

口径
----
``raw`` 已按 ``Σw²`` 归一化 (见 :mod:`smlm.trace.signal`), 所以:

    bg_flux   = bg · sw / sw2          # 孔径里的背景光子数
    intensity = raw − bg_flux          # 扣背景后的信号

SNR 用"孔径内总光子"的 Poisson 噪声加背景的方差::

    total    = intensity · sw2
    area_eff = sw² / sw2               # 等效像素数
    noise²   = total + area_eff · bg_std²
    snr      = total / noise

注意 ``bg_std`` 是**标准差** (与 ``bg`` 同量纲), 所以这里要平方; 只有纯 Poisson
且增益为 1 时才有 ``bg_std² ≈ bg``。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fill import TraceFiller


@dataclass(frozen=True)
class IntensityFromRawBg(TraceFiller):
    """``intensity = raw − bg·sw/sw2``, 同时给出 ``bg_flux = bg·sw/sw2``。"""

    name = 'intensity'
    provides = ('intensity', 'bg_flux')
    requires = ('raw', 'bg', 'sw', 'sw2')

    def fill(self, trace, ctx, box) -> None:
        sw = float(trace['sw'])
        sw2 = float(trace['sw2'])
        bg_flux = np.asarray(trace['bg'], float) * (sw / sw2)
        trace['bg_flux'] = bg_flux
        trace['intensity'] = np.asarray(trace['raw'], float) - bg_flux

    def describe(self) -> str:
        return 'intensity = raw - bg*sw/sw2'


@dataclass(frozen=True)
class PoissonSNR(TraceFiller):
    """孔径总光子 / (Poisson + 背景方差)^½, 非正信号处记 0。"""

    name = 'poisson_snr'
    provides = ('snr',)
    requires = ('intensity', 'bg_std', 'sw', 'sw2', 'area')

    def fill(self, trace, ctx, box) -> None:
        sw = float(trace['sw'])
        sw2 = float(trace['sw2'])
        # 等效像素数: 权重核越集中 (σ 越小), area_eff 越接近 1
        area_eff = float(sw ** 2 / sw2) if sw2 > 0 else float(trace['area'])

        total = np.asarray(trace['intensity'], float) * sw2
        bg_std = np.asarray(trace['bg_std'], float)
        noise = np.sqrt(np.maximum(total, 0.0) + area_eff * bg_std ** 2)
        trace['snr'] = np.where((total > 0) & (bg_std > 0),
                                total / np.maximum(noise, 1e-12), 0.0)

    def describe(self) -> str:
        return 'snr = total / sqrt(total + area_eff*bg_std^2)'


def standard_fillers(spec=None):
    """把 ``ApertureSpec`` 翻译成一条完整的阶段链 (信号 → 背景 → 强度 → SNR)。

    只覆盖 ``bg_ref='frames'``、``bg_smooth=0`` 的情形 —— 这是
    :class:`smlm.trace.annulus.AnnulusRing` 目前实现的背景; 其余口径请自己组合阶段。
    """
    from .aperture import ApertureSpec
    from .annulus import AnnulusRing
    from .signal import ApertureRaw

    spec = spec if isinstance(spec, ApertureSpec) else ApertureSpec(**(spec or {}))
    if spec.bg_ref != 'frames' or int(spec.bg_smooth) > 1:
        raise ValueError("standard_fillers 只支持逐帧背景 (bg_ref='frames', "
                         f"bg_smooth<=1); 得到 bg_ref={spec.bg_ref!r}, "
                         f"bg_smooth={spec.bg_smooth}")
    return [ApertureRaw(mode=spec.mode, k_sigma=spec.k_sigma, radius=spec.radius,
                        min_radius=spec.min_radius, max_radius=spec.max_radius),
            AnnulusRing(gap=spec.bg_gap, width=spec.bg_width,
                        excl_radius=spec.excl_radius, clip_n=spec.clip_n,
                        clip_iter=spec.clip_iter),
            IntensityFromRawBg(), PoissonSNR()]
