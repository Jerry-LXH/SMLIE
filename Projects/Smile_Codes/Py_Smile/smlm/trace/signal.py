"""信号阶段的实现: 从一个裁切块里取出该分子的原始积分强度。

``raw`` 的口径
-------------
``raw = Σ w·I / Σ w²``, 即**已经按权重平方和归一化**的加权和。之所以不存原始
``Σ w·I``, 是为了让 ``intensity = raw − bg·sw/sw2`` 可以直接由 ``raw`` 与 ``bg``
算出, 不必回头看图像。

量纲: 只有 ``mode='hard'`` 时 ``raw`` 才是孔径内的光子数; ``mode='matched'`` 下
权重核归一化到 Σg=1, 所以 ``raw ≈ 峰值 A · 2πσ²`` (a.u.), 下游的比值与 SNR
都在同一口径内自洽, 但不要把它当光子数去和 MLE 的 photons 直接比较。

因此这个阶段同时写 ``sw`` (=Σw) 与 ``sw2`` (=Σw²): 下游
:class:`smlm.trace.calculate.IntensityFromRawBg` 与
:class:`smlm.trace.calculate.PoissonSNR` 都要用它们做换算与噪声估计。
换一种 raw 估计方式 (高斯拟合、质心加权……) 时, 只要把这三个量按同一口径写出来,
下游一行都不用改。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .aperture import aperture_kernel, radius_from_sigma
from .fill import TraceFiller


@dataclass(frozen=True)
class ApertureRaw(TraceFiller):
    """孔径加权的原始强度。

    Parameters
    ----------
    mode : 'matched' | 'hard'
        ``matched`` 匹配滤波 (默认, 收集效率与 σ 无关); ``hard`` 硬圆盘。
    k_sigma : float
        半径 = ``k_sigma·σ``(逐分子)。
    radius : int, optional
        显式固定半径, 给了就忽略 ``k_sigma``。
    min_radius, max_radius : int
        半径上下限。
    """

    name = 'aperture_raw'
    provides = ('sigma', 'radius', 'mode', 'sw', 'sw2', 'area', 'raw')
    requires = ()

    mode: str = 'matched'
    k_sigma: float = 2.0
    radius: int | None = None
    min_radius: int = 2
    max_radius: int = 7

    def __post_init__(self):
        if self.mode not in ('matched', 'hard'):
            raise ValueError(f"mode 必须是 'matched' 或 'hard', 得到 {self.mode!r}")

    def plan(self, ctx) -> int:
        ctx.radius = radius_from_sigma(ctx.sigma, k_sigma=self.k_sigma,
                                       radius=self.radius,
                                       min_radius=self.min_radius,
                                       max_radius=self.max_radius)
        return ctx.radius

    def fill(self, trace, ctx, box) -> None:
        w, sw, sw2, area, sigma = aperture_kernel(
            self.mode, ctx.radius, ctx.sigma, ctx.center_frac)
        trace['sigma'] = float(sigma)
        trace['radius'] = int(ctx.radius)
        trace['mode'] = self.mode
        trace['sw'] = float(sw)
        trace['sw2'] = float(sw2)
        trace['area'] = float(area)
        trace['raw'] = np.einsum('tij,ij->t', box.signal(ctx.radius), w) / sw2

    def describe(self) -> str:
        if self.radius is not None:
            r = f'r={self.radius}(固定)'
        else:
            r = f'r={self.k_sigma:g}σ'
        return f'{self.mode}, {r}, [{self.min_radius}, {self.max_radius}] px'
