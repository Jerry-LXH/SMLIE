"""背景阶段: 在分子周围一个甜甜圈 (annulus) 里逐帧估计背景。

和 :class:`smlm.trace.signal.ApertureRaw` 对称 —— 一个负责"信号框里有多少光",
一个负责"信号框外一圈有多少背景"。

为什么环要避开信号框
-------------------
环的内边从 ``半径 + gap`` 起算, 保证不碰到 PSF 的裙边; 外边界是
``半径 + gap + width``。环里**其他分子**的核心会被挖掉 (``excl_radius``),
否则它们的外周会把背景抬高, 于是 trace 偏低、E 偏高。

输出口径
--------
``bg`` 是**每像素**背景 (ADU/px), 逐帧一个值; ``bg_std`` 是同一批环像素的逐帧
**标准差** —— 与 ``bg`` 同量纲, 但不是方差, 因此通常远小于 ``bg``。只有纯 Poisson
且增益为 1 时才有 ``bg_std² ≈ bg``; 真实 EMCCD 数据一般是 ``bg_std ≈ √bg`` (以光子计)。
SNR 用的是 ``bg_std²`` (方差), 见 :class:`smlm.trace.calculate.PoissonSNR`。
要把它们换成"孔径里的背景光子数"是
:class:`smlm.trace.calculate.IntensityFromRawBg` 的事, 不在这里做。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .aperture import ApertureSpec
from .background import AnnulusBackground
from .fill import TraceFiller


@dataclass(frozen=True)
class AnnulusRing(TraceFiller):
    """逐帧 sigma-clip 中位数, 估计每个分子周围的背景水平。

    单位是 **ADU/像素**(环里每个像素的背景), 不是孔径内的总量 —— 折算进信号核
    是下一个阶段的事: :class:`~smlm.trace.calculate.IntensityFromRawBg` 会写
    ``bg_flux = bg·sw/sw2``, 那个才与 ``raw`` / ``intensity`` 同单位
    (``intensity = raw − bg_flux``)。画图要和信号放在一起看时用 ``bg_flux``。

    Parameters
    ----------
    gap : int
        信号框外让出多少 px 才开始算背景。
    width : int
        环的宽度。
    excl_radius : int
        环里每个邻居分子被挖掉的半径 (只挖核心, 外周仍进估计)。设为 0 表示不挖 ——
        此时 ``neighbours`` 也记 0 (没有核心被移除)。
    bg_ref : 'frames' | 'tail'
        ``frames`` 逐帧估计(默认, 跟得住时间漂移); ``tail`` 只取最后 ``bg_tail``
        帧聚合成常数(噪声最低, 但要求背景不随时间变)。
    bg_smooth : int
        逐帧背景的**居中滑动平均**窗宽(帧), 0 = 不平滑。压掉环中位数的帧间抖动
        (√w 倍), 代价是轻微抹平真实的快变。
    bg_tail : int
        ``bg_ref='tail'`` 时参与聚合的帧数。
    clip_n, clip_iter : float, int
        逐帧 sigma-clip 的倍数与迭代次数。
    """

    name = 'annulus'
    provides = ('bg', 'bg_std', 'n_bg_pixels', 'neighbours')
    requires = ('radius',)          # 环的大小由信号半径决定

    gap: int = 2
    width: int = 4
    excl_radius: int = 2
    bg_ref: str = 'frames'
    bg_smooth: int = 0
    bg_tail: int = 200
    clip_n: float = 3.0
    clip_iter: int = 3

    def __post_init__(self):
        if int(self.gap) < 0:
            raise ValueError(f'gap 不能为负, 得到 {self.gap}')
        if int(self.width) <= 0:
            raise ValueError(f'width 必须为正, 得到 {self.width}')
        if int(self.excl_radius) < 0:
            raise ValueError(f'excl_radius 不能为负, 得到 {self.excl_radius}')
        if int(self.bg_smooth) < 0:
            raise ValueError(f'bg_smooth 不能为负, 得到 {self.bg_smooth}')

    @property
    def spec(self) -> ApertureSpec:
        """低层 :class:`AnnulusBackground` 仍按 spec 取环参数, 这里只填背景那一半。"""
        return ApertureSpec(bg_gap=int(self.gap), bg_width=int(self.width),
                            excl_radius=int(self.excl_radius),
                            bg_ref=str(self.bg_ref),
                            bg_smooth=int(self.bg_smooth),
                            bg_tail=int(self.bg_tail),
                            clip_n=float(self.clip_n),
                            clip_iter=int(self.clip_iter))

    def plan(self, ctx) -> int:
        if not ctx.radius:
            raise ValueError('AnnulusRing 需要一个先算好的信号半径; '
                             '请把 ApertureRaw 放在它前面')
        return int(ctx.radius) + int(self.gap) + int(self.width)

    def fill(self, trace, ctx, box) -> None:
        half = box.half
        nrel, _ = ctx.neighbours(half, self.excl_radius)
        ring = AnnulusBackground(self.spec, ctx.radius, half,
                                 neighbours_rel=nrel, valid=box.valid(),
                                 center_frac=ctx.center_frac)
        bg, info = ring.estimate(box.data)
        trace['bg'] = np.asarray(bg, float)
        trace['bg_std'] = np.asarray(info['bg_std'], float)
        trace['n_bg_pixels'] = int(info['n_pixels'])
        # 只数"核心真的被挖掉"的那些; excl_radius=0 时一个都没挖
        removed = nrel is not None and int(self.excl_radius) > 0
        trace['neighbours'] = int(len(nrel)) if removed else 0

    def describe(self) -> str:
        extra = ''
        if int(self.bg_smooth) > 0:
            extra = f', smooth {self.bg_smooth}f'
        if self.bg_ref != 'frames':
            extra += f', ref {self.bg_ref}'
        return (f'ring +{self.gap}..+{self.gap + self.width} px, '
                f'excl {self.excl_radius}, clip {self.clip_n:g}x{self.clip_iter}'
                f'{extra}')
