from __future__ import annotations
import numpy as np
from .base import EmitterSource, sample_size
from ..datatypes import EmitterSet


class SimpleRing2d(EmitterSource):
    """在视场中心生成居中的环形 (内圈/外圈) 发射体, 全部落在 z=0。

    环形定义 (nm):
      - radius   : 环的中心线半径 (即环心到圆环中线的距离);
      - thickness: 环壁的径向总厚度;
      因此内圈半径 r_inner = radius - thickness/2, 外圈半径 r_outer = radius + thickness/2,
      r_inner 会被截断到 >= 0。

    数量: 由线密度 (数量/nm) 乘环中心线周长决定,
        n ~ Poisson(line_density_per_nm * 2π * radius)

    分布: 在环上随机均匀分布 —— 角度 θ ~ U(0, 2π), 径向 r ~ U(r_inner, r_outer),
    因此每个发射体都以等角度间隔出现在环上, 线密度 ~ 数量/nm 近似恒定。

    size: 同 Random2D, (mean_nm, variance_nm²) 或 None。
    """
    def __init__(self, fov_size_nm: tuple[float, float],
                 radius: float, thickness: float,
                 line_density_per_nm: float, brightness: float = 5000.0,
                 size=None):
        self.fov_size = np.asarray(fov_size_nm, dtype=float)  # (W_nm, H_nm)
        self.radius = float(radius)
        self.thickness = float(thickness)
        self.line_density = line_density_per_nm               # 数量/nm
        self.brightness = brightness
        self.size = None if size is None else (float(size[0]), float(size[1]))

    def sample(self, rng):
        center = self.fov_size / 2.0                          # 视场中心
        r_inner = max(self.radius - self.thickness / 2.0, 0.0)
        r_outer = self.radius + self.thickness / 2.0
        circumference = 2 * np.pi * self.radius               # 中心线周长 nm
        n = rng.poisson(self.line_density * circumference)

        theta = rng.uniform(0.0, 2 * np.pi, size=n)           # 角度均匀
        r = rng.uniform(r_inner, r_outer, size=n)             # 径向均匀
        x = center[0] + r * np.cos(theta)
        y = center[1] + r * np.sin(theta)
        xyz = np.column_stack([x, y, np.zeros(n)])            # z = 0
        brightness = np.full(n, self.brightness)
        sizes = sample_size(rng, n, self.size)
        return EmitterSet(xyz=xyz, brightness=brightness, size=sizes)
