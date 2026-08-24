from __future__ import annotations
import numpy as np
from .base import EmitterSource, sample_size
from ..datatypes import EmitterSet


class Random3D(EmitterSource):
    """在玻底上方 z ∈ [0, z_max] 范围内均匀随机撒点 (底部固定 z=0)。

    - xy 与 Random2D 相同: 均匀分布, 数量由 density_per_um2 (面密度) 决定;
    - z 在 [0, z_max] 内均匀采样 (nm), 用于 3D / 离焦 / 像散等轴向依赖模拟。
    size: 同 Random2D, (mean_nm, variance_nm²) 或 None。
    """
    def __init__(self, fov_size_nm: tuple[float, float],
                 density_per_um2: float, brightness: float = 5000.0,
                 z_max: float = 500.0, size=None):
        self.fov_size = np.asarray(fov_size_nm, dtype=float)  # (W_nm, H_nm)
        self.density = density_per_um2
        self.brightness = brightness
        self.z_max = float(z_max)                    # 轴向范围上限 (nm), 底部为 0
        self.size = None if size is None else (float(size[0]), float(size[1]))

    def sample(self, rng):
        area_um2 = np.prod(self.fov_size) / 1e6      # nm^2 -> um^2
        n = rng.poisson(self.density * area_um2)
        xy = rng.uniform([0, 0], self.fov_size, size=(n, 2))
        z = rng.uniform(0.0, self.z_max, size=n)     # z ∈ [0, z_max]
        xyz = np.column_stack([xy, z])
        brightness = np.full(n, self.brightness)
        sizes = sample_size(rng, n, self.size)
        return EmitterSet(xyz=xyz, brightness=brightness, size=sizes)
