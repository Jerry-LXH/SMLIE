from __future__ import annotations
import numpy as np
from .base import EmitterSource, sample_size
from ..datatypes import EmitterSet


class Random2D(EmitterSource):
    """在玻底 (z=0) 均匀随机撒点, 密度可调。

    size: (mean_nm, variance_nm²) 或 None; 提供时每个粒子的物理尺寸由
    高斯随机采样得到(截断到非负), 不提供时 EmitterSet.size = None。
    """
    def __init__(self, fov_size_nm: tuple[float, float],
                 density_per_um2: float, brightness: float = 5000.0,
                 size=None):
        self.fov_size = np.asarray(fov_size_nm, dtype=float)  # (W_nm, H_nm)
        self.density = density_per_um2
        self.brightness = brightness
        self.size = None if size is None else (float(size[0]), float(size[1]))

    def sample(self, rng):
        area_um2 = np.prod(self.fov_size) / 1e6   # nm^2 -> um^2
        n = rng.poisson(self.density * area_um2)
        xy = rng.uniform([0, 0], self.fov_size, size=(n, 2))
        xyz = np.column_stack([xy, np.zeros(n)])
        brightness = np.full(n, self.brightness)
        sizes = sample_size(rng, n, self.size)
        return EmitterSet(xyz=xyz, brightness=brightness, size=sizes)
