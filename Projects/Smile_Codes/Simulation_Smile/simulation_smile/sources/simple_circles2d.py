from __future__ import annotations
import numpy as np
from .base import EmitterSource, sample_size
from ..datatypes import EmitterSet


class SimpleCircles2d(EmitterSource):
    """在视场中心生成一组同心圆上的发射体, 全部落在 z=0。

    radii: 每个圆的半径 (nm) 组成的数组, 其长度即圆的个数。
    每个圆上独立采样:
        n_i ~ Poisson(line_density_per_nm * 2π * r_i)
        θ   ~ U(0, 2π)
    因此每个圆上的线密度 ~ 数量/nm 近似恒定, 发射体恰好落在圆周线上
    (与 SimpleRing2d 在环带内分布不同)。

    size: 同 Random2D, (mean_nm, variance_nm²) 或 None。
    """
    def __init__(self, fov_size_nm: tuple[float, float],
                 radii, line_density_per_nm: float,
                 brightness: float = 5000.0, size=None):
        self.fov_size = np.asarray(fov_size_nm, dtype=float)  # (W_nm, H_nm)
        self.radii = np.asarray(radii, dtype=float)           # (n_circles,) nm
        self.line_density = line_density_per_nm               # 数量/nm
        self.brightness = brightness
        self.size = None if size is None else (float(size[0]), float(size[1]))

    def sample(self, rng):
        center = self.fov_size / 2.0                          # 视场中心
        xy_parts = []
        for r in self.radii:
            n = rng.poisson(self.line_density * 2 * np.pi * r)
            theta = rng.uniform(0.0, 2 * np.pi, size=n)       # 圆周上角度均匀
            if n:
                xy_parts.append(np.column_stack([
                    center[0] + r * np.cos(theta),
                    center[1] + r * np.sin(theta),
                ]))
        xy = np.vstack(xy_parts) if xy_parts else np.empty((0, 2))
        n = xy.shape[0]
        xyz = np.column_stack([xy, np.zeros(n)])              # z = 0
        brightness = np.full(n, self.brightness)
        sizes = sample_size(rng, n, self.size)
        return EmitterSet(xyz=xyz, brightness=brightness, size=sizes)
