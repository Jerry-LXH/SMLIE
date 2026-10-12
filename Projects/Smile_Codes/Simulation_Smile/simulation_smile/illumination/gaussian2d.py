from __future__ import annotations
import numpy as np
from .base import Illumination


class Gaussian2D(Illumination):
    """横向二维高斯激发，z 方向不衰减（Phase 2 再加 TIRF 等）。"""
    def __init__(self, center_xy, sigma: float, peak: float = 1.0):
        self.center = np.asarray(center_xy, dtype=float)  # (2,) nm
        self.sigma = sigma
        self.peak = peak

    def intensity(self, xyz, t=0.0):
        d2 = np.sum((xyz[:, :2] - self.center) ** 2, axis=1)
        return self.peak * np.exp(-d2 / (2 * self.sigma ** 2))
