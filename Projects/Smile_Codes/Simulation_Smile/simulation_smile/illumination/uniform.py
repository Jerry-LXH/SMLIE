from __future__ import annotations
import numpy as np
from .base import Illumination


class Uniform(Illumination):
    def __init__(self, level: float = 1.0):
        self.level = level

    def intensity(self, xyz, t=0.0):
        return np.full(xyz.shape[0], self.level, dtype=float)
