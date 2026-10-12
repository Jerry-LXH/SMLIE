from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np


class Illumination(ABC):
    @abstractmethod
    def intensity(self, xyz: np.ndarray, t: float = 0.0) -> np.ndarray:
        """给定位置 (N,3)，返回各点的相对激发强度 (N,)。"""

    @property
    def config(self) -> dict:
        """自动收集所有非下划线开头的实例属性,做 JSON 友好转换。"""
        out = {}
        for k, v in vars(self).items():
            if k.startswith("_"):
                continue
            if isinstance(v, np.ndarray):
                v = v.tolist()
            out[k] = v
        return out
