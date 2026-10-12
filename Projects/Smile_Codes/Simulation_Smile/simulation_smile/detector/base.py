from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np


class Detector(ABC):
    @abstractmethod
    def capture(self, photon_image: np.ndarray,
                rng: np.random.Generator) -> np.ndarray:
        """输入理想光子图 (C, H, W)，返回数字图像 (C, H, W)。"""

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
