from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from ..datatypes import EmitterSet, FieldOfView


class PSFModel(ABC):
    @abstractmethod
    def render(self, emitters: EmitterSet, photon_frame: np.ndarray,
               fov: FieldOfView) -> np.ndarray:
        """把某一帧的光子铺到像面，返回理想光子图 (C, H, W)。
        photon_frame: (N, C) 该帧每个 emitter 每通道的光子数。"""

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
