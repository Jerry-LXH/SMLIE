from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from ..datatypes import EmitterSet, PhotonStream


class PhotophysicsModel(ABC):
    @abstractmethod
    def simulate(self, emitters: EmitterSet, illumination,
                 n_frames: int, dt: float,
                 rng: np.random.Generator) -> PhotonStream: ...

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
