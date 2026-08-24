from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from ..datatypes import EmitterSet


def sample_size(rng, n: int, size):
    """按规格采样粒子物理尺寸 (nm)。

    size: (mean_nm, variance_nm²) 或 None。
      - 提供时: 每个粒子 size ~ N(mean, var), 物理上尺寸非负, 截断到 >= 0;
      - None 时: 返回 None (与 EmitterSet.size 默认一致)。
    返回 (n,) 数组或 None。
    """
    if size is None:
        return None
    mean, var = float(size[0]), float(size[1])
    s = rng.normal(loc=mean, scale=np.sqrt(max(var, 0.0)), size=n)
    return np.maximum(s, 0.0)


class EmitterSource(ABC):
    @abstractmethod
    def sample(self, rng: np.random.Generator) -> EmitterSet: ...

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
