from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
from ..datatypes import EmitterSet, FieldOfView


class PSFModel(ABC):
    @abstractmethod
    def render(self, emitters: EmitterSet, photon_frame: np.ndarray,
               fov: FieldOfView,
               rng: np.random.Generator | None = None) -> np.ndarray:
        """把某一帧的光子按 PSF 多项分布随机分配到像面, 返回 (C, H, W) 计数图。
        photon_frame: (N, C) 该帧每个 emitter 每通道的光子数(应为整数)。
        对每个 emitter 的 N 个光子, 以归一化 PSF 核为类别概率做多项分配
        (逐光子决定落点像素), 返回各像素计数直方图——shot noise 已包含。
        rng: 外部传入 Generator 以保证可复现; None 时内部新建 default_rng()。"""

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
