# smlm/imgprocess/base.py
"""imgprocess 的 Transform 抽象与流水线组合。

约定: 凡 ImageStack 进、ImageStack 出的处理单元都可以放进 `smlm.imgprocess`:
  - 纯预处理 (裁剪/选帧)         -> preprocess.py
  - 插值/上采样/下采样           -> interp.py
  - 去噪 (高斯等)                -> denoise.py
  - 反卷积 (RL / 稀疏 / Wiener)  -> deconv.py
  - SOFI / RIED 类超分辨核心     -> sofi.py  (Sofit 也是 Transform 子类)
  - 单纯可视化                   -> viz.py

Transform 是最小的可组合单元; Pipeline 把多个 Transform 串起来。
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from ..datatypes import ImageStack


class Transform(ABC):
    """ImageStack -> ImageStack 的抽象基类。"""
    @abstractmethod
    def __call__(self, stack: ImageStack) -> ImageStack: ...


class Pipeline(Transform):
    """把多个 Transform 依次作用在 ImageStack 上。"""
    def __init__(self, *transforms):
        self.transforms = transforms

    def __call__(self, stack):
        for t in self.transforms:
            stack = t(stack)
        return stack

    def __repr__(self):
        return f"Pipeline({', '.join(repr(t) for t in self.transforms)})"
