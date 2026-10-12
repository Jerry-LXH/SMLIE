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

**鸭子类型 (Option B)**: 每个 Transform 既可吃 ImageStack 也可直接吃 ndarray——
    stage(ImageStack) -> ImageStack   # 进 Pipeline, 自动维护 pixel_size / metadata
    stage(ndarray)    -> ndarray      # notebook 里当"函数"直接调用, 数组进数组出,
                                      # 便于缓存/复用中间结果、避免重复计算
子类只需实现 `_apply(data)->ndarray`(纯数组运算); 需要更新 pixel_size/metadata /
改变帧数结构时重写 `_wrap(stack, data)`; `__call__` 负责鸭子分发。
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import replace
import numpy as np
from ..datatypes import ImageStack


class Transform(ABC):
    """图像处理单元: ImageStack->ImageStack, 也可直接 ndarray->ndarray。

    - 进 ndarray : 走 `_apply`, 返回 ndarray (无 pixel/metadata 语义, 纯计算);
    - 进 ImageStack: `_apply` 后经 `_wrap` 保留/更新 pixel_size / metadata。
    """

    @abstractmethod
    def _apply(self, data: np.ndarray) -> np.ndarray: ...

    def _wrap(self, stack: ImageStack, data: np.ndarray) -> ImageStack:
        """默认: 形状/像素/元数据都不变, 仅替换 data。"""
        return replace(stack, data=data)

    def __call__(self, x):
        if isinstance(x, ImageStack):
            return self._wrap(x, self._apply(np.asarray(x.data)))
        return self._apply(np.asarray(x))

    def __repr__(self):
        return f"{type(self).__name__}()"


class Pipeline(Transform):
    """把多个 Transform 依次作用在 ImageStack 上 (ndarray 亦可直接喂入)。"""
    def __init__(self, *transforms):
        self.transforms = transforms

    def _apply(self, data):
        for t in self.transforms:          # 纯数组链 (__call__ 未接管 ndarray 时兜底)
            data = t(data)
        return data

    def __call__(self, x):
        if isinstance(x, ImageStack):
            for t in self.transforms:      # 逐个吃 ImageStack -> 保留 pixel/metadata 链
                x = t(x)
            return x
        return self._apply(x)              # ndarray -> 纯计算链, 返回 ndarray

    def __repr__(self):
        return f"Pipeline({', '.join(repr(t) for t in self.transforms)})"
