# smlm/imgprocess/preprocess.py
"""纯预处理 Transform: 空间裁剪 与 时间选帧。

原 `smlm/preprocess/{spatial,temporal}.py` 合并到这一个文件。
"""
from __future__ import annotations
from dataclasses import replace
from .base import Transform


class Crop(Transform):
    """按 (start, stop) 裁剪空间范围;None 表示该维度不裁。"""
    def __init__(self, y=None, x=None):
        self.y, self.x = y, x

    def __call__(self, stack):
        ys = slice(*self.y) if self.y else slice(None)
        xs = slice(*self.x) if self.x else slice(None)
        meta = dict(stack.metadata)
        # ⚠️ 记录裁剪偏移,定位坐标才能还原到相机原始坐标系
        meta['crop_offset'] = (self.y[0] if self.y else 0,
                               self.x[0] if self.x else 0)
        return replace(stack, data=stack.data[:, ys, xs], metadata=meta)


class SelectFrames(Transform):
    """选帧:start:stop:step。"""
    def __init__(self, start=0, stop=None, step=1):
        self.start, self.stop, self.step = start, stop, step

    def __call__(self, stack):
        data = stack.data[self.start:self.stop:self.step, :, :]
        return replace(stack, data=data)
