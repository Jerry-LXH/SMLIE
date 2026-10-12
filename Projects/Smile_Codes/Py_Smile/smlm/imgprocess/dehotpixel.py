# smlm/imgprocess/dehotpixel.py
"""sCMOS 热点像素去除 (对应 RIEDm.m 的第一步)。

sCMOS 传感器上少数像素有异常高的暗电流/增益 -> 逐帧出现"热点"。
RIED 原版做法: 整卷上 `data(data > thresh) = data(data > thresh) / divisor`。
"""
from __future__ import annotations
import numpy as np
from .base import Transform


def dehotpixel(data, thresh=5000.0, divisor=50.0):
    """把超过 thresh 的像素值除以 divisor (返回新数组, 不原地改)。

    data : (T, H, W) 或 (H, W) float 数组 (RIED 在整卷上一次性处理)。
    """
    data = np.asarray(data, dtype=np.float32)
    mask = data > thresh
    if not mask.any():
        return data
    out = data.copy()
    out[mask] = data[mask] / divisor
    return out


class Dehotpixel(Transform):
    """热点去除 Transform (ImageStack -> ImageStack)。

    thresh : 判定为热点的阈值 (默认 5000, 对齐 RIED 原版; 单位随数据, 对
             真实 ADU 数据用 ADU, 对光子数据用光子)。
    divisor: 热点值的衰减倍数 (默认 50, 对齐 RIED 原版)。
    """
    def __init__(self, thresh=5000.0, divisor=50.0):
        self.thresh = float(thresh)
        self.divisor = float(divisor)

    def _apply(self, data):
        return dehotpixel(data, self.thresh, self.divisor)

    def __repr__(self):
        return f"Dehotpixel(thresh={self.thresh:g}, divisor={self.divisor:g})"
