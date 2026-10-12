# smlm/imgprocess/denoise.py
"""去噪类 Transform。

RIED 的第一步就是高斯预滤波 (sigma ~ 0.7 px) 抑制高频采样噪声。
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.ndimage import gaussian_filter
from .base import Transform


class GaussianFilter(Transform):
    """逐帧高斯平滑 (对应 RIED 的高斯预滤波)。

    sigma: 高斯核标准差 (pixel)。Matlab `imfilter(data, gaussian(sigma), 'replicate')`
           等价于 scipy `gaussian_filter(mode='nearest')`。
    """
    def __init__(self, sigma=0.7):
        self.sigma = float(sigma)

    def __call__(self, stack):
        if self.sigma == 0:
            return stack
        data = np.asarray(stack.data, dtype=np.float32)
        out = np.empty_like(data)
        for i in range(data.shape[0]):
            out[i] = gaussian_filter(data[i], self.sigma, mode='nearest')
        return replace(stack, data=out)

    def __repr__(self):
        return f"GaussianFilter(sigma={self.sigma})"
