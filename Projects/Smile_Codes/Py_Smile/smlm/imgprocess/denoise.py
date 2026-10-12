# smlm/imgprocess/denoise.py
"""去噪类 Transform。

RIED 的第一步就是高斯预滤波 (sigma ~ 0.7 px) 抑制高频采样噪声。
"""
from __future__ import annotations
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

    def _apply(self, data):
        if self.sigma == 0:
            return np.asarray(data)              # sigma=0: 原样返回 (不cast)
        data = np.asarray(data, dtype=np.float32)
        if data.ndim == 2:
            return gaussian_filter(data, self.sigma, mode='nearest')
        out = np.empty_like(data)
        for i in range(data.shape[0]):
            out[i] = gaussian_filter(data[i], self.sigma, mode='nearest')
        return out

    def __repr__(self):
        return f"GaussianFilter(sigma={self.sigma})"
