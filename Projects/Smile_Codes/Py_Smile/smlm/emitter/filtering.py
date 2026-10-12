# =====================================================================
#  过滤器: 基类 + 可组合管线
# =====================================================================
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Iterable
import numpy as np

from ..datatypes import Emitter


class EmitterFilter(ABC):
    @abstractmethod
    def keep(self, emitters) -> np.ndarray: ...          # bool mask

    def __call__(self, emitters):
        m = self.keep(emitters)
        kept = [e for e, k in zip(emitters, m) if k]
        print(f'{type(self).__name__}: 保留 {len(kept)} / {len(emitters)}')
        return kept


class FilterPipeline(EmitterFilter):
    def __init__(self, *filters):
        self.filters = filters
    def keep(self, emitters):                            # 满足接口
        return np.ones(len(emitters), bool)
    def __call__(self, emitters):
        for filt in self.filters:
            emitters = filt(emitters)
        return emitters


class MinDurationFilter(EmitterFilter):
    """短寿命过滤。mode='total' 计总点数, 'consecutive' 计最长连续段。"""
    def __init__(self, min_frames, mode='consecutive'):
        self.min_frames = min_frames; self.mode = mode
    def _longest_run(self, fr):
        if fr.size == 0: return 0
        breaks = np.where(np.diff(fr) > 1)[0]
        seg = np.diff(np.concatenate([[-1], breaks, [fr.size - 1]]))
        return int(seg.max())
    def keep(self, emitters):
        if self.mode == 'consecutive':
            v = np.array([self._longest_run(e.frame) for e in emitters])
        else:
            v = np.array([e.n_locs for e in emitters])
        return v >= self.min_frames


class FirstFrameFilter(EmitterFilter):
    """剔除第一帧不存在的 emitter (只保留拍摄开始就已亮的)。"""
    def keep(self, emitters):
        if not emitters:
            return np.array([], bool)
        first = min(e.on_frame for e in emitters)
        return np.array([e.on_frame == first for e in emitters])


class EccentricityFilter(EmitterFilter):
    """用位置协方差特征值比剔除双-PSF 拟合伪影 (各向异性散点)。
    采用 Wishart 小样本校正, 与 MATLAB 版一致但向量化。"""
    def __init__(self, ecc_th=3.0, min_pts=8):
        self.ecc_th = ecc_th; self.min_pts = min_pts
    def keep(self, emitters):
        m = np.ones(len(emitters), bool)
        for k, e in enumerate(emitters):
            n = e.n_locs
            if n < self.min_pts:
                continue
            ecc_null = 1 + 2 * np.sqrt(2 / n)            # H0 下期望比
            if e.eccentricity() > self.ecc_th * ecc_null:
                m[k] = False
        return m