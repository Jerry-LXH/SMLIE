from __future__ import annotations
import numpy as np
from ..datatypes import Emitter

class _TrackBuilder:
    """连接过程中的可变轨迹, 内含速度估计以支持运动预测。"""
    __slots__ = ('frames', 'xs', 'ys', 'idxs',
                 'last_frame', 'last_xy', 'vel', 'confirmed_bleached')

    def __init__(self, f, x, y, idx):
        self.frames = [f]; self.xs = [x]; self.ys = [y]; self.idxs = [idx]
        self.last_frame = f
        self.last_xy = np.array([x, y], float)
        self.vel = np.zeros(2)
        self.confirmed_bleached = False

    def append(self, f, x, y, idx):
        new = np.array([x, y], float)
        dt = f - self.last_frame
        if dt > 0:                                    # 指数平滑速度
            self.vel = 0.5 * self.vel + 0.5 * (new - self.last_xy) / dt
        self.frames.append(f); self.xs.append(x)
        self.ys.append(y); self.idxs.append(idx)
        self.last_frame = f; self.last_xy = new

    def predict(self, f, motion) -> np.ndarray:
        if motion == 'linear':
            return self.last_xy + self.vel * (f - self.last_frame)
        return self.last_xy                            # 'static'

    def to_emitter(self) -> Emitter:
        fr = np.asarray(self.frames)
        o = np.argsort(fr)
        return Emitter(fr[o], np.asarray(self.xs)[o], np.asarray(self.ys)[o],
                       np.asarray(self.idxs)[o], self.confirmed_bleached)