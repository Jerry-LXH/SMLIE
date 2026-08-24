# =====================================================================
#  Linker 基类 + 带运动预测的门控连接器
# =====================================================================

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional
import numpy as np
from scipy.optimize import linear_sum_assignment   # 匈牙利匹配, gap closing 用
from scipy.spatial.distance import cdist

from ..datatypes import Emitter, Localizations
from .base import _TrackBuilder

_BIG = 1e9
class Linker(ABC):
    @abstractmethod
    def link(self, locs, n_frames: Optional[int] = None) -> list[Emitter]:
        ...


class GapClosingLinker(Linker):
    """逐帧全局最优指派 (Hungarian) + 间隙闭合。

    相比 MATLAB 版本的三点改进:
      1. motion='linear' 用速度预测下一帧位置 → 兼容慢速分子运动;
      2. drift_per_frame 在暗态期间按间隙线性放大搜索半径
         (暗越久, 位置越不确定) → 闪烁重亮不易断链;
      3. gap_penalty 让代价随暗帧数增加, 优先连"刚熄灭"的分子。

    这样大部分闪烁/漂移导致的断段在连接阶段就闭合, mergeEmitters 变成可选。
    """
    def __init__(self, search_radius=1.0, max_gap=30,
                 motion='linear', drift_per_frame=0.0, gap_penalty=0.0):
        self.r = float(search_radius)
        self.max_gap = int(max_gap)
        self.motion = motion
        self.drift_per_frame = float(drift_per_frame)
        self.gap_penalty = float(gap_penalty)

    def link(self, locs, n_frames=None):
        frame = np.asarray(locs.frame).round().astype(int)
        x = np.asarray(locs.x, float)
        y = np.asarray(locs.y, float)
        if frame.size == 0:
            return []

        # 一次排序分桶 (替代 MATLAB 每帧的 == f 扫描, O(N log N))
        order = np.argsort(frame, kind='stable')
        fsorted = frame[order]
        uniq = np.unique(fsorted)
        groups = np.split(order, np.searchsorted(fsorted, uniq[1:]))
        buckets = {int(f): g for f, g in zip(uniq, groups)}
        final_frame = int(uniq[-1]) if n_frames is None else int(n_frames)

        active: list[_TrackBuilder] = []
        finished: list[_TrackBuilder] = []

        for f in uniq.tolist():
            f = int(f)
            # --- 退休超时轨迹 (确认漂白) ---
            still = []
            for tb in active:
                if f - tb.last_frame > self.max_gap:
                    tb.confirmed_bleached = True
                    finished.append(tb)
                else:
                    still.append(tb)
            active = still

            idxs = buckets[f]
            pos = np.column_stack([x[idxs], y[idxs]])
            matched = np.zeros(len(idxs), bool)

            if active:
                pred = np.array([tb.predict(f, self.motion) for tb in active])
                gaps = np.array([f - tb.last_frame for tb in active], float)
                d = cdist(pos, pred)                       # [n_loc x n_active]
                radius = self.r + self.drift_per_frame * np.maximum(gaps - 1, 0)
                gate = d < radius[None, :]
                cost = np.where(gate,
                                d + self.gap_penalty * (gaps[None, :] - 1),
                                _BIG)
                ri, ci = linear_sum_assignment(cost)
                for i, j in zip(ri, ci):
                    if gate[i, j]:
                        li = int(idxs[i])
                        active[j].append(f, x[li], y[li], li)
                        matched[i] = True

            # --- 未匹配的建立新轨迹 ---
            for i in np.nonzero(~matched)[0]:
                li = int(idxs[i])
                active.append(_TrackBuilder(f, x[li], y[li], li))

        # 片尾仍活跃: 是否真漂白取决于剩余观察窗是否足够
        for tb in active:
            tb.confirmed_bleached = (final_frame - tb.last_frame) > self.max_gap
            finished.append(tb)

        return [tb.to_emitter() for tb in finished]