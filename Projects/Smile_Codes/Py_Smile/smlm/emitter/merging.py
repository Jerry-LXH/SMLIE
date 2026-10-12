# =====================================================================
#  合并器 (union-find, 比 MATLAB 的原地改 D 更清晰)
# =====================================================================
from __future__ import annotations
import numpy as np
from scipy.spatial.distance import cdist
from ..datatypes import Emitter

class EmitterMerger:
    """合并空间重叠但时间上几乎不重叠(闪烁断段)的 emitter。"""
    def __init__(self, r_merge=1.5, max_overlap=0.3):
        self.r_merge = r_merge
        self.max_overlap = max_overlap

    def __call__(self, emitters):
        N = len(emitters)
        if N < 2:
            return emitters
        cent = np.array([e.centroid for e in emitters])
        D = cdist(cent, cent)
        np.fill_diagonal(D, np.inf)

        parent = list(range(N))
        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]; a = parent[a]
            return a

        for i, j in zip(*np.where(np.triu(D < self.r_merge))):
            fi, fj = set(emitters[i].frame.tolist()), set(emitters[j].frame.tolist())
            ov = len(fi & fj) / len(fi | fj)
            if ov <= self.max_overlap:                    # 重叠少 → 同一分子闪烁
                parent[find(i)] = find(j)

        groups: dict[int, list[int]] = {}
        for k in range(N):
            groups.setdefault(find(k), []).append(k)

        merged = []
        for members in groups.values():
            if len(members) == 1:
                merged.append(emitters[members[0]]); continue
            fr = np.concatenate([emitters[m].frame for m in members])
            xs = np.concatenate([emitters[m].x for m in members])
            ys = np.concatenate([emitters[m].y for m in members])
            ix = np.concatenate([emitters[m].loc_idx for m in members])
            o = np.argsort(fr)
            merged.append(Emitter(fr[o], xs[o], ys[o], ix[o],
                                  all(emitters[m].bleached for m in members)))
        print(f'合并后剩余 {len(merged)} / {N} 个 emitter')
        return merged


