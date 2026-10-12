from __future__ import annotations
from typing import Literal
import numpy as np

from ..datatypes import Localizations, Emitter
from .linking import GapClosingLinker
from .merging import EmitterMerger
from .filtering import (FilterPipeline, MinDurationFilter,
                        FirstFrameFilter, EccentricityFilter)

def analyze_emitters(locs, n_frames, *,
                     search_radius=1.0, max_gap=30,
                     motion='linear', drift_per_frame=0.0,
                     min_frames=5, ecc_th=3.0,
                     do_merge=False, r_merge=None):
    linker = GapClosingLinker(search_radius, max_gap, motion, drift_per_frame)
    emitters = linker.link(locs, n_frames)
    print(f'连接得到 {len(emitters)} 个 emitter')

    if do_merge:
        emitters = EmitterMerger(r_merge or 1.5 * search_radius)(emitters)

    pipeline = FilterPipeline(
        MinDurationFilter(min_frames, 'consecutive'),
        FirstFrameFilter(),
        EccentricityFilter(ecc_th),
    )
    return pipeline(emitters), emitters