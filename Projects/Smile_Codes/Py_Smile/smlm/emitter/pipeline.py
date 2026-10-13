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
                     min_frames=5, first_frame_range=None, ecc_th=3.0,
                     do_merge=False, r_merge=None):
    """从 Localizations 连接出 emitter, 再跑一遍过滤管线。

    Parameters
    ----------
    locs : Localizations
        定位结果。
    n_frames : int | None
        帧数; None 表示按 locs 里出现的最大帧自动推断。
    search_radius, max_gap, motion, drift_per_frame
        传给 GapClosingLinker 的连接参数。
    min_frames : int
        MinDurationFilter 阈值(最长连续段长度)。
    first_frame_range : (start, end) | None
        FirstFrameFilter 的帧号窗口(含两端): 只保留首次出现帧 on_frame
        落在窗口内的 emitter; 任一端可为 None 表示不设限。
        None (默认) 保持旧行为, 即只保留在数据最早帧就亮的 emitter。
    ecc_th : float
        EccentricityFilter 阈值。
    do_merge : bool
        是否先用 EmitterMerger 合并。
    r_merge : float | None
        合并半径; None -> 1.5 * search_radius。

    Returns
    -------
    (kept, all_emitters) : tuple
        过滤后保留的 emitter 列表, 以及过滤前的全部 emitter。
    """
    linker = GapClosingLinker(search_radius, max_gap, motion, drift_per_frame)
    emitters = linker.link(locs, n_frames)
    print(f'连接得到 {len(emitters)} 个 emitter')

    if do_merge:
        emitters = EmitterMerger(r_merge or 1.5 * search_radius)(emitters)

    pipeline = FilterPipeline(
        MinDurationFilter(min_frames, 'total'),
        FirstFrameFilter(first_frame_range),
        EccentricityFilter(ecc_th),
    )
    return pipeline(emitters), emitters