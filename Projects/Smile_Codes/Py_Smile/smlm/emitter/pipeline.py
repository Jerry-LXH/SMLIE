from __future__ import annotations
from typing import Literal
import numpy as np

from ..datatypes import Localizations, Emitter
from .linking import GapClosingLinker
from .merging import EmitterMerger
from .filtering import (FilterPipeline, MinDurationFilter,
                        FirstFrameFilter, EccentricityFilter)


def analyze_emitters_sm(locs, n_frames, *,
                        search_radius=1.0, max_gap=30,
                        motion='static', drift_per_frame=0.0, gap_penalty=0.0,
                        min_frames=5, first_frame_range=None,
                        ecc_p=1e-3, ecc_min_pts=8, ecc_min=1.5,
                        do_merge=False, r_merge=None, max_overlap=0.3):
    """闪烁的**固定单分子**分析: 连接 -> 合并 -> 过滤。

    这条链的每一步都是为"一个不动的分子反复亮/灭、最后漂白"设计的:

    * **连接**: 位置不随时间漂移, 位置不确定性只来自定位噪声, 所以用
      ``motion='static'``。``motion='linear'`` 的速度是从逐帧差分平滑出来的,
      对静止分子完全是噪声 (量级 σ/dt), 而预测点还要外推 ``f - last_frame`` 帧 ——
      暗态越长错得越离谱 (100 帧暗态就能外推出十几像素), 直接断链。
      固定分子不要用 linear。
    * **合并**: 断链之后把"位置相同、时间几乎不重叠"的两段接回来。这是补救,
      不是常态 —— 先尽量让 ``max_gap`` 覆盖典型暗态。
    * **过滤**: 偏心率剔除"两个分子被拟合成一个质心"的假轨迹 (质心随相对亮度
      在两点间游走, 点云沿连线拉长)。判据是**精确零分布的 p 值 + 效应量下限**
      (``ecc_p`` 与 ``ecc_min``): 只看 p 值的话, 轨迹越长越容易被误杀 ——
      n=828 时轴比 1.30 就"显著", 而 0.3 px 的残余漂移就能拉到 1.35。
    * **不按距离丢弃近邻**: 真能分辨的近邻分子是真实存在的 (二聚体、FRET 双标记),
      按距离硬删只会误删它们。

    Parameters
    ----------
    locs : Localizations
        定位结果 (已做漂移校正)。
    n_frames : int | None
        帧数; None 表示按 locs 里出现的最大帧自动推断。
    search_radius : float
        逐帧连接的搜索半径 (px)。取定位抖动的数倍, 而不是 CRLB。
    max_gap : int
        能桥接的最长暗态 (帧)。取长于典型暗态; 超出它的断段交给 ``do_merge``。
    motion : 'static' | 'linear'
        固定分子用 ``'static'`` (默认)。
    drift_per_frame : float
        暗态期间搜索半径的扩张速度。漂移已在图像层面校正时应为 0。
    gap_penalty : float
        代价里按暗帧数加罚, 让"刚熄灭"的轨迹优先被连接, 减少错连。
    min_frames : int
        最少定位点数 (``MinDurationFilter`` 的总点数口径)。
    first_frame_range : (start, end) | None
        只保留首次出现帧落在窗口内的 emitter; None = 只保留最早出现的那些。
    ecc_p : float
        偏心率的 p 值阈值 (精确零分布, 默认 1e-3 = 0.1% 假阳性)。
    ecc_min_pts : int
        少于这么多定位点不做偏心率判定 (短轨迹检验功效低, 只是判不出来)。
    ecc_min : float | None
        偏心率的**效应量下限**: 轴比不超过它就不剔除 (默认 1.5)。p 值在长轨迹上
        会因样本量大而在"物理上无关紧要"的拉长处变显著 (实测 n=828 时轴比 1.30
        即 p<1e-3, 而 0.3 px 残余漂移就能达到), 所以必须再要求一定的实际拉长。
        ``None`` = 只看 p 值 (旧行为)。
    do_merge : bool
        是否用 ``EmitterMerger`` 接回断段。
    r_merge : float | None
        合并半径; None -> 1.5 * search_radius。按"两段质心之差"取; 静止分子一两像素
        就够, 不要为了"合并近邻"而放大。
    max_overlap : float
        合并时允许的帧集合 Jaccard 重叠上限; 超过它说明两条轨迹在同时亮,
        不是同一分子的断段。

    Returns
    -------
    (kept, all_emitters) : tuple
        过滤后保留的 emitter 列表, 以及过滤前的全部 emitter。
    """
    linker = GapClosingLinker(search_radius, max_gap, motion, drift_per_frame,
                              gap_penalty)
    emitters = linker.link(locs, n_frames)
    print(f'连接得到 {len(emitters)} 个 emitter')

    if do_merge:
        emitters = EmitterMerger(r_merge or 1.5 * search_radius,
                                 max_overlap)(emitters)

    pipeline = FilterPipeline(
        MinDurationFilter(min_frames, 'total'),
        FirstFrameFilter(first_frame_range),
        EccentricityFilter(ecc_p, min_pts=ecc_min_pts, ecc_min=ecc_min),
    )
    return pipeline(emitters), emitters
