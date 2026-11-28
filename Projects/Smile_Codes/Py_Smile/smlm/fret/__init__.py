# smlm/fret/__init__.py
"""双通道 smFRET: 取 trace -> 扫事件 -> 标定 γ/α/β -> 算 E -> 人工复核 -> 落盘。

风格与包内其他子包一致: **数据和算法用 dataclass / 纯函数**。

典型流程
--------
>>> import smlm.fret as fr
>>> import smlm.trace as tr
>>> fillers = [tr.ApertureRaw(mode='hard', radius=3),        # 硬圆盘不依赖 σ
...            tr.AnnulusRing(gap=2, width=4, excl_radius=3),
...            tr.IntensityFromRawBg()]
>>> tr_r = tr.TraceFillPipeline(corr_r.data, c_r, sigma=1.6, tag='r').run(fillers)
>>> tr_g = tr.TraceFillPipeline(corr_g.data, c_g, centers_for_bg=mask_g,
...                             sigma=1.6, tag='g').run(fillers)
>>> events = fr.scan_steps(tr_r, tr_g, fr.EventSpec())
>>> cal = fr.calibrate(events)
>>> res, stats = fr.fret_efficiency(events, tr_r, tr_g, cal, t_bleach=lambda j: j*0.1)
>>> fr.save_csv(res, 'x_fret.csv', cal)

约定
----
* ``trace_a`` 永远是**受体(红/acceptor)**, ``trace_d`` 是**供体(绿/donor)**。
* 孔径中心统一由**红通道定位**映射到各通道, 见 ``smlm.multiplex.traces``。
"""
from __future__ import annotations

# --- 事件 ---
from .events import EventSpec, find_step, scan_steps, red_state, RedState

# --- trace 对的取舍 ---
from .pairs import (PairCriteria, PairRecord, select_pairs, attach_efficiency,
                    noise_per_frame, pre_spread, donor_stats)

# --- 标定 ---
from .calibrate import (FretCalibration, calibrate, estimate_alpha,
                        estimate_gamma, estimate_beta, robust_line,
                        local_step, gamma_points)

# --- 效率 ---
from .efficiency import EfficiencySpec, fret_efficiency, window_mean

# --- 导出 ---
from .datatypes import StepEvent, FretResult, CRITERIA, CRITERIA_STRICT, passed
from .io import TracesByID, save_csv, save_traces_csv, load_traces_csv, TRACE_COLS
from .export import save_bundle, FILES, filenames

from . import viz

__all__ = [
    # 事件
    'EventSpec', 'find_step', 'scan_steps', 'StepEvent', 'passed',
    'red_state', 'RedState',
    'CRITERIA', 'CRITERIA_STRICT',
    # trace 对的取舍
    'PairCriteria', 'PairRecord', 'select_pairs', 'attach_efficiency',
    'noise_per_frame', 'pre_spread', 'donor_stats',
    # 标定
    'FretCalibration', 'calibrate', 'estimate_alpha', 'estimate_gamma',
    'estimate_beta', 'robust_line', 'local_step', 'gamma_points',
    # 效率
    'EfficiencySpec', 'fret_efficiency', 'window_mean', 'FretResult',
    # 导出
    'save_csv', 'viz', 'save_traces_csv', 'load_traces_csv', 'TracesByID',
    'TRACE_COLS',
    'save_bundle', 'FILES', 'filenames',
]
