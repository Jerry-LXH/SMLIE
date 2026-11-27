# smlm/fret/__init__.py
"""双通道 smFRET: 取 trace -> 扫事件 -> 标定 γ/α/β -> 算 E -> 人工复核 -> 落盘。

风格与包内其他子包一致: **数据和算法用 dataclass / 纯函数, 只有有状态的东西
(取 trace 的缓存、复核会话)才用类**。交互层单独隔离在 ``widgets``。

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

# --- 孔径 / 背景 ---
from .aperture import ApertureSpec, Aperture, crop_centered, sigma_per_molecule
from .background import (AnnulusBackground, ring_mask, clip_median,
                         valid_mask, ring_budget, moving_average,
                         model_weights, emitter_frame_photons, emitter_sigma)

# --- 事件 ---
from .events import EventSpec, find_step, scan_steps

# --- trace 对的取舍 ---
from .pairs import (PairCriteria, PairRecord, select_pairs, attach_efficiency,
                    noise_per_frame, pre_spread, donor_stats)

# --- 标定 ---
from .calibrate import (FretCalibration, calibrate, estimate_alpha,
                        estimate_gamma, estimate_beta, robust_line,
                        select_for_gamma, local_step, gamma_points)

# --- 效率 ---
from .efficiency import EfficiencySpec, fret_efficiency, window_mean

# --- 诊断 ---
from .qc import mle_track, mle_spatial, flux_ratio, flag

# --- 复核 / 导出 ---
from .review import ReviewSession
from .datatypes import StepEvent, FretResult, CRITERIA, CRITERIA_STRICT, passed
from .io import (save_csv, save_traces_csv, load_csv, TRACE_COLS)

from . import viz

__all__ = [
    # 孔径 / 背景
    'ApertureSpec', 'Aperture', 'crop_centered', 'sigma_per_molecule',
    'AnnulusBackground', 'ring_mask', 'clip_median', 'valid_mask', 'ring_budget',
    'moving_average', 'model_weights', 'emitter_frame_photons', 'emitter_sigma',
    # 事件
    'EventSpec', 'find_step', 'scan_steps', 'StepEvent', 'passed',
    'CRITERIA', 'CRITERIA_STRICT',
    # trace 对的取舍
    'PairCriteria', 'PairRecord', 'select_pairs', 'attach_efficiency',
    'noise_per_frame', 'pre_spread', 'donor_stats',
    # 标定
    'FretCalibration', 'calibrate', 'estimate_alpha', 'estimate_gamma',
    'estimate_beta', 'robust_line', 'select_for_gamma', 'local_step',
    'gamma_points',
    # 效率
    'EfficiencySpec', 'fret_efficiency', 'window_mean', 'FretResult',
    # 诊断
    'mle_track', 'mle_spatial', 'flux_ratio', 'flag',
    # 复核 / 导出
    'ReviewSession', 'save_csv', 'load_csv', 'viz',
    'save_traces_csv', 'TRACE_COLS',
]
