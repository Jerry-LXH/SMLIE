# smlm/fret/__init__.py
"""双通道 smFRET: 取 trace -> 扫事件 -> 标定 γ/α/β -> 算 E -> 人工复核 -> 落盘。

风格与包内其他子包一致: **数据和算法用 dataclass / 纯函数, 只有有状态的东西
(取 trace 的缓存、复核会话)才用类**。交互层单独隔离在 ``widgets``。

典型流程
--------
>>> import smlm.fret as fr
>>> spec = fr.ApertureSpec(mode='matched', k_sigma=2.0)      # 匹配滤波
>>> tr_r = fr.extract_channel(corr_r.data, c_r, spec=spec, locs=locs_r_c, tag='r')
>>> tr_g = fr.extract_channel(corr_g.data, c_g, spec=spec, locs=locs_g_c,
...                           mask_centers_yx=mask_g, tag='g')
>>> events = fr.scan_steps(tr_r, tr_g, fr.EventSpec())
>>> cal = fr.calibrate(events)
>>> res, stats = fr.fret_efficiency(events, tr_r, tr_g, cal, t_bleach=lambda j: j*0.1)
>>> fr.save_csv(res, 'x_fret.csv', cal)

约定
----
* ``trace_a`` 永远是**受体(红/acceptor)**, ``trace_d`` 是**供体(绿/donor)**。
* 孔径中心统一由**红通道定位**映射到各通道, 见 ``extract``。
"""
from __future__ import annotations

# --- 孔径 / 背景 ---
from .aperture import ApertureSpec, Aperture, crop_centered, sigma_per_molecule
from .background import (AnnulusBackground, ring_mask, clip_median,
                         valid_mask, ring_budget, moving_average,
                         model_weights, emitter_frame_photons, emitter_sigma)

# --- 取 trace ---
from .extract import extract_channel, PairedExtractor, ChannelTraces

# --- 事件 ---
from .events import EventSpec, find_step, scan_steps

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
from .io import save_csv, load_csv

from . import viz

__all__ = [
    # 孔径 / 背景
    'ApertureSpec', 'Aperture', 'crop_centered', 'sigma_per_molecule',
    'AnnulusBackground', 'ring_mask', 'clip_median', 'valid_mask', 'ring_budget',
    'moving_average', 'model_weights', 'emitter_frame_photons', 'emitter_sigma',
    # 取 trace
    'extract_channel', 'PairedExtractor', 'ChannelTraces',
    # 事件
    'EventSpec', 'find_step', 'scan_steps', 'StepEvent', 'passed',
    'CRITERIA', 'CRITERIA_STRICT',
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
]
