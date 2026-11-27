# smlm/multiplex/__init__.py
"""多通道 (双视角) 数据的配准工具。

分工原则
--------
* **图像**只做漂移校正(全局整数/亚像素平移, 光度中性), 不做形变重采样 ——
  空间相关的形变会让 PSF 插值核随位置变化, 从而引入逐分子的 kappa 抖动,
  这种抖动无法被一个全局 gamma 吸收, 会直接broaden FRET 分布。
* **坐标**才做完整标定(mov -> ref 的 2D 多项式), 它是点级别的, 没有插值损失。

典型流程
--------
标定(离体标准样品), 见 ``smTest_dual_calibration.ipynb``::

    from smlm.multiplex import (MultichannelROI, pair_multichannel_files,
                                fit_multichannel_pairs,
                                MultiChannelCalibrationBundle)
    roi = MultichannelROI(y=(0, 512), x=(10, 246),
                          offsets=((0, 0), (0, 256)), channels=('g', 'r'))
    pairs = pair_multichannel_files(['std2.sif', 'std.sif'], roi, reference='g')
    calibs = fit_multichannel_pairs(pairs, roi)
    MultiChannelCalibrationBundle.from_fit(calibs, roi, files).save('cal.npz')

分析(实验数据), 见 ``smTest_dual.ipynb``::

    from smlm.multiplex import (MultichannelROI, estimate_drift_multi,
                                correct_multi, analyze_multichannel_emitters,
                                extract_multichannel_traces, map_channel_xy)
"""
from __future__ import annotations

from .transform import DualViewCalibration, NPZ_VERSION
from .bundle import MultiChannelCalibrationBundle
from .pair import (PairSet, match_pairs, pair_channels, pair_multichannels,
                   merge_pairs, DEFAULT_TOLS)
from .fit import (poly_terms, n_params, design, fit_transform,
                  cv_select_order, MIN_PAIRS_PER_PARAM)
from .roi import MultichannelROI
from .pipeline import (localize_channel, mean_frame,
                       pair_multichannel_files, fit_multichannel_pairs,
                       multichannel_stack_from_calibration)
from .coords import DualViewGeometry, calib_windows_from_calibs, map_channel_xy
from .drift import (CropWindow, correct_multi, drift_crop_window,
                    drift_requirements, estimate_drift_multi,
                    merge_registered)
from .emitters import (analyze_multichannel_emitters, channel_locs,
                       channel_names, emitter_centers_xy)
from .traces import extract_multichannel_traces
from . import viz
from .viz import (plot_displacement_field, plot_residual_hist,
                  plot_residual_map, plot_overview,
                  plot_multichannel_localizations, plot_multichannel_roi,
                  plot_multichannel_crop, plot_multi_drift, plot_drift_fov,
                  plot_channel_traces)

__all__ = [
    # 变换
    'DualViewCalibration', 'NPZ_VERSION',
    'MultiChannelCalibrationBundle',
    # 配对
    'PairSet', 'match_pairs', 'pair_channels', 'pair_multichannels',
    'merge_pairs', 'DEFAULT_TOLS',
    # 拟合
    'poly_terms', 'n_params', 'design', 'fit_transform',
    'cv_select_order', 'MIN_PAIRS_PER_PARAM',
    # 取框
    'MultichannelROI',
    # 流程
    'localize_channel', 'mean_frame',
    'pair_multichannel_files', 'fit_multichannel_pairs',
    'multichannel_stack_from_calibration',
    # 坐标层
    'DualViewGeometry', 'calib_windows_from_calibs', 'map_channel_xy',
    # 多通道漂移
    'CropWindow', 'merge_registered', 'estimate_drift_multi', 'correct_multi',
    'drift_crop_window', 'drift_requirements',
    # 多通道 emitter 聚类
    'analyze_multichannel_emitters', 'channel_locs', 'channel_names',
    'emitter_centers_xy', 'extract_multichannel_traces',
    # 可视化
    'viz', 'plot_displacement_field', 'plot_residual_hist', 'plot_residual_map',
    'plot_overview', 'plot_multichannel_localizations',
    'plot_multichannel_roi', 'plot_multichannel_crop',
    'plot_multi_drift', 'plot_drift_fov', 'plot_channel_traces',
]
