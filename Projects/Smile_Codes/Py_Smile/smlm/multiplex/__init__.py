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
标定(离体标准样品)::

    from smlm.multiplex import calibrate_dual_view, DualViewROI
    calib, rep = calibrate_dual_view(['std2.sif', 'std.sif'],
                                     roi=DualViewROI(split_x=256,
                                                     row_range=(0, 512)))
    calib.save('dualcal_v2.npz')

分析(实验数据)::

    from smlm.multiplex import load_multi, DualViewCalibration
    ms = load_multi('sample.sif', margin=12)     # 两通道同一物理窗口
    calib = DualViewCalibration.load('dualcal_v2.npz')
    xy_ref = calib.transform(xy_red)             # red local -> green local
    locs_ref = calib.apply_locs(locs_r)          # 对 Localizations 一步到位
"""
from __future__ import annotations

from .transform import DualViewCalibration, NPZ_VERSION
from .bundle import MultiChannelCalibrationBundle
from .pair import (PairSet, match_pairs, pair_channels, pair_multichannels,
                   merge_pairs, DEFAULT_TOLS)
from .fit import (poly_terms, n_params, design, fit_transform, residual_report,
                  cv_select_order, MIN_PAIRS_PER_PARAM)
from .channels import (DualViewROI, split_dual_view, cut_pair, build_multi,
                       cut_channels, multi_from_crops,
                       refine_coarse_offset, check_drift_margin, seam_windows,
                       split_from_crops, windows_from_crops,
                       mirror_crop, crops_from_ref, analysis_crops,
                       DEFAULT_SPLIT_X, DEFAULT_ROW_RANGE, DEFAULT_SEAM,
                       SEAM_HINT)
from .roi import MultichannelROI
from .pipeline import (calibrate_dual_view, calib_one_file, load_multi,
                       localize_channel, mean_frame,
                       pair_multichannel_files, fit_multichannel_pairs,
                       multichannel_stack_from_calibration)
from .coords import (DualViewGeometry, geometry_from_multi, geometry_from_crops,
                     calib_windows_from_roi, window_from_range, map_channel_xy)
from .drift import (CropWindow, correct_multi, drift_crop_window,
                    drift_requirements, estimate_drift_multi,
                    merge_registered)
from .emitters import (analyze_multichannel_emitters, channel_locs,
                       channel_names, emitter_centers_xy)
from .traces import extract_multichannel_traces
from . import viz
from .viz import (plot_displacement_field, plot_residual_hist,
                  plot_residual_map, plot_overview, plot_channel_check,
                  plot_cut_windows, plot_multichannel_localizations,
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
    'poly_terms', 'n_params', 'design', 'fit_transform', 'residual_report',
    'cv_select_order', 'MIN_PAIRS_PER_PARAM',
    # 切通道 / 余量
    'DualViewROI', 'split_dual_view', 'cut_pair', 'build_multi',
    'MultichannelROI',
    'cut_channels', 'multi_from_crops',
    'refine_coarse_offset', 'check_drift_margin', 'seam_windows',
    'split_from_crops', 'windows_from_crops',
    'mirror_crop', 'crops_from_ref', 'analysis_crops',
    'DEFAULT_SPLIT_X', 'DEFAULT_ROW_RANGE', 'DEFAULT_SEAM', 'SEAM_HINT',
    # 流程
    'calibrate_dual_view', 'calib_one_file', 'load_multi', 'localize_channel',
    'mean_frame', 'pair_multichannel_files', 'fit_multichannel_pairs',
    'multichannel_stack_from_calibration',
    # 坐标层
    'DualViewGeometry', 'geometry_from_multi', 'geometry_from_crops',
    'calib_windows_from_roi', 'window_from_range', 'map_channel_xy',
    # 多通道漂移
    'CropWindow', 'merge_registered', 'estimate_drift_multi', 'correct_multi',
    'drift_crop_window', 'drift_requirements',
    # 多通道 emitter 聚类
    'analyze_multichannel_emitters', 'channel_locs', 'channel_names',
    'emitter_centers_xy', 'extract_multichannel_traces',
    # 可视化
    'viz', 'plot_displacement_field', 'plot_residual_hist', 'plot_residual_map',
    'plot_overview', 'plot_channel_check', 'plot_cut_windows',
    'plot_multichannel_localizations', 'plot_multichannel_crop',
    'plot_multi_drift', 'plot_drift_fov', 'plot_channel_traces',
]
