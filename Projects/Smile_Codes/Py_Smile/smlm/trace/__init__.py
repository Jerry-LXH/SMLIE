# smlm/trace/__init__.py
"""Trace extraction (aperture integration) and state analysis.

Both the plain single-channel flow and the multichannel / FRET flow use
:func:`extract_channel` from here; ``mask_centers_yx`` lets the caller separate
the aperture centres from the centres cut out of the background ring.
"""
from ..datatypes import Trace

# --- 孔径 / 背景 ---
from .aperture import (ApertureSpec, Aperture, crop_centered,
                       sigma_per_molecule, radius_from_sigma, aperture_kernel,
                       clean_sigma)
from .background import (AnnulusBackground, ring_mask, clip_median,
                         clip_median_std, valid_mask, ring_budget,
                         moving_average, model_weights, emitter_frame_photons,
                         emitter_sigma)

# --- 逐字段填充框架 ---
from .fill import (TRACE_FIELDS, FillContext, TraceBox, TraceFiller,
                   TraceFillPipeline, as_centers_yx, to_table, register_field)
from .signal import ApertureRaw
from .annulus import AnnulusRing
from .calculate import IntensityFromRawBg, PoissonSNR, standard_fillers
from . import changepoint
from .states import (OFF, ON, STATE_FIELDS, ThresholdStates, ChangepointStates,
                     summarize_states, runs_of, select_traces)

# --- 可视化 ---
from .viz import plot_trace_fields, plot_state_summary, show_trace_movies

__all__ = [
    # 数据类型
    "Trace",
    # 孔径 / 背景
    "ApertureSpec", "Aperture", "crop_centered", "sigma_per_molecule",
    "radius_from_sigma", "aperture_kernel", "clean_sigma",
    "AnnulusBackground", "ring_mask", "clip_median", "clip_median_std",
    "valid_mask", "ring_budget", "moving_average", "model_weights",
    "emitter_frame_photons", "emitter_sigma",
    # 逐字段填充
    "TRACE_FIELDS", "FillContext", "TraceBox", "TraceFiller",
    "TraceFillPipeline", "ApertureRaw", "AnnulusRing", "as_centers_yx",
    "to_table", "register_field", "IntensityFromRawBg", "PoissonSNR",
    "standard_fillers",
    # 状态判定
    "ThresholdStates", "ChangepointStates", "summarize_states",
    "select_traces", "runs_of", "STATE_FIELDS", "changepoint", "OFF", "ON",
    # 可视化
    "plot_trace_fields", "plot_state_summary", "show_trace_movies",
]
