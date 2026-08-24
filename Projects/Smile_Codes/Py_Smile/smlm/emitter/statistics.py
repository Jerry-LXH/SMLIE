# =====================================================================
#  统计 (向量化, 从源 Localizations 回查 photons/bg/sigma)
# =====================================================================
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from ..datatypes import Localizations, Emitter, EmitterStats

def collect_statistics(emitters, locs, ex_time, interval):
    N = len(emitters)
    one_frame = ex_time + interval
    ph  = None if locs.photons is None else np.asarray(locs.photons, float)
    sg  = None if locs.sigma   is None else np.asarray(locs.sigma, float)
    bg  = None if locs.bg      is None else np.asarray(locs.bg, float)
    g = lambda arr, e, f: f(arr[e.loc_idx]) if arr is not None else np.nan

    return EmitterStats(
        centroid       = np.array([e.centroid for e in emitters]) if N else np.empty((0, 2)),
        survival_sec   = np.array([e.span * one_frame for e in emitters]),
        n_locs         = np.array([e.n_locs for e in emitters]),
        duty_cycle     = np.array([e.duty_cycle for e in emitters]),
        bleached       = np.array([e.bleached for e in emitters]),
        eccentricity   = np.array([e.eccentricity() if e.n_locs >= 2 else np.nan
                                   for e in emitters]),
        photons_per_sec= np.array([g(ph, e, np.mean) / ex_time for e in emitters]),
        photons_sum    = np.array([g(ph, e, np.sum) for e in emitters]),
        sigma_mean     = np.array([g(sg, e, np.mean) for e in emitters]),
        bg_mean        = np.array([g(bg, e, np.mean) for e in emitters]),
    )
