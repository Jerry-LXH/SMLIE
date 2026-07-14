"""Per-emitter brightness / lifetime / (optional) spectral properties.

SKELETON -- signatures and docstrings only. Optional refinement layer on top of
the photophysics state machine; can be skipped for basic simulations.
"""
from __future__ import annotations

import numpy as np

from .config import FluorophoreConfig


def sample_burst_photons(
    n: int, cfg: FluorophoreConfig, rng: np.random.Generator
) -> np.ndarray:
    """Draw per-burst photon counts with mean ``photon_rate`` and CV ``photon_cv``.

    Models shot-to-shot brightness variability (e.g. log-normal).
    """
    raise NotImplementedError
