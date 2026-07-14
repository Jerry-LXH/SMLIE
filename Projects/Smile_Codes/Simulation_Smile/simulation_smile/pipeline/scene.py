"""Scene definition: emitters + illumination + acquisition settings.

SKELETON -- signatures and docstrings only.

A :class:`Scene` bundles everything the pipeline needs that is *not* an imaging
stage: the emitter set, the sampling grid, and the acquisition timing. The
pipeline holds the stages; the scene holds the subject.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.base import Emitters, Grid


@dataclass
class Scene:
    """The thing being imaged plus acquisition geometry/timing.

    Attributes
    ----------
    grid : Grid
        Sampling grid (FOV).
    emitters : Emitters
        Emitter set (positions + baseline photons + state).
    exposure : float
        Per-frame exposure time in seconds (drives the photophysics step).
    rng : numpy.random.Generator
        Master RNG (a fixed seed makes the whole acquisition reproducible).
    """

    grid: Grid
    emitters: Emitters
    exposure: float = 30e-3
    rng: np.random.Generator | None = None
