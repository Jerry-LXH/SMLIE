"""Emitter geometry generation.

SKELETON -- signatures and docstrings only.

Reference: ``Codes_revised_JL/+gen/gen_random_points.m`` for random sampling.
"""
from __future__ import annotations

import numpy as np

from ..core.base import Emitters, Grid
from .config import EmitterLayoutConfig, FluorophoreConfig


def make_emitters(
    layout: EmitterLayoutConfig,
    fluor: FluorophoreConfig,
    grid: Grid,
    rng: np.random.Generator,
) -> Emitters:
    """Build an :class:`Emitters` set per the layout config.

    Dispatches on ``layout.kind`` to the samplers below and seeds each emitter's
    baseline photon budget from ``fluor.photon_rate``.
    """
    raise NotImplementedError


def random_positions(
    n: int, grid: Grid, z_range: tuple[float, float], rng: np.random.Generator
) -> np.ndarray:
    """Uniform-random xyz within the FOV and z-range. Shape (n, 3)."""
    raise NotImplementedError


def grid_positions(n: int, grid: Grid) -> np.ndarray:
    """Regular lattice of positions (for PSF-calibration-style scenes)."""
    raise NotImplementedError


def structure_positions(name: str, grid: Grid, rng: np.random.Generator) -> np.ndarray:
    """Positions sampled along a named structure (e.g. nuclear pore, line)."""
    raise NotImplementedError
