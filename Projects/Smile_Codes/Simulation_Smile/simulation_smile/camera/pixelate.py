"""Pixelation: bin the oversampled optical field down to camera pixels.

SKELETON -- signatures and docstrings only.

Reference: ``Codes_revised_JL/+gen/downsample2D.m``.
"""
from __future__ import annotations

import numpy as np

from ..core.base import Grid


def downsample2d(field: np.ndarray, factor: int) -> np.ndarray:
    """Sum-bin a fine-grid image by integer ``factor`` (energy-conserving).

    Equivalent to integrating photons over each camera pixel footprint.
    """
    raise NotImplementedError


def pixelate(field: np.ndarray, grid: Grid) -> np.ndarray:
    """Bring an oversampled optical field to the camera pixel grid.

    Uses ``grid.oversampling``; no-op when oversampling == 1.
    """
    raise NotImplementedError
