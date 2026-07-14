"""Background and dark-current stages.

SKELETON -- signatures and docstrings only.

Reference: ``NonUniformBg`` in DeepSTORM3D's ``physics_utils.py`` for the
rotated super-Gaussian autofluorescence model.
"""
from __future__ import annotations

import numpy as np

from ..core.base import SimContext
from ..core.registry import register
from .config import NoiseConfig


@register("noise", "background")
class Background:
    """Add a uniform or non-uniform background photon offset.

    Runs in the photon domain before Poisson noise so background photons are
    themselves shot-noise limited.
    """

    name = "background"

    def __init__(self, cfg: NoiseConfig):
        self.cfg = cfg
        self._bg_map: np.ndarray | None = None

    def background_map(self, ctx: SimContext) -> np.ndarray:
        """Uniform value or a smooth rotated super-Gaussian map."""
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError


@register("noise", "dark_current")
class DarkCurrent:
    """Add dark electrons proportional to exposure time."""

    name = "dark_current"

    def __init__(self, cfg: NoiseConfig):
        self.cfg = cfg

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError
