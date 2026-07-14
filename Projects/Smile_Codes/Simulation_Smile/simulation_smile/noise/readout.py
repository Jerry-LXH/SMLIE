"""Gaussian read-out noise stage.

SKELETON -- signatures and docstrings only.

Reference: ``ReadNoise`` in DeepSTORM3D's ``physics_utils.py`` (uniform baseline
plus uniform or spatially varying std).
"""
from __future__ import annotations

import numpy as np

from ..core.base import SimContext
from ..core.registry import register
from .config import NoiseConfig


@register("noise", "readout")
class ReadNoise:
    """Add Gaussian read noise (electrons), optionally spatially non-uniform."""

    name = "readout"

    def __init__(self, cfg: NoiseConfig):
        self.cfg = cfg
        self._std_map: np.ndarray | None = None  # cached non-uniform std map

    def std_map(self, ctx: SimContext) -> np.ndarray:
        """Return the per-pixel read-noise std (scalar broadcast or a map)."""
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError
