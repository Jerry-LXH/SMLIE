"""Illumination spatial profiles as pipeline stages.

SKELETON -- signatures and docstrings only.

An illumination stage scales each emitter's photon output by the local
excitation intensity at its position (and, for TIRF, its depth). It runs before
the optics stage so that ``ctx.emitters.photons`` already reflects how brightly
each molecule is excited.

Reference for Gaussian-spot fitting: ``Codes_revised_JL/+optics/+illumi``.
"""
from __future__ import annotations

import numpy as np

from ..core.base import SimContext
from ..core.registry import register
from .config import IlluminationConfig


class _IlluminationBase:
    def __init__(self, cfg: IlluminationConfig):
        self.cfg = cfg

    def intensity_at(self, xyz: np.ndarray) -> np.ndarray:
        """Relative excitation intensity (0..power) at each emitter position."""
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        """Scale ``ctx.emitters.photons`` by the local intensity."""
        raise NotImplementedError


@register("illumination", "uniform")
class UniformField(_IlluminationBase):
    """Flat illumination across the FOV (widefield)."""

    name = "uniform"


@register("illumination", "gaussian")
class GaussianSpot(_IlluminationBase):
    """Gaussian excitation spot (typical for focused laser illumination)."""

    name = "gaussian"


@register("illumination", "tirf")
class TIRF(_IlluminationBase):
    """Total-internal-reflection field with exponential depth decay."""

    name = "tirf"
