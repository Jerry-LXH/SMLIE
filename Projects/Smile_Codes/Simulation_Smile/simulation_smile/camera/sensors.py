"""Camera sensor models as pipeline stages.

SKELETON -- signatures and docstrings only.

A camera stage takes the optical-domain photon image in ``ctx.field``,
pixelates it to the camera grid, applies QE / gain / baseline, and writes the
detector image to ``ctx.image`` (in electrons or ADU). Sensor-specific noise
(EM gain stochasticity, read noise) is delegated to the :mod:`..noise` stages
so the chain stays modular, but the deterministic conversion lives here.
"""
from __future__ import annotations

import numpy as np

from ..core.base import SimContext
from ..core.registry import register
from .config import CameraConfig


class _CameraBase:
    def __init__(self, cfg: CameraConfig):
        self.cfg = cfg

    def expose(self, photon_image: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """Convert a photon image to a detector image (electrons/ADU)."""
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        """Pixelate ctx.field, run :meth:`expose`, write ctx.image."""
        raise NotImplementedError


@register("camera", "ideal")
class IdealCamera(_CameraBase):
    """Unity QE, no gain, no quantisation -- photon image straight through."""

    name = "ideal"


@register("camera", "emccd")
class EMCCD(_CameraBase):
    """Electron-multiplying CCD: stochastic EM gain (excess-noise factor ~sqrt(2))."""

    name = "emccd"


@register("camera", "scmos")
class SCMOSCamera(_CameraBase):
    """sCMOS: per-pixel gain/offset, dominated by pixel-dependent read noise."""

    name = "scmos"
