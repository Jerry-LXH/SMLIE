"""PSF models and the optics ``Stage``.

SKELETON -- signatures and docstrings only (except trivial glue).

Two interchangeable PSF implementations share a common protocol so the pipeline
can swap them freely:

* :class:`PupilPSF`  -- Fourier optics, ``PSF = |FFT(pupil)|**2``. Physically
  faithful, supports aberrations and 3D. Port of ``gen_PSF3D.m`` /
  ``MaskPhasesToPSFs`` in DeepSTORM3D's ``physics_utils.py``.
* :class:`GaussianPSF` -- fast 2D/3D Gaussian approximation for quick prototyping.

Both act as a :class:`~..core.base.Stage`: they read ``ctx.emitters`` and write
the accumulated optical-domain intensity into ``ctx.field`` (and stash the raw
PSF in ``ctx.meta['psf']``).
"""
from __future__ import annotations

from typing import Protocol

import numpy as np

from ..core.base import Grid, SimContext
from ..core.registry import register
from .config import PSFConfig
from .pupil import Pupil


class PSFModel(Protocol):
    """Common interface for all PSF models."""

    def psf(self, z: float = 0.0) -> np.ndarray:
        """Single 2D PSF at axial position ``z`` (normalised to unit sum)."""
        ...

    def psf_stack(self, z: np.ndarray) -> np.ndarray:
        """3D PSF stack, one 2D PSF per entry of ``z``. Shape (len(z), H, W)."""
        ...

    def render(self, emitters, grid: Grid) -> np.ndarray:
        """Accumulate PSFs for all ON emitters into one intensity image."""
        ...


@register("optics", "pupil")
class PupilPSF:
    """Fourier-optics PSF from a (possibly aberrated) pupil function."""

    name = "pupil"

    def __init__(self, cfg: PSFConfig, grid: Grid | None = None):
        self.cfg = cfg
        self.grid = grid
        self._pupil = Pupil(grid, cfg) if grid is not None else None

    def psf(self, z: float = 0.0) -> np.ndarray:
        """``|FFT(pupil(z))|**2`` normalised to unit sum.

        TODO: port the FFT body of ``gen_PSF3D.m``.
        """
        raise NotImplementedError

    def psf_stack(self, z: np.ndarray) -> np.ndarray:
        """Stack of 2D PSFs over the given z positions."""
        raise NotImplementedError

    def render(self, emitters, grid: Grid) -> np.ndarray:
        """Place a PSF at each ON emitter (sub-pixel shift via pupil phase ramp).

        TODO: port ``EmittersToPhases`` -> apply lateral/axial phases per emitter
        then sum, scaling each by its photon count.
        """
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        """Stage entry point: emitters -> ctx.field."""
        raise NotImplementedError


@register("optics", "gaussian")
class GaussianPSF:
    """Gaussian-approximation PSF (fast, no aberrations)."""

    name = "gaussian"

    def __init__(self, cfg: PSFConfig, grid: Grid | None = None):
        self.cfg = cfg
        self.grid = grid

    def psf(self, z: float = 0.0) -> np.ndarray:
        """2D Gaussian with sigma from cfg (optionally z-dependent defocus)."""
        raise NotImplementedError

    def psf_stack(self, z: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def render(self, emitters, grid: Grid) -> np.ndarray:
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError
