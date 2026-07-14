"""Lens-level helpers: magnification, effective pixel size, phase masks.

SKELETON -- signatures and docstrings only.

Bridges :class:`..config.LensConfig` and the pupil/PSF, and provides a hook for
engineered PSFs (e.g. double-helix) via an SLM phase mask added to the pupil.
"""
from __future__ import annotations

import numpy as np

from .config import LensConfig, PSFConfig


def effective_pixel_size(camera_pixel: float, lens: LensConfig) -> float:
    """Camera pixel size projected back into sample space.

    ``camera_pixel / magnification``. Used to build the :class:`~..core.base.Grid`.
    """
    raise NotImplementedError


def back_focal_plane_diameter(lens: LensConfig) -> float:
    """Physical BFP diameter (metres). Port of the BFPdiam term in calc_bfp_grids."""
    raise NotImplementedError


class PhaseMask:
    """Optional engineered phase mask added to the pupil (e.g. SLM pattern).

    Lets the optics produce engineered PSFs. The mask phase (radians) is added
    to the pupil phase before the FFT.
    """

    def __init__(self, phase: np.ndarray | None = None):
        self.phase = phase

    @classmethod
    def double_helix(cls, cfg: PSFConfig) -> "PhaseMask":
        """Construct a double-helix phase mask. TODO."""
        raise NotImplementedError

    def phase_radians(self) -> np.ndarray:
        """Return the additive pupil phase in radians."""
        raise NotImplementedError
