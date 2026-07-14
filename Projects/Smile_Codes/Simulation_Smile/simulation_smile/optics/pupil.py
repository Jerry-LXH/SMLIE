"""Pupil (back-focal-plane) function construction.

SKELETON -- signatures and docstrings only.

Port target: ``calc_bfp_grids`` in
``ML_for_Smile/Paper_Codes/DeepSTORM3D-master/DeepSTORM3D/physics_utils.py``
(torch) -> numpy via ``backend.xp``. Also mirrors
``Codes_revised_JL/+gen/def_pupilcoor.m`` and ``defocus_wavefront.m``.

The pupil function is::

    P(kx, ky) = A(kx, ky) * exp( i * [ phi_aberration + phi_defocus(z) ] )

where ``A`` is the NA-limited circular aperture, ``phi_aberration`` comes from
:mod:`zernike`, and ``phi_defocus`` is the axial propagation phase that depends
on the immersion/sample refractive indices.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.base import Grid
from .config import PSFConfig


@dataclass
class PupilGrids:
    """Pre-computed back-focal-plane grids (cached per (grid, config)).

    Attributes
    ----------
    aperture : ndarray
        NA-limited circular aperture mask (1 inside, 0 outside).
    defocus_const : ndarray
        Axial phase per unit z (radians/metre) -- multiply by z and exponentiate
        to get the defocus phase. Equivalent to ``defocus_wavefront.m``.
    rho, theta : ndarray
        Normalised polar coords on the pupil (rho=1 at the aperture edge), for
        feeding Zernike modes.
    """

    aperture: np.ndarray
    defocus_const: np.ndarray
    rho: np.ndarray
    theta: np.ndarray


class Pupil:
    """Builds and caches the pupil function for a given grid + PSF config."""

    def __init__(self, grid: Grid, cfg: PSFConfig):
        self.grid = grid
        self.cfg = cfg
        self._grids: PupilGrids | None = None

    def grids(self) -> PupilGrids:
        """Return (and lazily compute) the back-focal-plane grids.

        TODO: port ``calc_bfp_grids`` -- build the angular grid from NA and
        immersion index, the circular aperture, and the defocus constant.
        """
        raise NotImplementedError

    def aberration_phase(self) -> np.ndarray:
        """Static aberration phase from the config's Zernike coefficients.

        Delegates to :func:`..zernike.wavefront` on the pupil polar grid.
        """
        raise NotImplementedError

    def function(self, z: float = 0.0) -> np.ndarray:
        """Complex pupil function at axial position ``z`` (metres).

        ``P = aperture * exp(i*(aberration + z*defocus_const))``.
        """
        raise NotImplementedError
