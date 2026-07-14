"""Zernike polynomials for wavefront aberrations.

SKELETON -- signatures and docstrings only.

Reference: ``Codes_revised_JL/+gen/gen_zern_coeffs.m`` for the indexing
convention used elsewhere in this repo.
"""
from __future__ import annotations

import numpy as np


def noll_to_nm(j: int) -> tuple[int, int]:
    """Map a Noll single index ``j`` to radial/azimuthal orders ``(n, m)``.

    TODO: implement the standard Noll ordering.
    """
    raise NotImplementedError


def zernike_radial(n: int, m: int, rho: np.ndarray) -> np.ndarray:
    """Radial Zernike polynomial ``R_n^m(rho)`` on the unit disk."""
    raise NotImplementedError


def zernike(j: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Zernike mode ``j`` (Noll) sampled on polar grid (rho, theta).

    ``rho`` should be normalised so the pupil edge is at rho=1.
    """
    raise NotImplementedError


def wavefront(
    coeffs: dict[int, float], rho: np.ndarray, theta: np.ndarray
) -> np.ndarray:
    """Sum of Zernike modes weighted by ``coeffs`` (Noll index -> waves).

    Returns the aberration phase in radians: ``2*pi * sum_j c_j Z_j``.
    Used by :class:`..pupil.Pupil` to build the aberrated pupil phase.
    """
    raise NotImplementedError
