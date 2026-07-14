"""Coordinate grids and sampling-theorem helpers.

Provides the cartesian/polar grids that the optics module needs (sample space
and Fourier/pupil space) plus a Nyquist check so notebooks get warned when the
camera pixel size under-samples the PSF.

These are foundational helpers used across modules, so they are implemented
(thin) rather than stubbed.
"""
from __future__ import annotations

import warnings

import numpy as np

from .base import Grid


def sample_grid(grid: Grid, fine: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Return centred cartesian coordinate arrays (X, Y) in metres.

    Parameters
    ----------
    fine : bool
        If True, use the oversampled optical grid (``grid.fine_size`` /
        ``grid.fine_pixel_size``); otherwise the camera-pixel grid.
    """
    n = grid.fine_size if fine else grid.size
    dx = grid.fine_pixel_size if fine else grid.pixel_size
    ax = (np.arange(n) - (n - 1) / 2.0) * dx
    return np.meshgrid(ax, ax)


def radial_grid(grid: Grid, fine: bool = True) -> np.ndarray:
    """Return the radial coordinate ``r = sqrt(X**2 + Y**2)`` in metres."""
    X, Y = sample_grid(grid, fine=fine)
    return np.hypot(X, Y)


def nyquist_pixel_size(wavelength: float, na: float) -> float:
    """Camera pixel size (sample space) that just satisfies Nyquist.

    The diffraction-limited cutoff frequency is ``2*NA/lambda``; Nyquist
    requires a sample spacing of ``lambda / (4*NA)``.
    """
    return wavelength / (4.0 * na)


def check_nyquist(grid: Grid, wavelength: float, na: float) -> bool:
    """Warn (and return False) if the camera grid under-samples the PSF.

    Returns True when ``grid.pixel_size <= lambda/(4*NA)``.
    """
    limit = nyquist_pixel_size(wavelength, na)
    ok = grid.pixel_size <= limit
    if not ok:
        warnings.warn(
            f"Camera pixel size {grid.pixel_size*1e9:.1f} nm exceeds the "
            f"Nyquist limit {limit*1e9:.1f} nm for NA={na}, "
            f"lambda={wavelength*1e9:.0f} nm. Increase oversampling or shrink "
            f"the pixel size to avoid aliasing.",
            stacklevel=2,
        )
    return ok
