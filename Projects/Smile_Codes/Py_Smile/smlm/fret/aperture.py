"""Deprecated location: aperture helpers now live in :mod:`smlm.trace.aperture`."""
from ..trace.aperture import (ApertureSpec, Aperture, crop_centered,  # noqa: F401
                              sigma_per_molecule)

__all__ = ['ApertureSpec', 'Aperture', 'crop_centered', 'sigma_per_molecule']
