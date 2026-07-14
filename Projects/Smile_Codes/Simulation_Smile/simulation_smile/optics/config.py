"""Optics configuration dataclasses (lens + PSF parameters).

Pure data, no behaviour. These are the knobs a notebook tunes; implementation
classes in :mod:`psf`, :mod:`pupil`, :mod:`lens` consume them.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class LensConfig:
    """Objective / relay-optics parameters.

    Attributes
    ----------
    na : float
        Numerical aperture of the objective.
    magnification : float
        Total system magnification (sample -> camera).
    n_immersion : float
        Refractive index of the immersion medium (e.g. 1.515 oil).
    n_sample : float
        Refractive index of the sample medium (e.g. 1.33 water). A mismatch
        with ``n_immersion`` introduces depth-dependent spherical aberration.
    tube_length : float
        Tube/relay focal length in metres (used for back-focal-plane scaling).
    """

    na: float = 1.49
    magnification: float = 100.0
    n_immersion: float = 1.515
    n_sample: float = 1.33
    tube_length: float = 0.2


@dataclass(frozen=True)
class PSFConfig:
    """Point-spread-function parameters.

    Attributes
    ----------
    wavelength : float
        Emission wavelength in metres.
    na : float
        Numerical aperture (mirrors LensConfig.na; kept here so a PSF can be
        built standalone for debugging).
    n_immersion, n_sample : float
        Refractive indices, as in LensConfig.
    zernike : dict[int, float]
        Aberration coefficients keyed by Noll index, in waves (RMS). Empty =
        unaberrated. e.g. ``{5: 0.3}`` adds 0.3 waves of astigmatism.
    nominal_focal_plane : float
        Axial position of the nominal focal plane in metres (z=0 reference).
    sigma_gaussian : float
        Std of the Gaussian PSF model, in metres (only used by GaussianPSF).
    """

    wavelength: float = 668e-9
    na: float = 1.49
    n_immersion: float = 1.515
    n_sample: float = 1.33
    zernike: dict[int, float] = field(default_factory=dict)
    nominal_focal_plane: float = 0.0
    sigma_gaussian: float = 130e-9
