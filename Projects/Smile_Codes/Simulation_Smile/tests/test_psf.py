"""PSF correctness tests (skeleton).

These assert physical properties the implementation must satisfy. Marked xfail
until the optics module is filled in.
"""
import numpy as np
import pytest

from simulation_smile import configs as C
from simulation_smile.optics import PupilPSF


@pytest.mark.xfail(reason="PupilPSF not implemented yet", strict=False)
def test_psf_unit_sum():
    """A normalised PSF integrates to 1 (energy conservation)."""
    grid = C.Grid(pixel_size=108e-9, size=64, oversampling=4)
    psf = PupilPSF(C.PSFConfig(na=1.49), grid=grid).psf(z=0.0)
    assert np.isclose(psf.sum(), 1.0, rtol=1e-3)


@pytest.mark.xfail(reason="PupilPSF not implemented yet", strict=False)
def test_psf_symmetric_unaberrated():
    """An unaberrated in-focus PSF is (approximately) centro-symmetric."""
    grid = C.Grid(pixel_size=108e-9, size=64, oversampling=4)
    psf = PupilPSF(C.PSFConfig(na=1.49, zernike={}), grid=grid).psf(z=0.0)
    assert np.allclose(psf, psf[::-1, ::-1], atol=1e-6)


@pytest.mark.xfail(reason="PupilPSF not implemented yet", strict=False)
def test_psf_fwhm_matches_diffraction_limit():
    """FWHM ~ 0.51 * lambda / NA (Abbe-ish sanity check, within tolerance)."""
    grid = C.Grid(pixel_size=20e-9, size=128, oversampling=2)
    cfg = C.PSFConfig(na=1.0, wavelength=500e-9)
    _ = PupilPSF(cfg, grid=grid).psf(z=0.0)
    # TODO: measure FWHM and compare to ~0.51*lambda/NA.
