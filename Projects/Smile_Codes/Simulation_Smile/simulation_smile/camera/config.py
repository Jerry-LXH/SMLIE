"""Camera configuration."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class CameraConfig:
    """Detector parameters.

    Attributes
    ----------
    kind : {"ideal", "emccd", "scmos"}
        Sensor model (selects implementation in :mod:`sensors`).
    pixel_size : float
        Physical camera pixel size in metres (sensor space).
    qe : float
        Quantum efficiency (0..1): photons -> photoelectrons.
    gain : float
        EM gain (EMCCD) or analog gain. e- multiplication factor.
    adu_per_electron : float
        Conversion gain: ADU per photoelectron.
    baseline : float
        Constant offset added to every pixel, in ADU.
    full_well : float
        Full-well capacity in electrons (saturation/clipping).
    bit_depth : int
        ADC bit depth (for quantisation and max ADU).
    exposure : float
        Exposure (integration) time per frame in seconds.
    """

    kind: Literal["ideal", "emccd", "scmos"] = "scmos"
    pixel_size: float = 6.5e-6
    qe: float = 0.95
    gain: float = 1.0
    adu_per_electron: float = 0.5
    baseline: float = 100.0
    full_well: float = 30000.0
    bit_depth: int = 16
    exposure: float = 30e-3
