"""Illumination configuration."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class IlluminationConfig:
    """Excitation light spatial profile.

    Attributes
    ----------
    kind : {"uniform", "gaussian", "tirf"}
        Profile shape. Selects the implementation in :mod:`profiles`.
    power : float
        Relative excitation intensity scale (multiplies emitter photon rates).
    waist : float
        1/e^2 Gaussian beam waist radius in metres (gaussian/tirf).
    center : tuple[float, float]
        Beam centre offset (x, y) in metres from the FOV centre.
    tirf_depth : float
        TIRF evanescent-field 1/e penetration depth in metres (tirf only).
    """

    kind: Literal["uniform", "gaussian", "tirf"] = "uniform"
    power: float = 1.0
    waist: float = 20e-6
    center: tuple[float, float] = (0.0, 0.0)
    tirf_depth: float = 100e-9
