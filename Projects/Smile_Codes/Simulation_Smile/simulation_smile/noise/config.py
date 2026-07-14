"""Noise configuration."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class NoiseConfig:
    """Noise-model parameters.

    Attributes
    ----------
    read_sigma : float
        Gaussian read-noise std in electrons.
    read_nonuniform : bool
        If True, draw a spatially varying read-noise std map (sCMOS-like).
    dark_current : float
        Dark electrons per pixel per second.
    background : float
        Uniform background photons per pixel per frame.
    background_kind : {"uniform", "nonuniform"}
        Uniform offset vs. a smooth rotated-super-Gaussian autofluorescence map.
    """

    read_sigma: float = 1.5
    read_nonuniform: bool = False
    dark_current: float = 0.0
    background: float = 5.0
    background_kind: Literal["uniform", "nonuniform"] = "uniform"
