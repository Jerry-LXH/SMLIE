"""Illumination module: excitation light spatial profiles."""
from __future__ import annotations

from .config import IlluminationConfig
from .profiles import GaussianSpot, TIRF, UniformField

__all__ = ["IlluminationConfig", "UniformField", "GaussianSpot", "TIRF"]
