"""Aggregated config namespace.

Import as ``from simulation_smile import configs as C`` and access every config
dataclass plus :class:`~core.base.Grid` from one place:

>>> C.PSFConfig(NA=1.49)
>>> C.Grid(pixel_size=108e-9, size=64)
"""
from __future__ import annotations

from .camera.config import CameraConfig
from .core.base import Grid
from .fluorophore.config import EmitterLayoutConfig, FluorophoreConfig
from .illumination.config import IlluminationConfig
from .noise.config import NoiseConfig
from .optics.config import LensConfig, PSFConfig

__all__ = [
    "Grid",
    "LensConfig",
    "PSFConfig",
    "IlluminationConfig",
    "FluorophoreConfig",
    "EmitterLayoutConfig",
    "CameraConfig",
    "NoiseConfig",
]
