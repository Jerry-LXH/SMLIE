"""Camera module: pixelation + sensor models (ideal / EMCCD / sCMOS)."""
from __future__ import annotations

from .config import CameraConfig
from .pixelate import downsample2d, pixelate
from .sensors import EMCCD, IdealCamera, SCMOSCamera

__all__ = [
    "CameraConfig",
    "pixelate",
    "downsample2d",
    "IdealCamera",
    "EMCCD",
    "SCMOSCamera",
]
