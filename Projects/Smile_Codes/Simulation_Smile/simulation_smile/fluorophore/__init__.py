"""Fluorophore module: emitter geometry + photophysics (blinking/bleaching)."""
from __future__ import annotations

from .config import EmitterLayoutConfig, FluorophoreConfig
from .emitters import make_emitters
from .photophysics import BLEACHED, OFF, ON, AlwaysOn, BlinkingModel

__all__ = [
    "FluorophoreConfig",
    "EmitterLayoutConfig",
    "make_emitters",
    "BlinkingModel",
    "AlwaysOn",
    "OFF",
    "ON",
    "BLEACHED",
]
