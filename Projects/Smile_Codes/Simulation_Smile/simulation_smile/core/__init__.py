"""Core layer: foundational interfaces shared by every component.

Re-exports the most commonly used names so other modules and notebooks can do
``from simulation_smile.core import SimContext, Stage, Grid``.
"""
from __future__ import annotations

from .backend import get_xp, set_backend, to_numpy, xp
from .base import BaseStage, Emitters, Grid, SimContext, Stage
from .coords import check_nyquist, nyquist_pixel_size, radial_grid, sample_grid
from .registry import available, categories, get, register

__all__ = [
    # backend
    "xp",
    "get_xp",
    "set_backend",
    "to_numpy",
    # base
    "Grid",
    "Emitters",
    "SimContext",
    "Stage",
    "BaseStage",
    # coords
    "sample_grid",
    "radial_grid",
    "nyquist_pixel_size",
    "check_nyquist",
    # registry
    "register",
    "get",
    "available",
    "categories",
]
