"""Visualisation helpers for debugging in notebooks.

SKELETON -- signatures and docstrings only. All functions accept host numpy
arrays (use ``core.backend.to_numpy`` first if on GPU).
"""
from __future__ import annotations

import numpy as np


def show_psf(psf: np.ndarray, ax=None):
    """Display a 2D PSF (or central slice of a 3D stack) with a colourbar."""
    raise NotImplementedError


def show_pupil(pupil: np.ndarray, ax=None):
    """Display the pupil amplitude and phase side by side."""
    raise NotImplementedError


def show_frame(image: np.ndarray, ax=None):
    """Display a single rendered camera frame."""
    raise NotImplementedError


def animate_movie(movie: np.ndarray, interval: int = 100):
    """Return a matplotlib animation of a (T, H, W) movie."""
    raise NotImplementedError


def plot_intensity_trace(trace: np.ndarray, dt: float, ax=None):
    """Plot a per-emitter ON/OFF or photon-vs-time trace (blinking debug)."""
    raise NotImplementedError
