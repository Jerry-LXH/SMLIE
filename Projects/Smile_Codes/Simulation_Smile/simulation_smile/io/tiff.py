"""Multi-frame TIFF I/O.

SKELETON -- signatures and docstrings only. Thin wrapper over ``tifffile``.
"""
from __future__ import annotations

import numpy as np


def save_tiff(path: str, movie: np.ndarray, dtype="uint16") -> None:
    """Write a (T, H, W) stack as a multi-page TIFF (ImageJ-compatible)."""
    raise NotImplementedError


def load_tiff(path: str) -> np.ndarray:
    """Read a multi-page TIFF into a (T, H, W) array."""
    raise NotImplementedError
