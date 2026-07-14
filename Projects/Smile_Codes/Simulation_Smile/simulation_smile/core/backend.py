"""Numerical backend abstraction.

All numerical code in the package should import the array module through this
module (``from simulation_smile.core.backend import xp``) instead of importing
numpy directly. This lets the whole simulation switch between CPU (numpy) and
GPU (cupy) with a single call to :func:`set_backend`.

The default backend is numpy. GPU backends are entirely optional -- importing
this module never requires cupy/torch to be installed.

Usage
-----
>>> from simulation_smile.core import backend
>>> backend.set_backend("cupy")      # opt into GPU
>>> xp = backend.get_xp()
>>> a = xp.zeros((64, 64))
"""
from __future__ import annotations

from typing import Any, Literal

import numpy as _np

BackendName = Literal["numpy", "cupy"]

# Module-level handle to the active array module. Code may either call
# ``get_xp()`` each time, or import the ``xp`` name and re-read it. We keep a
# mutable module global so ``set_backend`` is globally effective.
_active_name: BackendName = "numpy"
xp: Any = _np


def set_backend(name: BackendName) -> Any:
    """Select the active array backend and return its array module.

    Parameters
    ----------
    name : {"numpy", "cupy"}
        Backend to activate. "cupy" requires the optional cupy install.

    Notes
    -----
    Modules that captured ``xp`` at import time should re-read
    ``backend.xp`` (or call :func:`get_xp`) after switching.
    """
    global xp, _active_name
    if name == "numpy":
        xp = _np
    elif name == "cupy":  # pragma: no cover - optional dependency
        import cupy  # noqa: F401  (imported lazily on purpose)

        xp = cupy
    else:
        raise ValueError(f"Unknown backend: {name!r}")
    _active_name = name
    return xp


def get_xp() -> Any:
    """Return the active array module (numpy or cupy)."""
    return xp


def get_backend_name() -> BackendName:
    """Return the name of the active backend."""
    return _active_name


def to_numpy(a: Any) -> _np.ndarray:
    """Bring an array back to host numpy regardless of the active backend.

    Useful before plotting / saving, where matplotlib and tifffile expect
    host arrays.
    """
    if isinstance(a, _np.ndarray):
        return a
    # cupy arrays expose .get(); fall back to np.asarray otherwise.
    get = getattr(a, "get", None)
    if callable(get):  # pragma: no cover - optional dependency path
        return get()
    return _np.asarray(a)
