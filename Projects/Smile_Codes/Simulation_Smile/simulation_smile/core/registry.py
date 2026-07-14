"""Component registry: look up implementations by name.

Enables config-driven assembly and parameter sweeps -- e.g. a YAML file can say
``optics: {kind: pupil, ...}`` and the pipeline builder resolves "pupil" to the
:class:`PupilPSF` class via :func:`get`.

Each category (optics, illumination, photophysics, camera, noise, ...) has its
own namespace so the same short name can mean different things per category.

Usage
-----
>>> from simulation_smile.core.registry import register, get

>>> @register("optics", "pupil")
... class PupilPSF:
...     ...

>>> cls = get("optics", "pupil")
"""
from __future__ import annotations

from typing import Callable, Type, TypeVar

T = TypeVar("T")

# category -> {name -> class}
_REGISTRY: dict[str, dict[str, type]] = {}


def register(category: str, name: str) -> Callable[[Type[T]], Type[T]]:
    """Class decorator that records ``cls`` under ``(category, name)``."""

    def _decorator(cls: Type[T]) -> Type[T]:
        _REGISTRY.setdefault(category, {})[name] = cls
        return cls

    return _decorator


def get(category: str, name: str) -> type:
    """Return the class registered under ``(category, name)``."""
    try:
        return _REGISTRY[category][name]
    except KeyError as exc:  # pragma: no cover - trivial
        available = sorted(_REGISTRY.get(category, {}))
        raise KeyError(
            f"No '{name}' registered under category '{category}'. "
            f"Available: {available}"
        ) from exc


def available(category: str) -> list[str]:
    """List registered names for a category."""
    return sorted(_REGISTRY.get(category, {}))


def categories() -> list[str]:
    """List all known categories."""
    return sorted(_REGISTRY)
