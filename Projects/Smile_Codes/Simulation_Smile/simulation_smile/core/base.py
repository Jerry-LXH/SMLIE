"""Core data structures and the ``Stage`` protocol.

This module defines the *contract* every simulation component obeys. It is the
foundation that makes components composable, so unlike the per-stage modules
(which are skeletons for now) this file is fully implemented.

Mental model
------------
A simulation is a sequence of :class:`Stage` objects applied to a single mutable
:class:`SimContext`. Each stage reads the fields it cares about and writes the
fields it produces:

    emitters  --(optics)-->  field/psf
    field     --(camera)-->  image (photons)
    image     --(noise)-->   image (ADU / electrons)

Because every stage shares the same ``ctx`` object, stages stay decoupled: you
can reorder, drop, or swap any of them as long as the field it needs is present.
Intermediate products are stashed in ``ctx.meta`` for inspection in notebooks.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Protocol, runtime_checkable

import numpy as np


# ---------------------------------------------------------------------------
# Sampling grid
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Grid:
    """Pixel / physical sampling grid for the field of view.

    Parameters
    ----------
    pixel_size : float
        Physical size of one camera pixel projected into sample space, in
        metres (e.g. ``108e-9`` for a 108 nm effective pixel).
    size : int
        Number of pixels along each side of the (square) field of view.
    oversampling : int
        Optical-domain oversampling factor. The PSF/field is computed on a
        finer grid (``size * oversampling``) and binned down to camera pixels
        to respect the sampling theorem; see :mod:`..core.coords`.
    """

    pixel_size: float
    size: int
    oversampling: int = 1

    @property
    def fine_size(self) -> int:
        """Number of samples per side on the oversampled optical grid."""
        return self.size * self.oversampling

    @property
    def fine_pixel_size(self) -> float:
        """Sample spacing on the oversampled optical grid, in metres."""
        return self.pixel_size / self.oversampling

    @property
    def fov(self) -> float:
        """Field-of-view side length, in metres."""
        return self.pixel_size * self.size


# ---------------------------------------------------------------------------
# Emitters (fluorescent point sources)
# ---------------------------------------------------------------------------
@dataclass
class Emitters:
    """A set of point emitters and their per-emitter state.

    Arrays are parallel (one row per emitter). Positions are in metres in
    sample space, with the grid centre at the origin.

    Attributes
    ----------
    xyz : ndarray, shape (N, 3)
        Emitter positions [x, y, z] in metres.
    photons : ndarray, shape (N,)
        Expected detected photons per frame when fully ON.
    state : ndarray, shape (N,)
        Photophysical state code per emitter (see ``fluorophore.photophysics``:
        0=OFF, 1=ON, 2=BLEACHED). Filled by the photophysics stage.
    ids : ndarray, shape (N,), optional
        Stable identifiers for tracking emitters across frames.
    """

    xyz: np.ndarray
    photons: np.ndarray
    state: Optional[np.ndarray] = None
    ids: Optional[np.ndarray] = None

    def __len__(self) -> int:
        return int(self.xyz.shape[0])


# ---------------------------------------------------------------------------
# The context that flows through the pipeline
# ---------------------------------------------------------------------------
@dataclass
class SimContext:
    """Mutable state passed between stages.

    A stage receives a ``SimContext``, mutates the fields it owns, and returns
    it (returning is by convention so chains read nicely; mutation is in place).

    Fields
    ------
    grid : Grid
        Sampling geometry. Always present.
    rng : numpy.random.Generator
        The single source of randomness. Stages MUST draw from this (never from
        the global numpy RNG) so a fixed seed makes a run fully reproducible.
    emitters : Emitters | None
        Current emitters. Set by the scene / fluorophore stages.
    field : ndarray | None
        Optical-domain quantity (PSF stack or accumulated intensity on the fine
        grid). Produced by the optics stage.
    image : ndarray | None
        Detector-domain image (photons, then electrons/ADU after camera+noise).
    frame_index : int
        Index of the frame currently being rendered (for movies).
    time : float
        Absolute time of the current frame, in seconds.
    meta : dict
        Free-form scratch space for intermediate products (e.g. ``meta['psf']``,
        ``meta['background']``) to support notebook visualisation/debugging.
    """

    grid: Grid
    rng: np.random.Generator
    emitters: Optional[Emitters] = None
    field: Optional[np.ndarray] = None
    image: Optional[np.ndarray] = None
    frame_index: int = 0
    time: float = 0.0
    meta: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# The Stage protocol -- the universal component interface
# ---------------------------------------------------------------------------
@runtime_checkable
class Stage(Protocol):
    """Every imaging component implements this single method.

    Implementations are typically constructed from a config dataclass and hold
    no per-frame state (or reset it via ``reset``). ``apply`` is called once per
    frame with the shared context.
    """

    def apply(self, ctx: SimContext) -> SimContext:
        """Read/transform fields on ``ctx`` and return it."""
        ...


class BaseStage:
    """Optional convenience base class for stages.

    Provides a ``name`` (used by the registry / logging) and a no-op ``reset``.
    Subclasses implement :meth:`apply`. Using this base is not required -- any
    object with an ``apply(ctx)`` method satisfies :class:`Stage`.
    """

    #: Short identifier; defaults to the class name. Used by the registry.
    name: str = "stage"

    def reset(self, ctx: SimContext) -> None:
        """Reset any per-run internal state. Called once before a movie starts."""

    def apply(self, ctx: SimContext) -> SimContext:  # pragma: no cover - abstract
        raise NotImplementedError
