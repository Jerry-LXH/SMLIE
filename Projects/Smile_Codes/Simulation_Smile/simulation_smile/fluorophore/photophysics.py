"""Photophysical state machine: blinking and photobleaching.

SKELETON -- signatures and docstrings only.

This is the engine of the time series. Each emitter is a continuous-time Markov
process over states {OFF, ON, BLEACHED}. Per frame we advance the chain by the
exposure time and integrate the fraction of the exposure spent ON to get that
frame's photon count.

The stage reads/writes ``ctx.emitters.state`` and sets ``ctx.emitters.photons``
to the per-frame emitted photons for the current frame.
"""
from __future__ import annotations

import numpy as np

from ..core.base import Emitters, SimContext
from ..core.registry import register
from .config import FluorophoreConfig

# State codes (also used by Emitters.state)
OFF = 0
ON = 1
BLEACHED = 2


@register("photophysics", "blinking")
class BlinkingModel:
    """ON/OFF blinking with irreversible photobleaching.

    Parameters come from :class:`FluorophoreConfig`. Holds per-emitter state
    that persists across frames within a movie.
    """

    name = "blinking"

    def __init__(self, cfg: FluorophoreConfig):
        self.cfg = cfg
        self._state: np.ndarray | None = None  # per-emitter state vector

    def reset(self, ctx: SimContext) -> None:
        """Initialise per-emitter state from ``cfg.initial_state``.

        Called once before a movie starts (by the timeline engine).
        """
        raise NotImplementedError

    def step(self, dt: float, rng: np.random.Generator) -> np.ndarray:
        """Advance the chain by ``dt`` seconds; return per-emitter ON fraction.

        Uses the continuous-time transition probabilities (matrix exponential or
        exponential waiting times) for OFF<->ON and ON->BLEACHED. The returned
        ON fraction (0..1) times ``photon_rate*dt`` gives expected photons.
        """
        raise NotImplementedError

    def apply(self, ctx: SimContext) -> SimContext:
        """Advance one exposure and write per-frame photons into ctx.emitters.

        Uses ``ctx.meta['exposure']`` for dt and ``ctx.rng`` for randomness.
        """
        raise NotImplementedError


@register("photophysics", "always_on")
class AlwaysOn:
    """Trivial model: every emitter emits at its full rate every frame.

    Useful for static single-frame snapshots (no time dynamics).
    """

    name = "always_on"

    def __init__(self, cfg: FluorophoreConfig):
        self.cfg = cfg

    def reset(self, ctx: SimContext) -> None:  # no state
        pass

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError
