"""Fluorophore photophysics and emitter-geometry configuration."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class FluorophoreConfig:
    """Photophysical rate constants for a blinking/bleaching dye.

    Rates are in 1/s. The state machine is ON <-> OFF -> BLEACHED (bleaching is
    absorbing). Per-frame photon counts integrate the ON time over the exposure.

    Attributes
    ----------
    k_on : float
        OFF -> ON switching rate (recovery from dark state).
    k_off : float
        ON -> OFF switching rate (entry into dark state).
    k_bleach : float
        ON -> BLEACHED irreversible photobleaching rate.
    photon_rate : float
        Detected photons per second while ON (before illumination scaling).
    photon_cv : float
        Coefficient of variation for per-burst photon count (shot-to-shot
        brightness variability), 0 = deterministic.
    initial_state : {"off", "on"}
        State all emitters start in.
    """

    k_on: float = 0.5
    k_off: float = 20.0
    k_bleach: float = 0.01
    photon_rate: float = 5000.0
    photon_cv: float = 0.0
    initial_state: Literal["off", "on"] = "off"


@dataclass(frozen=True)
class EmitterLayoutConfig:
    """How emitter positions are generated for a scene.

    Attributes
    ----------
    kind : {"random", "grid", "from_array", "structure"}
        Sampling strategy (see :mod:`emitters`).
    n : int
        Number of emitters (random/grid).
    density : float
        Alternative to ``n``: emitters per square metre (random).
    z_range : tuple[float, float]
        Axial extent in metres for 3D scenes.
    structure : str
        Named structure for ``kind="structure"`` (e.g. "nuclear_pore", "line").
    """

    kind: Literal["random", "grid", "from_array", "structure"] = "random"
    n: int = 200
    density: float = 0.0
    z_range: tuple[float, float] = (0.0, 0.0)
    structure: str = ""
