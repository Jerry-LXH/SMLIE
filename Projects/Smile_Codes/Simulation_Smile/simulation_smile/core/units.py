"""Physical units and constants used throughout the simulator.

Convention: **SI base units internally** -- lengths in metres, time in
seconds, wavelengths in metres. Helper converters are provided so that
notebooks can write the more natural nm / um / ms and convert at the boundary.

Keeping a single canonical unit system avoids the most common class of bugs in
optics simulations (mixing nm and um in the same expression).
"""
from __future__ import annotations

# --- length ---------------------------------------------------------------
NM = 1e-9   # metres per nanometre
UM = 1e-6   # metres per micrometre
MM = 1e-3   # metres per millimetre

# --- time -----------------------------------------------------------------
MS = 1e-3   # seconds per millisecond
US = 1e-6   # seconds per microsecond


def nm(x: float) -> float:
    """Nanometres -> metres."""
    return x * NM


def um(x: float) -> float:
    """Micrometres -> metres."""
    return x * UM


def ms(x: float) -> float:
    """Milliseconds -> seconds."""
    return x * MS


def to_nm(x_m: float) -> float:
    """Metres -> nanometres."""
    return x_m / NM


def to_um(x_m: float) -> float:
    """Metres -> micrometres."""
    return x_m / UM
