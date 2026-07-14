"""Config serialisation for reproducible experiments.

SKELETON -- signatures and docstrings only.

Dumps/loads the dataclass configs (optics/camera/noise/...) to YAML or JSON so a
notebook run can be saved alongside its output and replayed exactly.
"""
from __future__ import annotations

from typing import Any


def save_config(path: str, configs: dict[str, Any]) -> None:
    """Serialise a dict of named config dataclasses to YAML/JSON."""
    raise NotImplementedError


def load_config(path: str) -> dict[str, Any]:
    """Load configs back into their dataclasses (dispatch by registry name)."""
    raise NotImplementedError
