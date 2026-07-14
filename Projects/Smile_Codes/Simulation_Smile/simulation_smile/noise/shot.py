"""Shot (Poisson) noise stage.

SKELETON -- signatures and docstrings only.

Reference: ``poisson_noise_numpy`` in DeepSTORM3D's ``physics_utils.py``.
"""
from __future__ import annotations

from ..core.base import SimContext
from ..core.registry import register


@register("noise", "poisson")
class Poisson:
    """Replace ctx.image with a Poisson draw (photon shot noise).

    Should run on the photon-domain image (before/around camera gain). Draws
    from ``ctx.rng`` for reproducibility.
    """

    name = "poisson"

    def apply(self, ctx: SimContext) -> SimContext:
        raise NotImplementedError
