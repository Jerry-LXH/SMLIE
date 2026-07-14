"""Preset scenes: fully-configured pipelines for common modalities.

SKELETON -- signatures and docstrings only.

Each preset wires up a :class:`~pipeline.pipeline.SimulationPipeline` with
sensible default configs so a notebook can start in one line:

>>> pipe = presets.smlm_storm()
>>> movie = pipe.render_movie(200)

Presets are just example assemblies -- copy one and tweak to build your own.
"""
from __future__ import annotations

from .pipeline.pipeline import SimulationPipeline


def smlm_storm() -> SimulationPipeline:
    """STORM/PALM: sparse blinking emitters, sCMOS, Poisson+read noise.

    TODO: assemble PupilPSF + BlinkingModel + UniformField + SCMOSCamera +
    [Background, Poisson, ReadNoise] over a Nyquist-sampled grid.
    """
    raise NotImplementedError


def widefield() -> SimulationPipeline:
    """Conventional widefield: all emitters on, uniform illumination."""
    raise NotImplementedError


def confocal() -> SimulationPipeline:
    """Point-scanning confocal approximation (Gaussian spot illumination)."""
    raise NotImplementedError
