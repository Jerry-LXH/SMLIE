"""simulation_smile -- a composable microscope-imaging simulator.

Architecture: **config + component + pipeline** (see README.md).

* Each imaging stage (optics/PSF, illumination, fluorophore photophysics,
  camera, noise) is a component implementing the :class:`~core.base.Stage`
  protocol and built from a config dataclass.
* :class:`~pipeline.pipeline.SimulationPipeline` composes stages in a canonical
  order and renders single frames or time-series movies.
* Components register themselves by name (:mod:`core.registry`) for
  config-driven assembly and parameter sweeps.

Quick start
-----------
>>> from simulation_smile import presets
>>> movie = presets.smlm_storm().render_movie(200)

Or assemble manually:

>>> from simulation_smile import SimulationPipeline, configs as C
>>> from simulation_smile.optics import PupilPSF
>>> pipe = SimulationPipeline(grid=C.Grid(pixel_size=108e-9, size=64))
>>> pipe = pipe.with_optics(PupilPSF(C.PSFConfig(na=1.49)))
"""
from __future__ import annotations

from . import configs, presets
from .core import set_backend
from .pipeline import Scene, SimulationPipeline

__version__ = "0.1.0"

__all__ = [
    "SimulationPipeline",
    "Scene",
    "configs",
    "presets",
    "set_backend",
    "__version__",
]
