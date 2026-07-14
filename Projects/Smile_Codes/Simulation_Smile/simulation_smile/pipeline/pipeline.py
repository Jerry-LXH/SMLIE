"""The composable simulation pipeline.

SKELETON -- the builder API is sketched; render bodies are stubbed.

``SimulationPipeline`` is an ordered list of stages plus the scene. The builder
methods (``with_optics`` etc.) are thin wrappers that append a stage in the
canonical order and return ``self`` for chaining. ``render_frame`` applies every
stage to a fresh :class:`SimContext`; ``render_movie`` advances the photophysics
across frames via :mod:`timeline`.

The order stages run is fixed and physically meaningful:

    photophysics -> illumination -> optics -> camera -> noise

but any stage can be omitted (e.g. drop photophysics for a static snapshot).
"""
from __future__ import annotations

from typing import Sequence

import numpy as np

from ..core.base import Grid, SimContext, Stage
from .scene import Scene


class SimulationPipeline:
    """Builder + executor for a microscope-imaging simulation."""

    def __init__(self, grid: Grid | None = None, scene: Scene | None = None):
        self.grid = grid if grid is not None else (scene.grid if scene else None)
        self.scene = scene
        # Named slots, applied in this canonical order during render.
        self.photophysics: Stage | None = None
        self.illumination: Stage | None = None
        self.optics: Stage | None = None
        self.camera: Stage | None = None
        self.noise: list[Stage] = []

    # --- builder API ------------------------------------------------------
    def with_scene(self, scene: Scene) -> "SimulationPipeline":
        self.scene = scene
        self.grid = scene.grid
        return self

    def with_photophysics(self, stage: Stage) -> "SimulationPipeline":
        self.photophysics = stage
        return self

    def with_illumination(self, stage: Stage) -> "SimulationPipeline":
        self.illumination = stage
        return self

    def with_optics(self, stage: Stage) -> "SimulationPipeline":
        self.optics = stage
        return self

    def with_camera(self, stage: Stage) -> "SimulationPipeline":
        self.camera = stage
        return self

    def with_noise(self, stages: Sequence[Stage]) -> "SimulationPipeline":
        self.noise = list(stages)
        return self

    # --- ordered stage list ----------------------------------------------
    def stages(self) -> list[Stage]:
        """Flatten the named slots into the canonical execution order."""
        ordered: list[Stage] = []
        for s in (self.photophysics, self.illumination, self.optics, self.camera):
            if s is not None:
                ordered.append(s)
        ordered.extend(self.noise)
        return ordered

    # --- execution --------------------------------------------------------
    def render_frame(self, frame_index: int = 0) -> np.ndarray:
        """Render a single frame: apply every stage to a fresh context.

        TODO: build a SimContext from the scene, run each stage's ``apply``,
        return ``ctx.image``.
        """
        raise NotImplementedError

    def render_movie(self, n_frames: int) -> np.ndarray:
        """Render ``n_frames`` advancing photophysics between frames.

        TODO: delegate to :func:`..timeline.run_movie`, which resets stateful
        stages once then steps the photophysics each frame. Returns a stack of
        shape (n_frames, H, W).
        """
        raise NotImplementedError
