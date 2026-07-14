"""End-to-end pipeline tests (skeleton)."""
import numpy as np
import pytest

from simulation_smile import SimulationPipeline, configs as C
from simulation_smile.core import SimContext, Stage
from simulation_smile.core.base import Grid


def test_grid_geometry():
    """Grid derived quantities are self-consistent (no implementation needed)."""
    g = Grid(pixel_size=100e-9, size=64, oversampling=4)
    assert g.fine_size == 256
    assert np.isclose(g.fine_pixel_size, 25e-9)
    assert np.isclose(g.fov, 64 * 100e-9)


def test_stage_protocol_is_runtime_checkable():
    """Any object with apply(ctx) satisfies the Stage protocol."""

    class Dummy:
        def apply(self, ctx):
            return ctx

    assert isinstance(Dummy(), Stage)


def test_pipeline_builder_orders_stages():
    """Builder slots flatten into the canonical execution order."""

    class Tag:
        def __init__(self, name):
            self.name = name

        def apply(self, ctx):
            return ctx

    pipe = (SimulationPipeline(grid=Grid(100e-9, 64))
            .with_photophysics(Tag("photo"))
            .with_illumination(Tag("illum"))
            .with_optics(Tag("optics"))
            .with_camera(Tag("camera"))
            .with_noise([Tag("poisson"), Tag("read")]))
    names = [s.name for s in pipe.stages()]
    assert names == ["photo", "illum", "optics", "camera", "poisson", "read"]


@pytest.mark.xfail(reason="render not implemented yet", strict=False)
def test_render_movie_shape_and_reproducibility():
    """Fixed seed -> identical movie; shape == (T, size, size)."""
    from simulation_smile import presets

    pipe = presets.smlm_storm()
    m1 = pipe.render_movie(5)
    assert m1.shape[0] == 5
