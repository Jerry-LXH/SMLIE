"""Photophysics tests (skeleton)."""
import numpy as np
import pytest

from simulation_smile import configs as C
from simulation_smile.fluorophore import BlinkingModel


@pytest.mark.xfail(reason="BlinkingModel not implemented yet", strict=False)
def test_bleaching_is_exponential():
    """Total ON population decays ~ exp(-k_bleach * t)."""
    rng = np.random.default_rng(0)
    cfg = C.FluorophoreConfig(k_on=0.0, k_off=0.0, k_bleach=1.0, initial_state="on")
    model = BlinkingModel(cfg)
    # TODO: init N emitters ON, step dt repeatedly, fit decay, compare to k_bleach.


@pytest.mark.xfail(reason="BlinkingModel not implemented yet", strict=False)
def test_on_duty_cycle():
    """Steady-state ON fraction ~ k_on / (k_on + k_off) without bleaching."""
    cfg = C.FluorophoreConfig(k_on=1.0, k_off=4.0, k_bleach=0.0)
    model = BlinkingModel(cfg)
    # TODO: long run, measure ON fraction ~ 0.2.
