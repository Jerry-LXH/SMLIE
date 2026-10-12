# smlm/camera/emccd.py
from .base import Camera
import numpy as np

class EMCCD(Camera):
    def __init__(self, conv_factor=5.75, em_gain=None, offset=0.0):
        self.conv_factor = conv_factor      # e-/ADU (你的 5.75)
        self.em_gain = em_gain
        self.offset = offset

    def adu_to_photons(self, images: np.ndarray) -> np.ndarray:
        if self.em_gain is None:
            raise ValueError("EMCCD 需要 em_gain,请在构造时传入或从 raw.em_gain 取")
        return (images - self.offset) * self.conv_factor / self.em_gain