# smlm/camera/scmos.py
from .base import Camera
class sCMOS(Camera):
    """sCMOS 需要逐像素校准图,不是单一标量。"""
    def __init__(self, gain_map, offset_map, var_map=None):
        self.gain_map = gain_map            # (H, W) e-/ADU
        self.offset_map = offset_map        # (H, W) 暗电平
        self.var_map = var_map              # (H, W) 读噪声方差, 定位时用

    def adu_to_photons(self, images):
        # images 约定 (T, H, W); 校准图 (H, W) 加前导帧轴广播到每一帧
        return (images - self.offset_map[None, ...]) / self.gain_map[None, ...]