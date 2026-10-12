# 探测器模型: 基类 + 具体相机
from .base import Detector
from .simple_cmos import SimpleCMOS
from .emccd import SimpleEMCCD

__all__ = ["Detector", "SimpleCMOS", "SimpleEMCCD"]
