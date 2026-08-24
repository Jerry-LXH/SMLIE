# 照明模型: 基类 + 具体激发分布
from .base import Illumination
from .uniform import Uniform
from .gaussian2d import Gaussian2D

__all__ = ["Illumination", "Uniform", "Gaussian2D"]
