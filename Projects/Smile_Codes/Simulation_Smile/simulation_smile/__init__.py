"""simulation_smile -- 可组合的显微镜成像模拟器。

结构对齐 py_smile/smlm:
  - datatypes : 核心数据类型与全局约定 (nm / s / 光子)
  - sources   : 发射体分布
  - illumination : 激发分布
  - photophysics : 光物理 (闪烁/漂白, 帧内 Gillespie)
  - optics    : PSF 模型 (局部窗口渲染 + 核缓存)
  - detector  : 相机模型
  - pipeline  : SimulationPipeline 编排整条管线

随机数约定: 一律使用外部传入的 np.random.Generator, 模块内部绝不 seed。
"""
from .datatypes import EmitterSet, PhotonStream, SimulationResult, FieldOfView
from .pipeline import SimulationPipeline
from .sources import (EmitterSource, Random2D, Random3D, SimpleRing2d,
                     SimpleCircles2d, SimpleLines, PairGrid)
from .illumination import Illumination, Uniform, Gaussian2D
from .photophysics import PhotophysicsModel, Static, Blinking
from .optics import (PSFModel, GaussianPSF, CircularAperturePSF,
                     CalibratedPSF)
from .detector import Detector, SimpleCMOS, SimpleEMCCD

__version__ = "0.1.0"

__all__ = [
    "EmitterSet", "PhotonStream", "SimulationResult", "FieldOfView",
    "SimulationPipeline",
    "EmitterSource", "Random2D", "Random3D", "SimpleRing2d", "SimpleCircles2d", "SimpleLines",
    "PairGrid",
    "Illumination", "Uniform", "Gaussian2D",
    "PhotophysicsModel", "Static", "Blinking",
    "PSFModel", "GaussianPSF", "CircularAperturePSF", "CalibratedPSF",
    "Detector", "SimpleCMOS", "SimpleEMCCD",
    "__version__",
]
