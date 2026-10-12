# smlm/imgprocess -- 图像处理流水线 (ImageStack 进 -> ImageStack 出)
"""`smlm.imgprocess` 取代旧的 `smlm/preprocess` + `smlm/sofi`。

按功能划分子模块:
  - base       : Transform / Pipeline 抽象
  - preprocess : 纯预处理 (Crop / SelectFrames)
  - temporal   : 时间维度处理 (光漂白校正 / 漂移校正)
  - interp     : 傅里叶插值 / 上采样
  - denoise    : 去噪 (高斯预滤波)
  - dehotpixel : sCMOS 热点像素去除
  - deconv     : 反卷积 (RL3D / SparseDeconv / OtfReweight)
  - core       : SOFI 累积量 + RIED 熵加权相关的纯数学
  - sofi       : SOFI/RIED 估计器 (Transform 子类)
  - viz        : 结果可视化

用法:
    from smlm.imgprocess import (Pipeline, Crop, SelectFrames, BleachCorrect,
        FourierUpsample, GaussianFilter, RL3D, RiedSofit, SparseDeconv)

    out = Pipeline(BleachCorrect("mono_exp"),
                   GaussianFilter(0.7),
                   FourierUpsample(2),
                   RL3D(psf, iterations=4),
                   RiedSofit(subfactor=0.5, entropy=True),
                   SparseDeconv(psf2, fidelity=50, sparsity=0.5, iterations=5),
                   )(stack)      # -> ImageStack (单帧)
"""
from .base import Transform, Pipeline
from .preprocess import Crop, SelectFrames
from .temporal import (BleachCorrect, DriftCorrector,
                       bleach_correct, apply_drift, TemporalMean)
from .interp import (FourierUpsample, fourier_zoom, fourier_upsample,
                     fourier_interp_2d, fourier_lateral)
from .denoise import GaussianFilter
from .dehotpixel import Dehotpixel, dehotpixel
from .deconv import (RL3D, SparseDeconv, OtfReweight,
                     rl_deconv, sparse_deconv, otf_reweight,
                     scalar_psf, otf_components, psf2otf,
                     forward_diff, back_diff)
from . import core
from .sofi import Sofit, CumulantSofit, RiedSofit
from . import viz

__all__ = [
    "Transform", "Pipeline",
    "Crop", "SelectFrames",
    "BleachCorrect", "DriftCorrector",
    "bleach_correct", "apply_drift",
    "FourierUpsample", "fourier_zoom", "fourier_upsample",
    "fourier_interp_2d", "fourier_lateral", "TemporalMean",
    "GaussianFilter",
    "Dehotpixel", "dehotpixel",
    "RL3D", "SparseDeconv", "OtfReweight",
    "rl_deconv", "sparse_deconv", "otf_reweight",
    "scalar_psf", "otf_components", "psf2otf",
    "forward_diff", "back_diff",
    "core",
    "Sofit", "CumulantSofit", "RiedSofit",
    "viz",
]
