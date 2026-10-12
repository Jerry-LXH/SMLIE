"""two_emitter.psf —— 直接拿**仿真器真正使用**的 PSF 核当模型模板。

背景
----
`level1` 默认用解析圆孔径 Airy 当模型 PSF。它与仿真器（``CircularAperturePSF``，
光瞳 FFT + 过采样核 + 双线性重采样）之间存在 **0.53% RMS of peak** 的形状差，
表现为"模型 PSF 偏宽 ~2–3%"（峰值占比 0.1678 vs 0.1784）。

已在 ``test_two_emitter.ipynb`` §7.6/§8.4 中定量证实：这个失配是残余间距偏差地板
（−1.4…−3 nm）的主因 —— 判决性实验是"用模型自己造数据"时偏差降到 ≤0.12 nm。

本模块的做法
------------
不重新推导 PSF，而是**复用仿真器的过采样核与同一套双线性插值**，
渲染 ROI 内的像素化模板。因此模板与仿真器实际投到像素上的权重**逐点一致**
（差别只来自 float 精度），失配被彻底消除。

用法
----
>>> from two_emitter import psf as te_psf
>>> psf_obj, K = te_psf.get_kernel(meta["psf"])          # 核只算一次(FFT)
>>> P = te_psf.roi_template(psf_obj, K, r_nm, roi_origin, h, w, pixel_size)
通常不必直接调用 —— ``build_pair_hmm(..., psf_model='template')`` 会自动接好。
"""
from __future__ import annotations

import os
import sys

import numpy as np

#: 仿真包根目录（与 dataset.DATASET_ROOT 同样是"项目内固定路径 + 环境变量可覆盖"的风格）
SIM_ROOT = os.environ.get(
    "SIMULATION_SMILE_ROOT",
    r"E:\CodingProj\SMLIE\Projects\Smile_Codes\Simulation_Smile")

#: 缓存：metadata 描述 -> (CircularAperturePSF, 过采样核)
_PSF_CACHE: dict = {}
_KERNEL_CACHE: dict = {}


def _ensure_sim_path() -> None:
    """保证 simulation_smile 可导入（分析 notebook 未加 SIM_ROOT 时也能用）。"""
    if SIM_ROOT not in sys.path:
        sys.path.insert(0, SIM_ROOT)


def _meta_key(meta_psf: dict, z_nm: float) -> tuple:
    return (float(meta_psf.get("na", 0.0)),
            float(meta_psf.get("wavelength", 0.0)),
            int(meta_psf.get("pupil_size", 128)),
            float(meta_psf.get("oversample", 4.0)),
            float(meta_psf.get("focus", 0.0)),
            repr(meta_psf.get("zernike")),
            float(z_nm))


def get_psf(meta_psf: dict):
    """按数据集 metadata 复现仿真用的 ``CircularAperturePSF``（带缓存）。"""
    key = _meta_key(meta_psf, 0.0)
    if key not in _PSF_CACHE:
        _ensure_sim_path()
        from simulation_smile.optics import CircularAperturePSF
        _PSF_CACHE[key] = CircularAperturePSF(
            na=float(meta_psf["na"]),
            wavelength=float(meta_psf["wavelength"]),
            zernike=meta_psf.get("zernike"),
            pupil_size=int(meta_psf.get("pupil_size", 128)),
            oversample=float(meta_psf.get("oversample", 4.0)),
            focus=float(meta_psf.get("focus", 0.0)))
    return _PSF_CACHE[key]


def get_kernel(meta_psf: dict, z_nm: float = 0.0):
    """返回 ``(psf_obj, K)``：过采样强度核 ``K``（间距 ``psf.kernel_spacing`` nm, 求和=1）。

    一次 FFT 的体积计算，之后复用；``K`` 形状 ``(n_pad, n_pad)``。
    """
    key = _meta_key(meta_psf, z_nm)
    if key not in _KERNEL_CACHE:
        psf = get_psf(meta_psf)
        _KERNEL_CACHE[key] = (psf, psf.psf_at(float(z_nm)))
    return _KERNEL_CACHE[key]


def roi_template(psf, K: np.ndarray, r_nm, roi_origin, h: int, w: int,
                 pixel_size_nm: float) -> np.ndarray:
    """ROI 内的像素化 PSF 模板 ``(h, w)``，与仿真器 ``render`` 逐点一致。

    复刻 ``CircularAperturePSF._build_render_kernels`` 的算法:
    在像素中心处对过采样核做双线性插值（同时完成到 ``pixel_size`` 的重采样与
    亚像素平移），再乘通量因子 ``(px/h)²``。

    约定与 ``level1._airy_pixels`` 相同：总能量归一（核求和=1），截到 ROI，
    因此返回值的和 = **ROI 内圈入比例**（ROI 外的光子不进似然）。

    参数
    ----
    r_nm        : 发射体位置 (x, y)（全视场, nm）
    roi_origin  : ROI 在全视场中的起点 (r0, c0)（像素索引, 与 dataset.roi 同约定）
    """
    n_pad = K.shape[0]
    hsp = float(psf.kernel_spacing)
    kc = (np.arange(n_pad) - n_pad // 2) * hsp              # 核样点偏移 (nm)
    cols = (np.arange(int(roi_origin[1]), int(roi_origin[1]) + w) + 0.5) * pixel_size_nm
    rows = (np.arange(int(roi_origin[0]), int(roi_origin[0]) + h) + 0.5) * pixel_size_nm
    dy = rows - float(r_nm[1])
    dx = cols - float(r_nm[0])
    k = psf._bilinear(K, kc, dy, dx) * (pixel_size_nm / hsp) ** 2
    return k


def roi_template_from_meta(meta_psf: dict, r_nm, roi_origin, h: int, w: int,
                           pixel_size_nm: float, z_nm: float = 0.0) -> np.ndarray:
    """``get_kernel`` + ``roi_template`` 的便捷封装（适合单次调用/自检）。"""
    psf, K = get_kernel(meta_psf, z_nm=z_nm)
    return roi_template(psf, K, r_nm, roi_origin, h, w, pixel_size_nm)
