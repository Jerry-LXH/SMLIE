# smlm/localize/loc/gaussian2d.py
"""2D 各向同性高斯 MLE 定位 (对应 MATLAB locSingle.m / locMolecules.m)。

模型 theta = [x0, y0, logN, logb, logs]:
  - use_pixel=False: 连续高斯  b + N/(2πσ²) exp(-r²/2σ²)
  - use_pixel=True : 像素积分高斯 (逐像素 erf 积分)

产物 Localizations 的 z / z_lp 保持 None (2D 定位)。
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from scipy.special import erf

from ..base import Localizer
from ...datatypes import Localizations
from .common import extract_rois, run_pool
from .fisher import fisher_from_jacobian

SQRT2 = np.sqrt(2.0)


# ----------------------------------------------------------------------
# 模型 + 解析导数 (供 NLL 梯度用)
# ----------------------------------------------------------------------
def _model_and_deriv(theta, xx, yy, use_pixel):
    """返回 (model, dmodel) ; dmodel shape=(5,)+roi.shape, 对 theta 各分量求导。"""
    x0, y0 = theta[0], theta[1]
    N = np.exp(theta[2]); b = np.exp(theta[3]); s = np.exp(theta[4])
    r2 = (xx - x0) ** 2 + (yy - y0) ** 2

    if use_pixel:
        Ex = 0.5 * (erf((xx - x0 + 0.5) / (SQRT2 * s)) - erf((xx - x0 - 0.5) / (SQRT2 * s)))
        Ey = 0.5 * (erf((yy - y0 + 0.5) / (SQRT2 * s)) - erf((yy - y0 - 0.5) / (SQRT2 * s)))
        E = Ex * Ey
        model = b + N * E

        gx_p = np.exp(-((xx - x0 + 0.5) ** 2) / (2 * s ** 2))
        gx_m = np.exp(-((xx - x0 - 0.5) ** 2) / (2 * s ** 2))
        gy_p = np.exp(-((yy - y0 + 0.5) ** 2) / (2 * s ** 2))
        gy_m = np.exp(-((yy - y0 - 0.5) ** 2) / (2 * s ** 2))
        dEx_dx = (gx_m - gx_p) / (np.sqrt(2 * np.pi) * s)
        dEy_dy = (gy_m - gy_p) / (np.sqrt(2 * np.pi) * s)

        d_dx = N * dEx_dx * Ey
        d_dy = N * Ex * dEy_dy
        d_dlogN = N * E
        # σ 用中心高斯近似 (与你 .m 一致)
        NG = N * np.exp(-r2 / (2 * s ** 2)) / (2 * np.pi * s ** 2)
        d_dlogs = NG * (r2 / s ** 2 - 2.0)
    else:
        G = np.exp(-r2 / (2 * s ** 2)) / (2 * np.pi * s ** 2)
        NG = N * G
        model = b + NG
        d_dx = NG * (xx - x0) / s ** 2
        d_dy = NG * (yy - y0) / s ** 2
        d_dlogN = NG
        d_dlogs = NG * (r2 / s ** 2 - 2.0)

    d_dlogb = b * np.ones_like(model)
    dmodel = np.stack([d_dx, d_dy, d_dlogN, d_dlogb, d_dlogs], axis=0)
    return model, dmodel


def _nll_and_grad(theta, roi, xx, yy, F_EM, use_pixel):
    """负对数似然及其解析梯度 (对应 poissonApproxNLL_varSigma.m)。"""
    model, dmodel = _model_and_deriv(theta, xx, yy, use_pixel)
    v = F_EM ** 2 * (model + 1e-12)                       # 有效方差
    diff = roi - model
    L = 0.5 * np.sum(diff ** 2 / v + np.log(v))
    # dL/dmodel : 注意 v 依赖 model, 故有三项
    dL_dm = 0.5 * (-2 * diff / v - diff ** 2 * F_EM ** 2 / v ** 2 + F_EM ** 2 / v)
    grad = np.array([np.sum(dL_dm * dmodel[j]) for j in range(5)])
    return L, grad


def _nll(theta, roi, xx, yy, F_EM, use_pixel):
    model, _ = _model_and_deriv(theta, xx, yy, use_pixel)
    v = F_EM ** 2 * (model + 1e-12)
    return 0.5 * np.sum((roi - model) ** 2 / v + np.log(v))


# ----------------------------------------------------------------------
# Fisher 信息 (像素积分模型; 噪声机制复用共享的 fisher_from_jacobian)
# ----------------------------------------------------------------------
def fisher_info_2d(theta, xx, yy, F_EM):
    """5x5 Fisher 信息矩阵 (对应 fisherInfo_approx_varSigma.m)。"""
    model, dmodel = _model_and_deriv(theta, xx, yy, use_pixel=True)
    return fisher_from_jacobian(dmodel, model, F_EM)


# ----------------------------------------------------------------------
# 单 ROI 定位 —— 模块级函数, 可被多进程 pickle
# 返回 [y0, x0, N, b, sigma, sigma_loc] (坐标相对 ROI 中心)
# ----------------------------------------------------------------------
def localize_single(roi, F_EM=SQRT2, use_pixel=False, use_gradient=True):
    ny = (roi.shape[0] - 1) / 2
    nx = (roi.shape[1] - 1) / 2
    yy, xx = np.meshgrid(np.arange(-ny, ny + 1), np.arange(-nx, nx + 1), indexing='ij')

    b0 = max(np.median(roi), 1.0)
    N0 = max(np.sum(np.clip(roi - b0, 0, None)), 10.0)
    theta0 = np.array([0.0, 0.0, np.log(N0), np.log(b0), np.log(1.2)])

    bounds = [(-ny, ny), (-nx, nx),      # 位置限制在 ROI 内
            (np.log(1.0), np.log(1e7)), # 光子数 1 ~ 1e7
            (np.log(1e-3), np.log(1e5)),# 背景
            (np.log(0.3), np.log(10.0))]# sigma 0.3~10 px

    if use_gradient:
        res = minimize(_nll_and_grad, theta0, args=(roi, xx, yy, F_EM, use_pixel),
                    method='L-BFGS-B', jac=True, bounds=bounds)
    else:
        res = minimize(_nll, theta0, args=(roi, xx, yy, F_EM, use_pixel),
                    method='L-BFGS-B', bounds=bounds)
    th = res.x

    x0, y0 = th[0], th[1]
    Np, b, sigma = np.exp(th[2]), np.exp(th[3]), np.exp(th[4])

    try:
        C = np.linalg.inv(fisher_info_2d(th, xx, yy, F_EM))
        sigma_loc = np.sqrt(0.5 * (C[0, 0] + C[1, 1]))
    except np.linalg.LinAlgError:
        sigma_loc = np.nan

    return np.array([y0, x0, Np, b, sigma, sigma_loc])


# ----------------------------------------------------------------------
# Localizer 接口 (对应 locMolecules.m)
# ----------------------------------------------------------------------
class MLELocalizer(Localizer):
    """2D 各向同性高斯 MLE 亚像素定位。

    stack.data 约定 (T, H, W); cand 需有 .frame/.x(col)/.y(row)。

    edge        : ROI 半宽, ROI 尺寸 = (2*edge+1)。
    F_EM        : EMCCD 过噪因子 (非 EMCCD 设 1)。
    use_pixel   : True=像素积分高斯, False=连续高斯。
    use_gradient: True 用解析梯度 (快很多)。
    n_jobs      : >1 启用多进程。

    产物为 2D 定位: z / z_lp 均为 None。
    """

    def __init__(self, edge=4, F_EM=SQRT2, use_pixel=False,
                 use_gradient=True, n_jobs=1, verbose=True):
        self.edge = int(edge)
        self.F_EM = F_EM
        self.use_pixel = use_pixel
        self.use_gradient = use_gradient
        self.n_jobs = int(n_jobs)
        self.verbose = verbose

    def localize(self, stack, cand):
        rois, meta = extract_rois(stack, cand, self.edge)
        if self.verbose:
            print(f"{len(rois)} / {len(cand)} candidates in-range, fitting...")

        ests = run_pool(localize_single, rois,
                        (self.F_EM, self.use_pixel, self.use_gradient),
                        n_jobs=self.n_jobs, desc="MLE fitting",
                        verbose=self.verbose)

        # 过滤 + 转绝对坐标 (0-indexed 像素中心, 与 detect/viz 同一约定)
        out = []
        for (f, rr, cc), est in zip(meta, ests):
            if est is None:
                continue
            y_off, x_off, Np, b, sigma, sloc = est
            if not (np.isfinite(Np) and Np > 0 and
                    np.isfinite(sigma) and sigma > 0 and
                    np.isfinite(sloc) and sloc > 0):
                continue
            out.append([f, cc + x_off, rr + y_off, Np, b, sigma, sloc])

        out = np.array(out) if out else np.empty((0, 7))
        if self.verbose:
            print(f"MLE done, {out.shape[0]} localizations kept.")

        return Localizations(
            frame=out[:, 0].astype(int),
            x=out[:, 1], y=out[:, 2],
            photons=out[:, 3], bg=out[:, 4],
            sigma=out[:, 5], lp=out[:, 6],
            z=None, z_lp=None,
        )
