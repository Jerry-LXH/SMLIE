# smlm/localize/loc/astigmatic.py
"""3D 像散 (柱面镜) 定位: 各向异性椭圆高斯模型。

模型 theta = [x0, y0, logN, logb, logsx, logsy]:
    model = b + N * Gx(x - x0, σx) * Gy(y - y0, σy)

像散方案的 z 信息编码在 (σx, σy) 的椭圆度里:
  - 本类负责拟合出 (x0, y0, N, b, σx, σy) 及横向 Fisher 精度;
  - 若传入 z_calibration (σx, σy -> z, z_lp), 输出 Localizations.z / z_lp;
  - 否则 z / z_lp 为 None, σx / σy 保存在 extra['sigma_x'] / extra['sigma_y']。

与 2D 方案共享: ROI 提取、多进程 worker 池、Fisher 矩阵机制。
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
# 椭圆高斯模型 + 解析导数
# ----------------------------------------------------------------------
def _model_and_deriv(theta, xx, yy, use_pixel):
    """返回 (model, dmodel); dmodel shape=(6,)+roi.shape。"""
    x0, y0 = theta[0], theta[1]
    N = np.exp(theta[2]); b = np.exp(theta[3])
    sx = np.exp(theta[4]); sy = np.exp(theta[5])
    A = (xx - x0) ** 2
    B = (yy - y0) ** 2

    if use_pixel:
        Ex = 0.5 * (erf((xx - x0 + 0.5) / (SQRT2 * sx)) - erf((xx - x0 - 0.5) / (SQRT2 * sx)))
        Ey = 0.5 * (erf((yy - y0 + 0.5) / (SQRT2 * sy)) - erf((yy - y0 - 0.5) / (SQRT2 * sy)))
        E = Ex * Ey
        model = b + N * E

        gx_p = np.exp(-((xx - x0 + 0.5) ** 2) / (2 * sx ** 2))
        gx_m = np.exp(-((xx - x0 - 0.5) ** 2) / (2 * sx ** 2))
        gy_p = np.exp(-((yy - y0 + 0.5) ** 2) / (2 * sy ** 2))
        gy_m = np.exp(-((yy - y0 - 0.5) ** 2) / (2 * sy ** 2))
        dEx_dx = (gx_m - gx_p) / (np.sqrt(2 * np.pi) * sx)
        dEy_dy = (gy_m - gy_p) / (np.sqrt(2 * np.pi) * sy)

        d_dx = N * dEx_dx * Ey
        d_dy = N * Ex * dEy_dy
        d_dlogN = N * E
        # σ 用中心高斯近似 (与 2D 方案一致)
        NG = N * np.exp(-(A / (2 * sx ** 2) + B / (2 * sy ** 2))) / (2 * np.pi * sx * sy)
        d_dlogsx = NG * (A / sx ** 2 - 1.0)
        d_dlogsy = NG * (B / sy ** 2 - 1.0)
    else:
        G = np.exp(-(A / (2 * sx ** 2) + B / (2 * sy ** 2))) / (2 * np.pi * sx * sy)
        NG = N * G
        model = b + NG
        d_dx = NG * (xx - x0) / sx ** 2
        d_dy = NG * (yy - y0) / sy ** 2
        d_dlogN = NG
        d_dlogsx = NG * (A / sx ** 2 - 1.0)
        d_dlogsy = NG * (B / sy ** 2 - 1.0)

    d_dlogb = b * np.ones_like(model)
    dmodel = np.stack([d_dx, d_dy, d_dlogN, d_dlogb, d_dlogsx, d_dlogsy], axis=0)
    return model, dmodel


def _nll_and_grad(theta, roi, xx, yy, F_EM, use_pixel):
    """负对数似然及其解析梯度 (与 2D 方案同一噪声机制)。"""
    model, dmodel = _model_and_deriv(theta, xx, yy, use_pixel)
    v = F_EM ** 2 * (model + 1e-12)
    diff = roi - model
    L = 0.5 * np.sum(diff ** 2 / v + np.log(v))
    dL_dm = 0.5 * (-2 * diff / v - diff ** 2 * F_EM ** 2 / v ** 2 + F_EM ** 2 / v)
    grad = np.array([np.sum(dL_dm * dmodel[j]) for j in range(6)])
    return L, grad


def _nll(theta, roi, xx, yy, F_EM, use_pixel):
    model, _ = _model_and_deriv(theta, xx, yy, use_pixel)
    v = F_EM ** 2 * (model + 1e-12)
    return 0.5 * np.sum((roi - model) ** 2 / v + np.log(v))


# ----------------------------------------------------------------------
# Fisher 信息 (像素积分模型; 复用共享机制)
# ----------------------------------------------------------------------
def fisher_info_astigmatic(theta, xx, yy, F_EM):
    """6x6 Fisher 信息矩阵: [x0, y0, logN, logb, logsx, logsy]。"""
    model, dmodel = _model_and_deriv(theta, xx, yy, use_pixel=True)
    return fisher_from_jacobian(dmodel, model, F_EM)


# ----------------------------------------------------------------------
# 单 ROI 定位 —— 模块级函数, 可被多进程 pickle
# 返回 [y0, x0, N, b, sx, sy, sigma_lateral] (坐标相对 ROI 中心)
# ----------------------------------------------------------------------
def localize_single(roi, F_EM=SQRT2, use_pixel=False, use_gradient=True):
    ny = (roi.shape[0] - 1) / 2
    nx = (roi.shape[1] - 1) / 2
    yy, xx = np.meshgrid(np.arange(-ny, ny + 1), np.arange(-nx, nx + 1), indexing='ij')

    b0 = max(np.median(roi), 1.0)
    N0 = max(np.sum(np.clip(roi - b0, 0, None)), 10.0)
    theta0 = np.array([0.0, 0.0, np.log(N0), np.log(b0), np.log(1.2), np.log(1.2)])

    bounds = [(-ny, ny), (-nx, nx),      # 位置限制在 ROI 内
              (np.log(1.0), np.log(1e7)),# 光子数 1 ~ 1e7
              (np.log(1e-3), np.log(1e5)),# 背景
              (np.log(0.3), np.log(10.0)),# σx 0.3~10 px
              (np.log(0.3), np.log(10.0))]# σy 0.3~10 px

    if use_gradient:
        res = minimize(_nll_and_grad, theta0, args=(roi, xx, yy, F_EM, use_pixel),
                       method='L-BFGS-B', jac=True, bounds=bounds)
    else:
        res = minimize(_nll, theta0, args=(roi, xx, yy, F_EM, use_pixel),
                       method='L-BFGS-B', bounds=bounds)
    th = res.x

    x0, y0 = th[0], th[1]
    Np, b = np.exp(th[2]), np.exp(th[3])
    sx, sy = np.exp(th[4]), np.exp(th[5])

    try:
        C = np.linalg.inv(fisher_info_astigmatic(th, xx, yy, F_EM))
        sigma_lateral = np.sqrt(0.5 * (C[0, 0] + C[1, 1]))
    except np.linalg.LinAlgError:
        sigma_lateral = np.nan

    return np.array([y0, x0, Np, b, sx, sy, sigma_lateral])


# ----------------------------------------------------------------------
# Localizer 接口
# ----------------------------------------------------------------------
class AstigmaticLocalizer(Localizer):
    """3D 像散 (柱面镜) MLE 定位, 各向同性高斯的各向异性推广。

    参数
    ----
    edge / F_EM / use_pixel / use_gradient / n_jobs: 同 MLELocalizer。
    z_calibration : 可选。恢复 z 的标定:
        - 可调用对象 calib(sx, sy) -> (z, z_lp), z 单位 nm;
        - 或 dict {'sigma_x': [...], 'sigma_y': [...], 'z': [...],
          'z_lp': [...]}, 按 (σx, σy) 最近邻查表。
      提供时 Localizations.z / z_lp 被填充; 否则为 None,
      σx / σy 保存在 extra['sigma_x'] / extra['sigma_y']。
    """

    def __init__(self, edge=4, F_EM=SQRT2, use_pixel=False,
                 use_gradient=True, n_jobs=1, verbose=True, z_calibration=None):
        self.edge = int(edge)
        self.F_EM = F_EM
        self.use_pixel = use_pixel
        self.use_gradient = use_gradient
        self.n_jobs = int(n_jobs)
        self.verbose = verbose
        self.z_calibration = z_calibration

    def localize(self, stack, cand):
        rois, meta = extract_rois(stack, cand, self.edge)
        if self.verbose:
            print(f"{len(rois)} / {len(cand)} candidates in-range, fitting (astigmatic)...")

        ests = run_pool(localize_single, rois,
                        (self.F_EM, self.use_pixel, self.use_gradient),
                        n_jobs=self.n_jobs, desc="Astigmatic fitting",
                        verbose=self.verbose)

        # 过滤 + 转绝对坐标
        out = []
        for (f, rr, cc), est in zip(meta, ests):
            if est is None:
                continue
            y_off, x_off, Np, b, sx, sy, sloc = est
            if not (np.isfinite(Np) and Np > 0 and
                    np.isfinite(sx) and sx > 0 and
                    np.isfinite(sy) and sy > 0 and
                    np.isfinite(sloc) and sloc > 0):
                continue
            out.append([f, cc + x_off, rr + y_off, Np, b, sx, sy, sloc])

        out = np.array(out) if out else np.empty((0, 8))
        if self.verbose:
            print(f"Astigmatic done, {out.shape[0]} localizations kept.")

        # 标定: (σx, σy) -> (z, z_lp)
        z, z_lp = None, None
        if out.shape[0] and self.z_calibration is not None:
            z = np.empty(out.shape[0])
            z_lp = np.empty(out.shape[0])
            for i in range(out.shape[0]):
                zi, lpi = self._calibrate(out[i, 5], out[i, 6])
                z[i], z_lp[i] = zi, lpi

        return Localizations(
            frame=out[:, 0].astype(int),
            x=out[:, 1], y=out[:, 2],
            photons=out[:, 3], bg=out[:, 4],
            sigma=None, lp=out[:, 7],
            z=z, z_lp=z_lp,
            extra={'sigma_x': out[:, 5], 'sigma_y': out[:, 6]},
        )

    # ------------------------------------------------------------------
    def _calibrate(self, sx, sy):
        cal = self.z_calibration
        if callable(cal):
            return cal(sx, sy)
        # dict 查表: 最近邻 (σx, σy)
        sx_arr = np.asarray(cal['sigma_x'], dtype=float)
        sy_arr = np.asarray(cal['sigma_y'], dtype=float)
        d2 = (sx_arr - sx) ** 2 + (sy_arr - sy) ** 2
        i = int(np.argmin(d2))
        return float(cal['z'][i]), float(cal['z_lp'][i])
