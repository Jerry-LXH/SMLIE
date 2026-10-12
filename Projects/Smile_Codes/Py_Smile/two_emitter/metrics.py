"""two_emitter.metrics —— 位置 / 间距误差指标 (小工具集合)。"""
from __future__ import annotations

import numpy as np


def abs_error(a, b) -> np.ndarray:
    """逐元素绝对误差。"""
    return np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))


def bias(a, b) -> float:
    """平均偏差 mean(a − b); 忽略 nan。"""
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.nanmean(d))


def rmse(a, b) -> float:
    """均方根误差 sqrt(mean((a − b)²)); 忽略 nan。"""
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.nanmean(d ** 2)))


def std(a) -> float:
    """标准差 (忽略 nan)。"""
    return float(np.nanstd(np.asarray(a, dtype=float)))


def position_error_2d(xy_est, xy_true) -> np.ndarray:
    """二维点集逐点欧氏距离 (nm)。xy_est / xy_true 形状 (N, 2)。"""
    e = np.asarray(xy_est, float) - np.asarray(xy_true, float)
    return np.hypot(e[:, 0], e[:, 1])


def separation_error_nm(d_est, d_true) -> np.ndarray:
    """间距误差 (带符号)。"""
    return np.asarray(d_est, float) - np.asarray(d_true, float)


def summarize(err_nm) -> dict:
    """把一组误差汇总成 mean / std / rmse / median / p95 (nm)。"""
    e = np.asarray(err_nm, dtype=float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return {"n": 0, "mean_nm": np.nan, "std_nm": np.nan,
                "rmse_nm": np.nan, "median_nm": np.nan, "p95_nm": np.nan}
    return {
        "n": int(e.size),
        "mean_nm": float(e.mean()),
        "std_nm": float(e.std(ddof=1)) if e.size > 1 else 0.0,
        "rmse_nm": float(np.sqrt(np.mean(e ** 2))),
        "median_nm": float(np.median(e)),
        "p95_nm": float(np.percentile(e, 95)),
    }
