"""two_emitter.level0 —— level0(已知态) 的理论极限与 oracle 估计器。

level0 问题
-----------
**假设每帧哪些分子在发光都已经知道**, 在给定的单帧光子数 / PSF / 像素尺寸 / 背景下,
两个静态发射体最多能靠多近还能被分开?

思路
----
1. 单帧单分子定位精度用 Thompson–Larson–Webb 的 CRLB 公式给出::

       σ_μ² = (s² + a²/12) / N  +  8π·s⁴·b² / (a²·N²)

   其中 ``s``=PSF 等效高斯标准差, ``a``=像素尺寸, ``N``=信号光子数,
   ``b``=每像素背景光子数(单位与 s、a 一致, 即 nm)。

2. 已知态后, 可以只用"只有它亮"的帧 (isolated frames) 定位它 —— 相邻帧独立,
   定位精度提升 ``√M`` 倍 (M = isolated 帧数)。

3. 两个发射体各自独立定位, 间距精度 ``σ_d = √(σ_A² + σ_B²)``。

4. 定义"可分辨"判据 ``d_min = k · σ_d`` (默认 k=2, 即 2σ)。

结论预告: 因为 isolated 帧数随电影长度线性增长, level0 的最小可分辨间距
可以做到 **远低于衍射极限**(~ 亚 nm 量级), 即"信息本来就在闪烁里"。

除理论外, 本模块还提供 **oracle 估计器**: 直接使用真值 ``states`` 挑 isolated 帧,
求和后取强度质心, 作为 level0 的实测量, 用来和 CRLB 对照。
"""
from __future__ import annotations

import numpy as np

#: EMCCD 过噪因子 (电子倍增的多余噪声)。
SQRT2 = float(np.sqrt(2.0))

# --------------------------------------------------------------------------- #
# PSF 等效宽度
# --------------------------------------------------------------------------- #
def airy_sigma_nm(na: float, wavelength_nm: float) -> float:
    """圆孔径 Airy 斑的"高斯等效"标准差 (nm)。

    经验换算: 半高全宽 FWHM ≈ 0.51·λ/NA, 高斯里 σ = FWHM/2.355,
    即 σ ≈ 0.2166·λ/NA。
    """
    return 0.2166 * float(wavelength_nm) / float(na)


def kernel_sigma_nm(kernel: np.ndarray, spacing_nm: float) -> float:
    """从 PSF 核的二阶矩估计等效 σ (nm)。含重尾, 会比 FWHM 等效值偏大。"""
    K = np.asarray(kernel, dtype=float)
    n = K.shape[0]
    ax = (np.arange(n) - n // 2) * float(spacing_nm)
    tot = K.sum()
    if tot <= 0:
        return np.nan
    cx = float((K.sum(axis=0) * ax).sum() / tot)
    cy = float((K.sum(axis=1) * ax).sum() / tot)
    vx = float((K.sum(axis=0) * (ax - cx) ** 2).sum() / tot)
    vy = float((K.sum(axis=1) * (ax - cy) ** 2).sum() / tot)
    return float(np.sqrt(0.5 * (vx + vy)))


# --------------------------------------------------------------------------- #
# CRLB
# --------------------------------------------------------------------------- #
def thompson_crlb_nm(sigma_psf_nm: float, pixel_size_nm: float, n_photons: float,
                     bg_per_px: float = 0.0,
                     excess_noise_factor: float = 1.0) -> float:
    """单帧单分子定位精度 σ_μ (nm)。Thompson–Larson–Webb (2002) 公式。

    excess_noise_factor: EMCCD 的电子倍增引入多余噪声 (F=√2),
        以 σ 上的乘性因子近似 (等价于信号光子数折半)。
    """
    s = float(sigma_psf_nm)
    a = float(pixel_size_nm)
    N = float(n_photons)
    b = float(bg_per_px)
    if N <= 0 or not np.isfinite(N):
        return np.inf
    var = (s ** 2 + a ** 2 / 12.0) / N + (8.0 * np.pi * s ** 4 * b ** 2) / (a ** 2 * N ** 2)
    return float(np.sqrt(var) * float(excess_noise_factor))


def effective_sigma_nm(sigma_single_nm: float, n_effective_frames: float) -> float:
    """把单帧精度按独立帧数 M 平均: σ_single / √M。"""
    m = np.maximum(np.asarray(n_effective_frames, dtype=float), 1e-12)
    return float(np.asarray(sigma_single_nm, dtype=float) / np.sqrt(m))


def sigma_separation_nm(sigma_a_nm: float, sigma_b_nm: float) -> float:
    """两发射体间距的精度 σ_d = √(σ_A² + σ_B²)。"""
    return float(np.hypot(float(sigma_a_nm), float(sigma_b_nm)))


def min_resolvable_separation_nm(sigma_d_nm: float, criterion: float = 2.0) -> float:
    """可分辨判据: d_min = criterion · σ_d (默认 2σ)。"""
    return float(criterion) * float(sigma_d_nm)


def level0_budget(photons_per_on_frame: float, sigma_psf_nm: float,
                  pixel_size_nm: float, bg_per_px: float, n_frames: int,
                  p_on: float, excess_noise_factor: float = 1.0,
                  criterion: float = 2.0) -> dict:
    """由光物理 / 光学参数一次性算出 level0 的定位精度与最小可分辨间距。

    p_on: 单个发射体的亮态占空比。一对独立发射体中,
          "只有一个亮" 的帧占比 = 2·p_on·(1−p_on)。
    """
    sigma1 = thompson_crlb_nm(sigma_psf_nm, pixel_size_nm, photons_per_on_frame,
                              bg_per_px, excess_noise_factor)
    n_iso = float(n_frames) * 2.0 * float(p_on) * (1.0 - float(p_on))
    sigma_eff = effective_sigma_nm(sigma1, n_iso)
    sigma_d = sigma_separation_nm(sigma_eff, sigma_eff)
    return {
        "sigma_single_frame_nm": sigma1,
        "n_isolated_frames": n_iso,
        "n_isolated_per_emitter": n_iso / 2.0,
        "sigma_emitter_nm": sigma_eff,
        "sigma_separation_nm": sigma_d,
        "d_min_nm": min_resolvable_separation_nm(sigma_d, criterion),
    }


# --------------------------------------------------------------------------- #
# oracle (使用真值 states)
# --------------------------------------------------------------------------- #
def _centroid_result(sub_img: np.ndarray, roi_slices, pixel_size_nm: float,
                     bg_percentile: float = 20.0, n_frames: int = 0):
    """强度质心 (兜底方法)。

    注意坐标约定: 像素**索引** i 的中心在 ``(i + 0.5)·px``,
    所以质心索引要 +0.5 才是物理坐标。
    """
    sub = np.asarray(sub_img, dtype=float)
    if sub.size == 0 or sub.sum() <= 0:
        return None
    bg = np.percentile(sub, bg_percentile) if bg_percentile > 0 else 0.0
    w = np.clip(sub - bg, 0.0, None)
    tot = w.sum()
    if tot <= 0:
        return None
    ys, xs = roi_slices
    yy, xx = np.mgrid[ys.start:ys.stop, xs.start:xs.stop]
    cy = float((w * yy).sum() / tot) + 0.5
    cx = float((w * xx).sum() / tot) + 0.5
    return {"x_nm": cx * pixel_size_nm, "y_nm": cy * pixel_size_nm,
            "photons": float(tot), "bg": float(bg), "sigma_px": np.nan,
            "lp_nm": np.nan, "method": "centroid",
            "n_isolated_frames": int(n_frames)}


def mle_localize_image(img, pixel_size_nm: float = 1.0, half_px: int | None = None,
                       F_EM: float = SQRT2, use_pixel: bool = True):
    """对单张 2D 图像做像素积分高斯 MLE 定位 —— 复用 ``smlm.localize.MLELocalizer``。

    输入图像尺寸须为 ``(2·half_px+1)²``; 像素 ``[i, j]`` 的中心位于
    ``((j+0.5)·px, (i+0.5)·px)``。返回 dict(含 ``x_nm/y_nm/photons/lp_nm``)
    或 ``None``(拟合失败)。
    """
    from smlm.datatypes import ImageStack, Candidates
    from smlm.localize import MLELocalizer

    img = np.asarray(img, dtype=float)
    h, w = img.shape
    if half_px is None:
        half_px = min((h - 1) // 2, (w - 1) // 2)
    half_px = int(half_px)
    if (h, w) != (2 * half_px + 1, 2 * half_px + 1):
        raise ValueError(f"img 尺寸应为 ({2 * half_px + 1}, {2 * half_px + 1}), "
                         f"得到 {img.shape}")

    stack = ImageStack(data=img[None].astype(np.float32),
                       pixel_size=float(pixel_size_nm))
    cand = Candidates(frame=np.array([0]),
                      y=np.array([half_px]), x=np.array([half_px]))
    loc = MLELocalizer(edge=half_px, F_EM=float(F_EM), use_pixel=bool(use_pixel),
                       use_gradient=True, verbose=False).localize(stack, cand)
    if len(loc) == 0:
        return None
    x_px, y_px = float(loc.x[0]), float(loc.y[0])
    return {
        "x_px": x_px, "y_px": y_px,
        "x_nm": (x_px + 0.5) * float(pixel_size_nm),
        "y_nm": (y_px + 0.5) * float(pixel_size_nm),
        "photons": float(loc.photons[0]),
        "bg": float(loc.bg[0]),
        "sigma_px": float(loc.sigma[0]),
        "lp_nm": float(loc.lp[0]) * float(pixel_size_nm),
        "method": "mle",
    }


def isolated_sum_image(ds, pair_idx: int, emitter_k: int, half_px: int = 6,
                       mode: str = "counts", min_photons: int = 1):
    """把"只有该发射体发光"的帧在 ROI 内求和 (独立泊松之和仍是泊松)。

    ★ 默认用 ``counts`` 判断(而非帧末 ``states``): 帧内 on→dark 跳变会让帧末状态
    变成 dark 但该帧仍有光子, 用 states 会把搭档的光子混进来, 造成明显偏置。

    返回 ``(summed, roi_slices, n_isolated)``; 无 isolated 帧时 ``summed=None``。
    """
    m = ds.isolated_masks(pair_idx, mode=mode, min_photons=min_photons)[emitter_k]
    roi = ds.roi(pair_idx, half_px=half_px)
    if not np.any(m):
        return None, roi, 0
    summed = ds.photons()[m][:, roi[0], roi[1]].sum(axis=0)
    return summed, roi, int(m.sum())


def oracle_positions_nm(ds, pair_idx: int, half_px: int = 6, mode: str = "counts",
                        min_photons: int = 1, F_EM: float = SQRT2,
                        method: str = "mle") -> list:
    """level0 oracle: 只用"只有它发光"的帧**求和**后定位 (默认像素积分高斯 MLE)。

    返回 ``list[dict | None]``, 与 ``ds.pair_emitters(pair_idx)`` 同序;
    每项含 ``x_nm / y_nm / photons / bg / sigma_px / lp_nm / n_isolated_frames``,
    坐标已换算到全视场 nm 坐标系。
    """
    out = []
    for k in range(len(ds.pair_emitters(pair_idx))):
        summed, roi, n_iso = isolated_sum_image(ds, pair_idx, k, half_px=half_px,
                                                mode=mode, min_photons=min_photons)
        if summed is None:
            out.append(None)
            continue

        res = None
        if method == "mle":
            try:
                res = mle_localize_image(summed, ds.pixel_size, F_EM=F_EM)
            except Exception:
                res = None
            if res is not None:
                res["x_nm"] += roi[1].start * ds.pixel_size
                res["y_nm"] += roi[0].start * ds.pixel_size
        if res is None:
            res = _centroid_result(summed, roi, ds.pixel_size)

        if res is not None:
            res["n_isolated_frames"] = n_iso
        out.append(res)
    return out


def oracle_separation_nm(ds, pair_idx: int, half_px: int = 6, **kw) -> float:
    """oracle 估计出的 pair 间距 (nm); 任一发射体定位失败返回 nan。"""
    pos = oracle_positions_nm(ds, pair_idx, half_px=half_px, **kw)
    if len(pos) < 2 or pos[0] is None or pos[1] is None:
        return float("nan")
    return float(np.hypot(pos[0]["x_nm"] - pos[1]["x_nm"],
                          pos[0]["y_nm"] - pos[1]["y_nm"]))


def oracle_emitter_error_nm(ds, pair_idx: int, half_px: int = 6, **kw) -> np.ndarray:
    """oracle 定位相对真值的位移误差 (nm), 每个发射体一个数。"""
    idx = ds.pair_emitters(pair_idx)
    pos = oracle_positions_nm(ds, pair_idx, half_px=half_px, **kw)
    errs = []
    for k, p in enumerate(pos):
        if p is None:
            errs.append(np.nan)
            continue
        gt = ds.gt_xyz[idx[k], :2]
        errs.append(float(np.hypot(p["x_nm"] - gt[0], p["y_nm"] - gt[1])))
    return np.asarray(errs, dtype=float)


# --------------------------------------------------------------------------- #
# 批量汇总
# --------------------------------------------------------------------------- #
def summarize_oracle(datasets, half_px: int = 6, n_pairs: int | None = None,
                     **kw) -> list[dict]:
    """对一组数据集逐 pair 跑 oracle, 汇总成表 (每个数据集一行)。

    含: separation_nm / n_pairs / n_pair_used / n_valid / sep_est_mean_nm /
    sep_est_std_nm / sep_bias_nm / sep_rmse_nm / emitter_err_mean_nm /
    emitter_err_rmse_nm / lp_mean_nm (MLE 给出的理论定位精度)。
    注意 ``sep_bias`` 在 d→0 时天然为正 (无噪声时 |估| 也会因噪声被抬高)。

    参数
    ----
    n_pairs : 只用前 n_pairs 个 pair (None = 全部), 方便和 level1 用同一子集对比。
    """
    rows = []
    for ds in datasets:
        npair = ds.n_pairs if n_pairs is None else min(int(n_pairs), ds.n_pairs)
        est, truth, eerr, lp = [], [], [], []
        for p in range(npair):
            pos = oracle_positions_nm(ds, p, half_px=half_px, **kw)
            idx = ds.pair_emitters(p)
            truth.append(ds.pair_gt_separation_nm(p))

            if len(pos) >= 2 and pos[0] is not None and pos[1] is not None:
                est.append(float(np.hypot(pos[0]["x_nm"] - pos[1]["x_nm"],
                                          pos[0]["y_nm"] - pos[1]["y_nm"])))
            else:
                est.append(np.nan)

            e = [float(np.hypot(q["x_nm"] - ds.gt_xyz[idx[k], 0],
                                q["y_nm"] - ds.gt_xyz[idx[k], 1]))
                 for k, q in enumerate(pos) if q is not None]
            eerr.append(float(np.mean(e)) if e else np.nan)

            lps = [q["lp_nm"] for q in pos
                   if q is not None and np.isfinite(q.get("lp_nm", np.nan))]
            lp.append(float(np.mean(lps)) if lps else np.nan)

        est = np.asarray(est, float)
        truth = np.asarray(truth, float)
        eerr = np.asarray(eerr, float)
        lp = np.asarray(lp, float)
        ok = np.isfinite(est)
        diff = est[ok] - truth[ok]
        rows.append({
            "kind": ds.kind,
            "separation_nm": ds.separation_nm,
            "n_pairs": ds.n_pairs,
            "n_pair_used": npair,
            "n_valid": int(ok.sum()),
            "sep_est_mean_nm": float(np.mean(est[ok])) if ok.any() else np.nan,
            "sep_est_std_nm": float(np.std(est[ok])) if ok.any() else np.nan,
            "sep_bias_nm": float(np.mean(diff)) if ok.any() else np.nan,
            "sep_rmse_nm": float(np.sqrt(np.mean(diff ** 2))) if ok.any() else np.nan,
            "emitter_err_mean_nm": (float(np.nanmean(eerr))
                                    if np.any(np.isfinite(eerr)) else np.nan),
            "emitter_err_rmse_nm": (float(np.sqrt(np.nanmean(eerr ** 2)))
                                    if np.any(np.isfinite(eerr)) else np.nan),
            "lp_mean_nm": float(np.nanmean(lp)) if np.any(np.isfinite(lp)) else np.nan,
        })
    return rows
