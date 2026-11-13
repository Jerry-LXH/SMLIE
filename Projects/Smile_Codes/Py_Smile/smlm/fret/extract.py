# smlm/fret/extract.py
"""在分子位置取双通道 trace: 孔径积分 + 圆环背景。

积分的**中心**约定
------------------
孔径中心一律来自**参考通道(红/受体)的定位**, 再经配准搬到其他通道::

    红 emitter 质心(红图像坐标) --geom.to_reference(·,'r')--> 绿图像坐标

也就是"用红光定位映射的绿光点"。这样两个通道的孔径指向**同一个物理位置**,
不会因为某一帧绿通道没检出就漂。

``sigma`` 逐分子、逐通道
------------------------
匹配滤波的核权重依赖 σ, 所以半径是 ``k_sigma·σ_i``。σ 取"离该中心最近的若干
定位点的中位 σ" —— 逐分子, 不是全幅常数。
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree

from ..datatypes import Trace
from .aperture import ApertureSpec, Aperture, crop_centered, sigma_per_molecule
from .background import (AnnulusBackground, valid_mask, model_weights,
                         emitter_frame_photons, emitter_sigma)


# ---------------------------------------------------------------------------
#  单通道
# ---------------------------------------------------------------------------
def extract_channel(image, centers_yx, *, spec: ApertureSpec, n_frames=None,
                    sigmas=None, locs=None, mask_centers_yx=None,
                    model_emitters=None, sigma_default=None, tag='',
                    verbose=False):
    """在一个通道上给每个中心取一条 trace。

    Parameters
    ----------
    image : (T, H, W) ndarray
        已做漂移校正的图像栈。
    centers_yx : (N, 2) ndarray
        ``[row(y), col(x)]`` 顺序的孔径中心, **该通道的图像坐标**。
    spec : ApertureSpec
    sigmas : (N,) ndarray, optional
        逐分子 PSF σ。不给时用 ``locs`` 就近估计。
    locs : Localizations, optional
        该通道的定位, 只用来估 σ。
    mask_centers_yx : (K,2) ndarray, optional
        背景环要挖掉的核心中心。``None`` = 用 ``centers_yx`` 自己(旧行为)。
    model_emitters : sequence, optional
        与 ``mask_centers_yx`` **一一对应**的 emitter(带 ``.loc_idx``)。
        给了之后, 除了挖掉邻居核心, 还会把每个邻居的 PSF 贡献**逐帧从环像素里
        减掉** —— 早期邻居还亮着时它们的翅膀会抬高环, 但那股光不在本分子孔径里,
        直接扣会把 I_D^pre 推成负数。

    Returns
    -------
    list[Trace]
        ``meta`` 里带 ``mode / radius / sigma / area / sw / sw2 / bg_flux /
        n_bg_pixels / bg_ref``, 复核界面直接画 ``bg_flux``。
    """
    image = np.asarray(image)
    F, H, W = image.shape
    n_frames = int(F if n_frames is None else min(int(n_frames), F))
    image = image[:n_frames]

    c = np.asarray(centers_yx, float).reshape(-1, 2)
    N = len(c)
    if N == 0:
        return []

    if sigmas is None:
        sigmas = sigma_per_molecule(locs, c, radius=2.0, default=sigma_default)
    sigmas = np.asarray(sigmas, float).reshape(-1)
    if sigmas.size != N:
        raise ValueError(f'sigmas 长度 {sigmas.size} 与中心数 {N} 不符')

    # 背景掩膜中心: 默认用自己(与旧 ApertureTraceExtractor 一致)
    self_mask = mask_centers_yx is None
    mc = c if self_mask else np.asarray(mask_centers_yx, float).reshape(-1, 2)
    tree = cKDTree(mc) if len(mc) else None

    # 邻居模型: 逐帧光子数 + 各自 σ
    model_ok = (model_emitters is not None and locs is not None
                and len(model_emitters) == len(mc))
    A_nb = sig_nb = None
    if model_ok:
        A_nb = emitter_frame_photons(model_emitters, locs, n_frames,
                                    n_expect=len(mc))
        sig_nb = emitter_sigma(model_emitters, locs,
                               default=(sigma_default or 1.7))

    radii = np.array([spec.radius_for(s) for s in sigmas], int)
    if verbose:
        print(f'  [{tag}] {N} 条 trace   {spec.report()}')
        print(f'  [{tag}] 半径 {radii.min()}~{radii.max()} px '
              f'(中位 {int(np.median(radii))}), '
              f'σ 中位 {np.median(sigmas):.2f} px, '
              f'背景中心 {len(mc)} 个'
              f'{"(外部屏蔽)" if not self_mask else "(自己)"}'
              f'{"  + 环内邻居模型扣除" if model_ok else ""}')

    traces = []
    for i in range(N):
        cy, cx = c[i]
        ir, ic = int(round(cy)), int(round(cx))
        center_frac = (cy - ir, cx - ic)      # 块坐标 [half,half] 相对浮点中心的小数偏移

        ap = Aperture(spec, sigmas[i], center_frac=center_frac)
        half = ap.radius + spec.bg_radius

        box, r0, c0 = crop_centered(image, cy, cx, half)
        valid = valid_mask(half, r0, c0, (H, W))

        # 邻居在**块坐标**里的位置: 块左上角 = (ir - half, ic - half)
        nrel = None
        nidx = None
        if tree is not None and len(mc):
            qr = half + float(spec.excl_radius) + 1.0
            idx = np.asarray(tree.query_ball_point([cy, cx], qr), int)
            near = mc[idx].astype(float) if idx.size else np.zeros((0, 2))
            if len(near):
                # 去掉"自己"(最近的那个中心), 否则会把内圈环挖空
                d_self = np.hypot(near[:, 0] - cy, near[:, 1] - cx)
                j = int(np.argmin(d_self))
                if self_mask or d_self[j] <= 1.5:
                    near = np.delete(near, j, axis=0)
                    idx = np.delete(idx, j)
            if len(near):
                nidx = idx
                nrel = near + np.array([half - ir, half - ic], float)

        bg_obj = AnnulusBackground(spec, ap.radius, half,
                                   neighbours_rel=nrel, valid=valid,
                                   center_frac=center_frac)

        # 邻居模型扣除: 把它们的 PSF 贡献从环像素上减掉(而不是挖掉)
        subtract = None
        if model_ok and nidx is not None and len(nidx):
            yy, xx = np.nonzero(bg_obj.mask)
            # 注意: 不要叫 W —— 会覆盖 image.shape 解包出来的图像宽度
            kern = model_weights(np.column_stack([yy, xx]).astype(float),
                                 nrel, sig_nb[nidx])        # (M, K)
            subtract = A_nb[:, nidx].astype(float) @ kern.T  # (T, M)
        bg, bg_info = bg_obj.estimate(box, subtract=subtract)

        sl = slice(half - ap.radius, half + ap.radius + 1)
        sig_box = box[:, sl, sl]

        raw = np.einsum('tij,ij->t', sig_box, ap.w) / ap.sw2
        flux = (np.einsum('tij,ij->t', sig_box, ap.w) - bg * ap.sw) / ap.sw2
        bg_flux = bg * ap.sw / ap.sw2

        traces.append(Trace(
            emitter_id=i,
            frames=np.arange(n_frames),
            intensity=flux,
            raw=raw,
            bg=bg,
            x=float(cx), y=float(cy),
            meta={'mode': spec.mode, 'radius': int(ap.radius),
                  'sigma': float(sigmas[i]), 'area': ap.area,
                  'sw': ap.sw, 'sw2': ap.sw2,
                  'bg_flux': bg_flux, 'bg_ref': spec.bg_ref,
                  'bg_smooth': int(spec.bg_smooth), 'bg_drift': bg_info.get('drift', {}),
                  'bg_model': bool(bg_info.get('model', False)),
                  'bg_sub': float(bg_info.get('sub_mean', 0.0)),
                  'bg_sub_rel': float(bg_info.get('sub_rel', 0.0)),
                  'n_bg_pixels': int(bg_info.get('n_pixels', 0)),
                  'neighbours': int(0 if nrel is None else len(nrel)),
                  'channel': tag}))

    return traces


# ---------------------------------------------------------------------------
#  多通道
# ---------------------------------------------------------------------------
@dataclass
class ChannelTraces:
    """一个通道的结果包。"""
    channel: str
    traces: list
    centers_yx: np.ndarray
    sigmas: np.ndarray

    @property
    def radii(self):
        return np.array([t.meta['radius'] for t in self.traces], int)


class PairedExtractor:
    """在一批分子位置上, 从多个通道各取一条 trace。

    中心只给**参考通道**的坐标, 其他通道由 ``geom.to_channel`` 换算。默认
    ``from_ch`` 是参考通道, 所以正常用法是"红通道定位 -> 各通道各自取"。

    Examples
    --------
    >>> ex = PairedExtractor({'r': corr_r.data, 'g': corr_g.data}, geom=geom_c)
    >>> out = ex.run(red_centers_yx, from_ch='r', locs={'g': locs_g_c, 'r': locs_r_c},
    ...              mask_centers={'g': green_own_centers})
    >>> out['g'][0].meta['radius']
    """

    def __init__(self, images, spec=None, geom=None, reference='g'):
        self.images = {str(k): np.asarray(v) for k, v in images.items()}
        self.spec = spec if isinstance(spec, ApertureSpec) else ApertureSpec(
            **(spec or {}))
        self.geom = geom
        self.reference = str(reference)

    # ---------- 坐标换算 ----------
    def map_centers(self, xy, from_ch, to_ch):
        """从 ``from_ch`` 的**图像坐标**换到 ``to_ch`` 的图像坐标。"""
        xy = np.asarray(xy, float).reshape(-1, 2)
        if from_ch == to_ch:
            return xy.copy()
        if self.geom is None:
            raise ValueError('跨通道换算需要 geom(DualViewGeometry)')
        p = (self.geom.to_reference_calib(xy, from_ch) if from_ch != self.reference
             else self.geom.to_calib_frame(xy, self.reference))
        if to_ch == self.reference:
            return self.geom.from_calib_frame(p, self.reference)
        return self.geom.from_calib_frame(
            self.geom._inv_poly(p, to_ch), to_ch)

    # ---------- 主入口 ----------
    def run(self, centers_yx, from_ch=None, *, n_frames=None, sigmas=None,
            locs=None, mask_centers=None, verbose=False):
        from_ch = from_ch or self.reference
        centers_yx = np.asarray(centers_yx, float).reshape(-1, 2)
        out = {}
        for ch, img in self.images.items():
            cyx = self.map_centers(centers_yx, from_ch, ch)
            out[ch] = extract_channel(
                img, cyx, spec=self.spec, n_frames=n_frames,
                sigmas=(sigmas or {}).get(ch) if isinstance(sigmas, dict) else sigmas,
                locs=(locs or {}).get(ch),
                mask_centers_yx=(mask_centers or {}).get(ch),
                tag=ch, verbose=verbose)
        return out
