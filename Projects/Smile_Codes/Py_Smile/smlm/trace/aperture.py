# smlm/fret/aperture.py
"""孔径积分: 匹配滤波(默认) 与 硬圆盘(旧行为)。

为什么用匹配滤波
----------------
正确的光子数估计应该"把 PSF 当权重", 而不是硬圆盘求和再除以收集率::

    g_p = exp(-r_p^2 / 2σ²) / (2πσ²)          # 归一化到 Σg = 1 的高斯
    Â   = Σ w·(I_p − b) / Σ w²,    w = g

关键性质: 分子和分母被**同一个框**截断, 截断因子互相抵消 —— 于是 ``Â`` 与框大小、
与 σ 都无关。这正是我们要的: **逐分子的收集效率恒等于 1**, κ 抖动从构造上消失。

对比硬圆盘: 收集率 f(R/σ) = 1 − exp(−R²/2σ²) 随 σ 变化。实测绿通道 σ 的 p16~p84
约在 1.72~2.46, 固定 R 时 f 会差 ±8%, 而 γ 只有一个全局值, 吸收不了它 —— 直接
变成 E 的虚假展宽。匹配滤波把这一项直接消掉。

噪声代价(整平面近似)
--------------------
    σ_Â = σ_n / sqrt(Σ w²) = 2σ·sqrt(π)·σ_n
半径 R 的框拿到的匹配滤波信息比例是 1 − exp(−R²/σ²):

    R = 1.5σ -> 89.5%
    R = 2.0σ -> 98.2%     <-- 默认
    R = 2.5σ -> 99.8%

所以 R = 2σ 就够, 框可以做得很小(σ=1.7 时 9x9), 比旧的 r=5 硬圆盘还小一圈。
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def clean_sigma(sigma, default=1.5) -> float:
    """非法 σ (nan/inf/<=0) 退回 ``default``。"""
    return float(sigma) if np.isfinite(sigma) and sigma > 0 else float(default)


def radius_from_sigma(sigma, *, k_sigma=2.0, radius=None, min_radius=2,
                      max_radius=7) -> int:
    """逐分子半径: 显式 ``radius``, 否则 ``ceil(k_sigma·σ)``, 再截到上下限。"""
    if radius is not None:
        r = int(radius)
    else:
        r = int(np.ceil(float(k_sigma) * clean_sigma(sigma)))
    return int(np.clip(r, int(min_radius), int(max_radius)))


def aperture_kernel(mode, radius, sigma, center_frac=(0.0, 0.0)):
    """一个分子的信号核。

    Returns
    -------
    w : (n, n) ndarray
        权重核; ``matched`` 是归一化高斯, ``hard`` 是 0/1 圆盘。
    sw, sw2 : float
        ``w.sum()`` 与 ``(w**2).sum()`` —— 扣背景与 SNR 都靠它们换算, 所以
        **换任何新的 raw 估计方式都必须同时给出这两个量**。
    area : float
        覆盖像素数 (硬圆盘就是圆内像素数)。
    sigma : float
        清洗后的 σ。
    """
    s = clean_sigma(sigma)
    r = int(radius)
    n = 2 * r + 1
    yy, xx = np.mgrid[0:n, 0:n]
    fy, fx = float(center_frac[0]), float(center_frac[1])

    if mode == 'matched':
        # 权重核必须建在**真实的亚像素中心**上, 否则等于把 PSF 整体平移了
        # 半个像素 —— σ 小时损失可观(实测 σ=1.2 会掉 4%)。
        rr = np.hypot(yy - (r + fy), xx - (r + fx))
        g = np.exp(-rr ** 2 / (2.0 * s ** 2)) / (2.0 * np.pi * s ** 2)
        cover = rr <= r + 0.5
        w = np.where(cover, g, 0.0)
        sw2 = float((w ** 2).sum())
    elif mode == 'hard':                       # 保持旧行为, 整数中心 + 不补偿截断
        rr2 = (yy - r) ** 2 + (xx - r) ** 2
        cover = rr2 <= r ** 2 + 1e-9
        w = cover.astype(float)
        sw2 = 1.0
    else:
        raise ValueError(f"mode 必须是 'matched' 或 'hard', 得到 {mode!r}")

    return w, float(w.sum()), sw2, float(cover.sum()), s


@dataclass(frozen=True)
class ApertureSpec:
    """孔径 + 背景环的规格。

    Parameters
    ----------
    mode : 'matched' | 'hard'
        ``matched`` 匹配滤波(推荐); ``hard`` 硬圆盘, 保持旧 ``ApertureTraceExtractor``
        的行为, 只用于对照/回归。
    k_sigma : float
        半径 = ``k_sigma * σ``(逐分子)。2.0 已能拿到 98% 的匹配滤波信息。
    radius : int, optional
        显式固定半径, 给了就忽略 ``k_sigma``。
    min_radius, max_radius : int
        半径上下限。``max_radius`` 防止 σ 异常时框变大(实测 σ p95 = 3.11)。
    bg_gap : int
        信号框外让出多少 px 才开始算背景环。
    bg_width : int
        背景环宽度。
    excl_radius : int
        背景环里每个邻居被**挖掉**的半径。只挖核心(2 px), 外周仍进背景估计 ——
        挖太狠等于把最能代表本分子位置背景的那部分光删掉。
    bg_ref : 'frames' | 'tail'
        ``frames`` 逐帧估计(跟得住时间漂移); ``tail`` 只取最后 ``bg_tail`` 帧
        聚合成常数(噪声最低, 但要求背景不随时间变)。
    bg_smooth : int
        逐帧背景的**居中滑动平均**窗宽(帧)。0 = 不平滑。
        这是 ``frames`` 与 ``tail`` 之间的中间态: 仍然逐帧扣(跟得住漂移),
        但把帧间抖动按 √w 压下去。代价是会轻微抹平真实的快变。
    """
    mode: str = 'matched'
    k_sigma: float = 2.0
    radius: int | None = None
    min_radius: int = 2
    max_radius: int = 7

    bg_gap: int = 2
    bg_width: int = 4
    excl_radius: int = 2
    bg_ref: str = 'frames'
    bg_tail: int = 200
    bg_smooth: int = 0

    clip_n: float = 3.0
    clip_iter: int = 3

    def __post_init__(self):
        if self.mode not in ('matched', 'hard'):
            raise ValueError(f"mode 必须是 'matched' 或 'hard', 得到 {self.mode!r}")
        if self.bg_ref not in ('frames', 'tail'):
            raise ValueError(f"bg_ref 必须是 'frames' 或 'tail', 得到 {self.bg_ref!r}")

    # ---------- 半径 ----------
    def radius_for(self, sigma) -> int:
        """逐分子半径: ``k_sigma·σ`` 截到 [min_radius, max_radius]。"""
        return radius_from_sigma(sigma, k_sigma=self.k_sigma, radius=self.radius,
                                 min_radius=self.min_radius,
                                 max_radius=self.max_radius)

    @property
    def bg_radius(self) -> int:
        """背景环外半径相对信号半径的增量。"""
        return int(self.bg_gap) + int(self.bg_width)

    def report(self) -> str:
        if self.radius is not None:
            r = f'r={self.radius}(固定)'
        else:
            r = f'r={self.k_sigma:g}σ (σ={self.min_radius/self.k_sigma:g}~' \
                f'{self.max_radius/self.k_sigma:g} px)'
        sm = f', smooth {self.bg_smooth} 帧' if self.bg_smooth else ''
        return (f'ApertureSpec({self.mode}, {r}; '
                f'bg ring +{self.bg_gap}..+{self.bg_radius} px, '
                f'excl {self.excl_radius}, ref={self.bg_ref}{sm})')


# ---------------------------------------------------------------------------
#  单分子积分核
# ---------------------------------------------------------------------------
class Aperture:
    """一个分子在一个通道上的积分核(尺寸随 σ 变, 所以逐分子各建一个)。

    ``flux(box, bg)`` 里 ``box`` 是 (T, n, n) 的裁剪块, 中心在 ``[r, r]``。
    """

    def __init__(self, spec: ApertureSpec, sigma, center_frac=(0.0, 0.0)):
        self.spec = spec
        self.sigma = clean_sigma(sigma)
        self.radius = spec.radius_for(self.sigma)
        self.center_frac = (float(center_frac[0]), float(center_frac[1]))
        (self.w, self.sw, self.sw2,
         self.area, _) = aperture_kernel(spec.mode, self.radius, self.sigma,
                                         self.center_frac)
        self.n = 2 * self.radius + 1
        self.signal = self.w > 0 if spec.mode == 'matched' else self.w.astype(bool)

    # ---------- 积分 ----------
    def flux(self, box, bg=0.0):
        """``box`` (T, n, n) -> 逐帧光子数 (T,)。``bg`` 逐像素背景, 标量或 (T,)。"""
        box = np.asarray(box, float)
        if box.shape[1:] != (self.n, self.n):
            raise ValueError(f'box 形状 {box.shape[1:]} 与核 {self.n}x{self.n} 不符')
        s = np.einsum('tij,ij->t', box, self.w)
        b = np.asarray(bg, float)
        return (s - b * self.sw) / self.sw2

    def __repr__(self):
        return (f'Aperture({self.spec.mode}, r={self.radius}, σ={self.sigma:.2f}, '
                f'sw={self.sw:.3f}, sw2={self.sw2:.3e})')


# ---------------------------------------------------------------------------
#  工具
# ---------------------------------------------------------------------------
def crop_centered(vol, cy, cx, half):
    """从 (T,H,W) 里裁 (2*half+1)² 的小块, 越界补 0。

    返回 ``(box, r0, c0)``; ``r0/c0`` 是块左上角在图上的位置(可能为负)。
    """
    vol = np.asarray(vol)
    F, H, W = vol.shape
    n = 2 * int(half) + 1
    r0, c0 = int(round(cy)) - int(half), int(round(cx)) - int(half)
    dt = vol.dtype if vol.dtype.kind == 'f' else np.float32
    out = np.zeros((F, n, n), dt)
    rr0, rr1 = max(0, r0), min(H, r0 + n)
    cc0, cc1 = max(0, c0), min(W, c0 + n)
    if rr1 > rr0 and cc1 > cc0:
        out[:, rr0 - r0:rr1 - r0, cc0 - c0:cc1 - c0] = vol[:, rr0:rr1, cc0:cc1]
    return out, r0, c0


def sigma_per_molecule(locs, centers_yx, *, radius=2.0, default=None):
    """每个中心附近定位点的 σ 中位数(逐分子、逐通道)。

    Parameters
    ----------
    locs : Localizations
        该通道的定位, 坐标与 ``centers_yx`` 同一坐标系。需要 ``.x/.y/.sigma``。
    centers_yx : (N, 2) ndarray
        ``[row(y), col(x)]`` 顺序的孔径中心。
    radius : float
        归属半径(px)。只取落在里面的定位点。
    default : float, optional
        一个点都没有时用的 σ; 默认取整幅的中位 σ。
    """
    from scipy.spatial import cKDTree

    c = np.asarray(centers_yx, float).reshape(-1, 2)
    if len(c) == 0:
        return np.zeros(0)
    if locs is None or len(locs) == 0 or getattr(locs, 'sigma', None) is None:
        d = float(default) if default is not None else 1.5
        return np.full(len(c), d)

    xy = np.column_stack([np.asarray(locs.y, float), np.asarray(locs.x, float)])
    sg = np.asarray(locs.sigma, float)
    ok = np.isfinite(sg) & (sg > 0)
    if not ok.any():
        d = float(default) if default is not None else 1.5
        return np.full(len(c), d)
    dflt = float(default) if default is not None else float(np.median(sg[ok]))

    tree = cKDTree(xy[ok])
    sg = sg[ok]
    out = np.full(len(c), dflt)
    for i, p in enumerate(c):
        ii = tree.query_ball_point(p, float(radius))
        if ii:
            out[i] = float(np.median(sg[ii]))
    return out
