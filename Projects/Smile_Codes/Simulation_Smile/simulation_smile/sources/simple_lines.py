from __future__ import annotations
import numpy as np
from .base import EmitterSource, sample_size
from ..datatypes import EmitterSet


def _hard_core_thin(s: np.ndarray, gap: float, rng) -> np.ndarray:
    """Matérn-II 型硬核细化(沿 1D 线): 任意两个保留点间距 >= gap。

    每条线上标记不可能无限堆叠: 真实的免疫荧光标记中, 一个表位/位点通常只能
    结合一个抗体, 两个标记不会落在几乎同一位置。这里给每个泊松候选一个随机
    优先级(越小越优先), 只保留"在其 gap 邻域内优先级最低"的点 —— 无偏、
    不偏向线段某端, 且保留点之间间距严格 >= gap。
    """
    n = s.size
    if n < 2 or gap <= 0:
        return s
    order = np.argsort(s)
    s_s = s[order]
    rank = rng.permutation(n)          # 随机优先级 0..n-1 (越小越先保留)
    keep = np.zeros(n, bool)
    for i in range(n):
        si, ri = s_s[i], rank[i]
        ok = True
        j = i - 1                       # 向左: 在 gap 邻域内找更低优先级
        while j >= 0 and si - s_s[j] <= gap:
            if rank[j] < ri:
                ok = False
                break
            j -= 1
        if ok:
            j = i + 1                   # 向右
            while j < n and s_s[j] - si <= gap:
                if rank[j] < ri:
                    ok = False
                    break
                j += 1
        keep[i] = ok
    return s_s[keep]


def _clip_t(a0, lo, hi, u):
    """一维参数化裁剪: a(t) = a0 + u·t 须落在 [lo, hi]。

    返回 (ok, t_lo, t_hi); |u|≈0(与轴平行)时由 a0 是否落在区间内决定。
    """
    if abs(u) < 1e-12:
        return (lo <= a0 <= hi), -np.inf, np.inf
    t1, t2 = (lo - a0) / u, (hi - a0) / u
    return True, min(t1, t2), max(t1, t2)


def _line_fov_interval(cx, cy, ux, uy, W, H):
    """无限直线 (cx,cy) + t·(ux,uy) 与 FOV 矩形 [0,W]×[0,H] 的交集 t 区间。

    u=(ux,uy) 为单位向量 → t 即沿线的弧长 (nm)。返回 (ok, t_lo, t_hi)。
    """
    okx, x_lo, x_hi = _clip_t(cx, 0.0, W, ux)
    if not okx:
        return False, 0.0, 0.0
    oky, y_lo, y_hi = _clip_t(cy, 0.0, H, uy)
    if not oky:
        return False, 0.0, 0.0
    t_lo, t_hi = max(x_lo, y_lo), min(x_hi, y_hi)
    return (t_lo <= t_hi), float(t_lo), float(t_hi)


class SimpleLines(EmitterSource):
    """穿过给定中心、**自动截断到 FOV** 的直线上的荧光标记, 全部落在 z=0。

    lines: (n_lines, 3) 或 (n_lines, 4), 每行
        [angle, cx, cy]           → width 默认 0 (细线);
        [angle, cx, cy, width]    → width: 标记在线条带内的横向随机范围 (nm)。
        angle : 线方向角, 从 +x 轴逆时针 (默认度; angle_units='rad' 用弧度);
        cx,cy : 直线经过的中心点 (nm)。

    ★ 不需要输入线长: 每条"直线"沿 angle 穿过 (cx,cy) 无限延伸, 自动与
      FOV 矩形 [0,W]×[0,H] 求交并截断, 只在被 FOV 截出的线段上放样标记
      (放样长度 = 线段在 FOV 内的长度; 中心在 FOV 外但直线仍穿过 FOV 的
      情况同样正确)。

    density_per_nm: 沿线的**名义**标记密度 (个/nm)。每条线独立,
        候选标记数 n_i ~ Poisson(density_per_nm · L_in_fov_i),
        沿线弧长均匀随机, 横向偏移 ~ U(-width/2, width/2)。

    ★ 标记建模(贴近真实免疫荧光标记):
      真实抗体标记是"随机点过程 + 位点饱和"——一个表位只能结合一个抗体,
      同一位置不可能被堆叠多个标记; 纯 Poisson/均匀放样会在局部出现
      "过密/同位点重叠"的伪影。因此默认加**硬核最小间距**:
        min_spacing_nm=None → 自动取 平均间距/8 ≈ 1/(8·density) nm,
            只排除"几乎同一点"的堆叠, 保留泊松式随机性(实现密度≈名义密度,
            误差 ~±15%);
        min_spacing_nm=显式值(如抗体 ~10 nm) → 严格保证同一条线上任意两个
            标记的弧长间距 >= 该值(实现密度会被硬核天然上限约束, 低于名义值)。
      细化用 Matérn-II 型随机优先级吸附(见 _hard_core_thin), 无偏向。

    size: 同 Random2D, (mean_nm, variance_nm²) 或 None。
    """
    def __init__(self, fov_size_nm: tuple[float, float],
                 lines, density_per_nm: float,
                 brightness: float = 5000.0, size=None,
                 min_spacing_nm: float | None = None,
                 angle_units: str = "deg"):
        fov = np.asarray(fov_size_nm, dtype=float)
        if fov.size < 2 or np.any(fov <= 0):
            raise ValueError("fov_size_nm 应为正的 (W_nm, H_nm)")
        self.fov_size = fov                                   # (W_nm, H_nm)
        arr = np.atleast_2d(np.asarray(lines, dtype=float))
        if arr.ndim != 2 or arr.shape[1] < 3:
            raise ValueError("lines 应为 (n_lines, 3): [angle, cx, cy]"
                             " 或 (n_lines, 4): [angle, cx, cy, width]")
        if arr.shape[1] >= 4:
            widths = np.maximum(arr[:, 3], 0.0)
        else:
            widths = np.zeros(arr.shape[0])                   # 细线
        if angle_units not in ("deg", "rad"):
            raise ValueError("angle_units 只能为 'deg' 或 'rad'")
        ang = np.deg2rad(arr[:, 0]) if angle_units == "deg" else arr[:, 0]
        self.angle_units = angle_units
        self.angles = ang                                     # 弧度 (n,)
        self.cx = arr[:, 1]
        self.cy = arr[:, 2]
        self.widths = widths
        self.n_lines = int(arr.shape[0])
        self.density = float(density_per_nm)                  # 个/nm
        self.min_spacing = None if min_spacing_nm is None else float(min_spacing_nm)
        self.brightness = brightness
        self.size = None if size is None else (float(size[0]), float(size[1]))

    def sample(self, rng):
        lam = self.density
        if lam <= 0.0 or self.n_lines == 0:
            return EmitterSet(xyz=np.empty((0, 3)), brightness=np.empty(0))
        # 硬核间距: 同一条线上两标记(弧长)的最小间隔
        gap = self.min_spacing if self.min_spacing is not None else (0.125 / lam)
        W, H = float(self.fov_size[0]), float(self.fov_size[1])

        xy_parts = []
        for i in range(self.n_lines):
            th, cx, cy, wi = self.angles[i], self.cx[i], self.cy[i], self.widths[i]
            ux, uy = np.cos(th), np.sin(th)                   # 轴向单位向量
            ok, t_lo, t_hi = _line_fov_interval(cx, cy, ux, uy, W, H)
            if not ok:
                continue                                      # 该直线不穿过 FOV
            Lseg = t_hi - t_lo                                # FOV 内线段长度
            n_cand = int(rng.poisson(lam * Lseg))
            if n_cand <= 0:
                continue
            s = rng.uniform(0.0, Lseg, size=n_cand)           # 沿线段弧长
            if gap > 0.0:
                s = _hard_core_thin(s, gap, rng)              # 去同位点堆叠
            if s.size == 0:
                continue
            # 横向偏移 ~ U(-w/2, w/2) —— 让标记分布在宽度内
            tt = (rng.uniform(-0.5, 0.5, size=s.size) * wi) if wi > 0 else np.zeros(s.size)
            x = cx + (t_lo + s) * ux - tt * uy
            y = cy + (t_lo + s) * uy + tt * ux
            xy_parts.append(np.column_stack([x, y]))
        xy = np.vstack(xy_parts) if xy_parts else np.empty((0, 2))
        n = xy.shape[0]
        xyz = np.column_stack([xy, np.zeros(n)])
        brightness = np.full(n, self.brightness)
        sizes = sample_size(rng, n, self.size)
        return EmitterSet(xyz=xyz, brightness=brightness, size=sizes)

    def __repr__(self):
        return (f"SimpleLines(n_lines={self.n_lines}, "
                f"density={self.density:.3g}/nm, "
                f"min_spacing={self.min_spacing if self.min_spacing is not None else 'auto'})")
