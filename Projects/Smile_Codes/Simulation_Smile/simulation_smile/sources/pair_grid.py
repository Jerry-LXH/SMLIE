from __future__ import annotations
import numpy as np
from .base import EmitterSource
from ..datatypes import EmitterSet


class PairGrid(EmitterSource):
    """把"双发射体(pair)"铺成规则网格的静态源 —— 供亚衍射极限分辨研究使用。

    每个 pair 由 ``n_per_pair`` 个**静态**发射体组成(默认 2), 同一 pair 内间距
    固定为 ``separation_nm``, 取向随机; pair 中心落在 ``grid × grid`` 网格上并加
    随机抖动, 使不同 pair 之间隔离充分(不串扰), 同时避免所有 pair 落在相同的
    像素相位上。

    ★ 几何在**构造时**用独立 ``seed`` 生成一次, :meth:`sample` 只做拷贝。
      因此位置与外部传入的 rng、光物理模型完全无关 —— 同一布局可以配上
      "闪烁 / 常亮" 等不同光物理模型做严格对照(只变一个变量)。

    参数
    ----
    fov_size_nm : (W_nm, H_nm)
        视场物理尺寸 (nm)。
    separation_nm : float
        pair 内两个发射体的间距 (nm); ``n_per_pair=1`` 时忽略(取 0)。
    grid : int
        网格边长, pair 数 = ``grid²``。
    n_per_pair : int
        每个 pair 的发射体数, 1 或 2。
    brightness : float
        内禀亮度 (光子/s), 所有发射体相同。
    jitter_frac : float
        pair 中心抖动幅度, 以网格单元尺寸为单位 (0 = 严格网格)。
    angle_deg : float | None
        固定取向 (deg); None = 每个 pair 随机取向 U[0, π)。
    seed : int
        几何专用随机种子。
    """

    def __init__(self, fov_size_nm, separation_nm: float = 0.0, grid: int = 4,
                 n_per_pair: int = 2, brightness: float = 5000.0,
                 jitter_frac: float = 0.12, angle_deg: float | None = None,
                 seed: int = 0):
        fov = np.asarray(fov_size_nm, dtype=float).ravel()
        if fov.size < 2 or np.any(fov <= 0):
            raise ValueError("fov_size_nm 应为正的 (W_nm, H_nm)")
        if int(n_per_pair) not in (1, 2):
            raise ValueError("n_per_pair 只能是 1 或 2")

        self.fov_size = fov
        self.separation_nm = float(separation_nm)
        self.grid = int(grid)
        self.n_per_pair = int(n_per_pair)
        self.brightness = float(brightness)
        self.jitter_frac = float(jitter_frac)
        self.angle_deg = None if angle_deg is None else float(angle_deg)
        self.seed = int(seed)
        self._build()

    # ------------------------------------------------------------------ 几何
    def _build(self):
        rng = np.random.default_rng(self.seed)
        g = self.grid
        W, H = float(self.fov_size[0]), float(self.fov_size[1])
        cell = min(W, H) / g
        jit = self.jitter_frac * cell
        xs = (np.arange(g) + 0.5) * W / g
        ys = (np.arange(g) + 0.5) * H / g
        # 注意遍历顺序与下面 pair_id/index 的对应关系: 先 y 后 x
        centers = np.array(
            [(x + rng.uniform(-jit, jit), y + rng.uniform(-jit, jit))
             for y in ys for x in xs], dtype=float)
        n_pairs = centers.shape[0]

        if self.angle_deg is None:
            ang = rng.uniform(0.0, np.pi, size=n_pairs)          # rad
        else:
            ang = np.full(n_pairs, np.deg2rad(self.angle_deg))

        h = 0.5 * self.separation_nm
        n = n_pairs * self.n_per_pair
        xyz = np.zeros((n, 3))
        for k in range(n_pairs):
            cx, cy = centers[k]
            base = k * self.n_per_pair
            for r in range(self.n_per_pair):
                # 2 个 -> ±h;  1 个 -> 0 (居中)
                s = (r - 0.5 * (self.n_per_pair - 1)) * self.separation_nm
                xyz[base + r] = (cx + s * np.cos(ang[k]),
                                 cy + s * np.sin(ang[k]), 0.0)

        if n and (xyz[:, 0].min() < 0 or xyz[:, 0].max() > W
                  or xyz[:, 1].min() < 0 or xyz[:, 1].max() > H):
            raise ValueError("生成的发射体超出视场: 请减小 separation_nm / jitter_frac "
                             "或增大 fov_size_nm")

        self.xyz = xyz
        self.pair_id = np.repeat(np.arange(n_pairs), self.n_per_pair)
        self.pair_role = np.tile(np.arange(self.n_per_pair), n_pairs)
        self.pair_center = centers
        self.pair_sep_nm = np.full(n_pairs, self.separation_nm)
        self.pair_angle_deg = np.rad2deg(ang)

    # ------------------------------------------------------------------ 接口
    @property
    def n(self) -> int:
        return int(self.xyz.shape[0])

    @property
    def n_pairs(self) -> int:
        return int(self.pair_center.shape[0])

    def sample(self, rng: np.random.Generator) -> EmitterSet:
        return EmitterSet(xyz=self.xyz.copy(),
                          brightness=np.full(self.n, self.brightness))

    def layout(self) -> dict:
        """真值布局 dict, 字段名与 ``two_emitter.dataset`` 的数据集 schema 对齐。"""
        return {
            "gt_xyz": self.xyz.copy(),
            "gt_brightness": np.full(self.n, self.brightness),
            "pair_id": self.pair_id.copy(),
            "pair_role": self.pair_role.copy(),
            "pair_center": self.pair_center.copy(),
            "pair_sep_nm": self.pair_sep_nm.copy(),
            "pair_angle_deg": self.pair_angle_deg.copy(),
        }

    def __repr__(self):
        return (f"PairGrid(n_pairs={self.n_pairs}, n_per_pair={self.n_per_pair}, "
                f"sep={self.separation_nm:.1f}nm, grid={self.grid})")
