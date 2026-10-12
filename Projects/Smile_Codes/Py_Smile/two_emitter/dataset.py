"""two_emitter.dataset —— 亚衍射极限双发射体数据集: 统一 schema 的保存 / 加载。

目录结构 (默认 :data:`DATASET_ROOT`)::

    two_emitter_level0/
        manifest.json                        # 顶层清单: 全局配置 + 每个数据集一行
        blink/pair_d040nm.tif|.npz|.json      # 主数据集 (闪烁 ON/OFF)
        single/single_grid.tif|.npz|.json     # 单发射体参考集 (同一光物理)

约定
----
- ``images``  : (T, H, W) uint16, 单位 **ADU**(帧在第一维, 与 ``smlm.io.RawDataset`` 一致);
  用 :meth:`TwoEmitterDataset.photons` 按 detector 参数转回光子数。
- ``gt_xyz``  : (N, 3) nm; 原点在像面左上角, x=列方向, y=行方向(与 ``FieldOfView`` 一致)。
- ``states``  : (N, T) int8; 0=on(亮), 1=dark(暗), 2=bleached(漂白)。
- ``counts``  : (N, T) 每个发射体每帧发射的光子数(理想值, 不含相机噪声)。
- 一个数据集由若干 **pair** 组成: ``pair_id`` 把发射体归属到 pair;
  ``kind='single'`` 时每个 pair 只有 1 个发射体(便于沿用同一套索引代码)。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
import os

import numpy as np
import tifffile

#: schema 版本号; 数据结构变化时递增。
SCHEMA = "two_emitter_level0/v1"

#: 数据集根目录, 可用环境变量 ``TWO_EMITTER_DATA_ROOT`` 覆盖。
DATASET_ROOT = os.environ.get(
    "TWO_EMITTER_DATA_ROOT",
    r"E:\CodingProj\Data\Simulated_Datasets\two_emitter_level0",
)

#: 亮态在 ``states`` 里的取值。
ON_STATE = 0

#: 单发射体在 ROI 内默认的半边长(像素)。4x4 网格时 cell ≈ 16 px, 取 6 可避免串到邻 pair。
DEFAULT_HALF_PX = 6


def _json_default(o):
    """numpy 类型 → 原生类型, 供 json 序列化。"""
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return str(o)


# --------------------------------------------------------------------------- #
# 保存
# --------------------------------------------------------------------------- #
def save_two_emitter_dataset(stem, *, kind, images, gt_xyz, gt_brightness,
                             pair_id, pair_role, pair_center, pair_sep_nm,
                             pair_angle_deg, states, counts, meta,
                             root=None):
    """把一个数据集写成 ``<root>/<kind>/<stem>.{tif,npz,json}`` 三件套。

    参数
    ----
    stem : str
        数据集名(不含扩展名), 例如 ``"pair_d040nm"``。
    images : (T, H, W) uint16
        ADU 图像栈。
    其余为真值 / 光物理数组, 见模块 docstring。
    meta : dict
        任意附加元信息(会与 schema/形状等自动字段合并)。

    返回 ``{"tif": ..., "npz": ..., "json": ...}``。
    """
    root = DATASET_ROOT if root is None else root
    outdir = os.path.join(root, kind)
    os.makedirs(outdir, exist_ok=True)

    tif_path = os.path.join(outdir, f"{stem}.tif")
    npz_path = os.path.join(outdir, f"{stem}.npz")
    json_path = os.path.join(outdir, f"{stem}.json")

    images = np.asarray(images)
    tifffile.imwrite(tif_path, images)

    np.savez_compressed(
        npz_path,
        gt_xyz=np.asarray(gt_xyz, dtype=np.float64),
        gt_brightness=np.asarray(gt_brightness, dtype=np.float64),
        pair_id=np.asarray(pair_id, dtype=np.int64),
        pair_role=np.asarray(pair_role, dtype=np.int64),
        pair_center=np.asarray(pair_center, dtype=np.float64).reshape(-1, 2),
        pair_sep_nm=np.asarray(pair_sep_nm, dtype=np.float64).ravel(),
        pair_angle_deg=np.asarray(pair_angle_deg, dtype=np.float64).ravel(),
        states=np.asarray(states, dtype=np.int8),
        counts=np.asarray(counts, dtype=np.int32),
    )

    meta = dict(meta)
    sep = np.asarray(meta.get("separation_nm", np.nan), dtype=float)
    meta.update({
        "schema": SCHEMA,
        "kind": kind,
        "stem": stem,
        "images_file": os.path.basename(tif_path),
        "arrays_file": os.path.basename(npz_path),
        "image_shape": list(images.shape),
        "image_dtype": str(images.dtype),
        "n_frames": int(images.shape[0]),
        "n_emitters": int(np.asarray(gt_xyz).shape[0]),
        "n_pairs": int(np.asarray(pair_sep_nm).size),
    })
    if np.isfinite(sep):
        meta["separation_nm"] = float(sep)

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False, default=_json_default)
    return {"tif": tif_path, "npz": npz_path, "json": json_path}


# --------------------------------------------------------------------------- #
# 加载
# --------------------------------------------------------------------------- #
@dataclass
class TwoEmitterDataset:
    """一组"亚衍射双发射体"仿真数据 (ADU 图像 + 真值 + 光物理状态)。"""

    kind: str
    images: np.ndarray                   # (T, H, W) uint16, ADU
    gt_xyz: np.ndarray                   # (N, 3) nm
    gt_brightness: np.ndarray            # (N,) 光子/s
    pair_id: np.ndarray                  # (N,) 每个发射体归属的 pair
    pair_role: np.ndarray                # (N,) pair 内角色 0/1
    pair_center: np.ndarray              # (P, 2) nm
    pair_sep_nm: np.ndarray              # (P,) 目标间距 nm
    pair_angle_deg: np.ndarray           # (P,) 取向角 (deg)
    states: np.ndarray                   # (N, T) int8
    counts: np.ndarray                   # (N, T) 理想光子数
    meta: dict = field(default_factory=dict)
    path: str | None = None              # 数据集 stem (不含扩展名)
    _cache: dict = field(default_factory=dict, repr=False, compare=False)

    # ---------------------------------------------------------------- 基本属性
    @property
    def n_frames(self) -> int:
        return int(self.images.shape[0])

    @property
    def n_emitters(self) -> int:
        return int(self.gt_xyz.shape[0])

    @property
    def n_pairs(self) -> int:
        return int(np.asarray(self.pair_sep_nm).size)

    @property
    def shape(self) -> tuple[int, int]:
        return (int(self.images.shape[1]), int(self.images.shape[2]))

    @property
    def pixel_size(self) -> float:
        return float(self.meta.get("pixel_size", 1.0))

    @property
    def dt(self) -> float:
        return float(self.meta.get("dt", 1.0))

    @property
    def separation_nm(self) -> float:
        """设计间距 (nm); 单发射体数据集或未记录时返回 ``nan``。"""
        v = self.meta.get("separation_nm", None)
        return float("nan") if v is None else float(v)

    @property
    def detector(self) -> dict:
        return dict(self.meta.get("detector", {}))

    # ---------------------------------------------------------------- 单位换算
    def photons(self, clip_negative: bool = True) -> np.ndarray:
        """ADU → 光子数 (float32)。按 detector 参数反解 ``adu = e·em_gain/sens + offset``。

        结果按 ``clip_negative`` 缓存 —— 分析循环里会对同一数据集反复取值。
        """
        key = ("photons", bool(clip_negative))
        if key in self._cache:
            return self._cache[key]

        det = self.detector
        em_gain = float(det.get("em_gain", 1.0) or 1.0)
        sens = float(det.get("sensitivity", det.get("conv_factor", 1.0)) or 1.0)
        offset = float(det.get("offset", 0.0) or 0.0)
        p = (self.images.astype(np.float32) - offset) * sens / em_gain
        if clip_negative:
            p = np.maximum(p, 0.0)
        self._cache[key] = p
        return p

    # ---------------------------------------------------------------- pair 几何
    def pair_emitters(self, pair_idx: int) -> np.ndarray:
        """返回该 pair 的发射体索引 (kind='single' 时长度为 1)。"""
        return np.flatnonzero(self.pair_id == pair_idx)

    def pair_gt_separation_nm(self, pair_idx: int) -> float:
        """由真值坐标实测的 pair 间距 (nm); 单发射体返回 0。"""
        idx = self.pair_emitters(pair_idx)
        if idx.size < 2:
            return 0.0
        p = self.gt_xyz[idx, :2]
        return float(np.hypot(*(p[0] - p[1])))

    def roi(self, pair_idx: int, half_px: int = DEFAULT_HALF_PX):
        """pair 中心的像素切片 ``(slice_y, slice_x)``, 保证不越界。"""
        cx, cy = self.pair_center[pair_idx]
        px = self.pixel_size
        col = cx / px - 0.5
        row = cy / px - 0.5
        h, w = self.shape
        r0 = int(np.floor(row)) - half_px
        r1 = int(np.floor(row)) + half_px + 1
        c0 = int(np.floor(col)) - half_px
        c1 = int(np.floor(col)) + half_px + 1
        return (slice(max(r0, 0), min(r1, h)), slice(max(c0, 0), min(c1, w)))

    # ---------------------------------------------------------------- 状态
    def is_on(self, on_state: int = ON_STATE) -> np.ndarray:
        """(N, T) bool: 每个发射体每帧**帧末**是否处于亮态。"""
        return self.states == on_state

    def emits(self, min_photons: int = 1) -> np.ndarray:
        """(N, T) bool: 该发射体本帧是否**真的发出光子** (``counts >= min_photons``)。

        ★ 为什么不能只看 ``states``: ``states[i, t]`` 记录的是**帧末**状态。
        一帧内完全可以发生 on→dark 跳变 —— 帧末是 dark, 但该帧前半段仍在发光。
        因此"本帧只有谁在发光"必须看 ``counts``, 用 ``states`` 会把这类帧
        误判为"搭档不发光", 从而把搭档的光子混进求和, 造成朝向 pair 中心的偏置。
        """
        return self.counts >= int(min_photons)

    def isolated_masks(self, pair_idx: int, mode: str = "counts",
                       min_photons: int = 1, on_state: int = ON_STATE):
        """level0 的"已知谁在发光"掩码: 每个发射体"只有它发光、搭档不发光"的帧。

        mode
            ``'counts'`` (默认, exact): 用真值光子数判断 —— level0 oracle 用这个;
            ``'states'``: 用帧末状态判断 (会混入帧内跳变帧, 引入朝向 pair 中心的
            偏置, 仅保留用于对比演示)。

        返回长度 = pair 内发射体数的 ``list[bool ndarray (T,)]``。
        """
        key = ("iso", int(pair_idx), mode, int(min_photons), int(on_state))
        if key in self._cache:
            return self._cache[key]

        idx = self.pair_emitters(pair_idx)
        if mode == "counts":
            on = [self.counts[i] >= int(min_photons) for i in idx]
        elif mode == "states":
            on = [self.states[i] == int(on_state) for i in idx]
        else:
            raise ValueError("mode 只能是 'counts' 或 'states'")

        masks = []
        for k, m in enumerate(on):
            others = np.zeros(self.n_frames, dtype=bool)
            for j, mj in enumerate(on):
                if j != k:
                    others |= mj
            masks.append(m & ~others)
        self._cache[key] = masks
        return masks

    def pair_emission_fraction(self, pair_idx: int) -> dict:
        """该 pair 三类帧的占比: 只有 A 亮 / 只有 B 亮 / 都不亮 / 都亮 (按 counts)。"""
        a, b = self.isolated_masks(pair_idx, mode="counts")
        on = [self.counts[i] > 0 for i in self.pair_emitters(pair_idx)]
        both = on[0] & on[1] if len(on) > 1 else np.zeros(self.n_frames, bool)
        neither = ~(on[0] | on[1]) if len(on) > 1 else ~on[0]
        n = float(self.n_frames)
        return {"only_0": a.sum() / n, "only_1": b.sum() / n,
                "both": both.sum() / n, "neither": neither.sum() / n}

    # ---------------------------------------------------------------- 构造
    @classmethod
    def load(cls, path, root=None) -> "TwoEmitterDataset":
        """从 stem(可带扩展名)或 json/npz/tif 路径加载一个数据集。"""
        stem = str(path)
        for ext in (".json", ".npz", ".tif", ".tiff"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]

        json_path = stem + ".json"
        with open(json_path, "r", encoding="utf-8") as f:
            meta = json.load(f)

        base = os.path.dirname(json_path)
        tif_path = os.path.join(base, meta.get("images_file", os.path.basename(stem) + ".tif"))
        npz_path = os.path.join(base, meta.get("arrays_file", os.path.basename(stem) + ".npz"))

        images = tifffile.imread(tif_path)
        arrays = np.load(npz_path)

        return cls(
            kind=meta.get("kind", "unknown"),
            images=images,
            gt_xyz=arrays["gt_xyz"],
            gt_brightness=arrays["gt_brightness"],
            pair_id=arrays["pair_id"],
            pair_role=arrays["pair_role"],
            pair_center=arrays["pair_center"],
            pair_sep_nm=arrays["pair_sep_nm"],
            pair_angle_deg=arrays["pair_angle_deg"],
            states=arrays["states"],
            counts=arrays["counts"],
            meta=meta,
            path=stem,
        )

    def __repr__(self) -> str:
        return (f"TwoEmitterDataset(kind={self.kind!r}, n_frames={self.n_frames}, "
                f"shape={self.shape}, n_pairs={self.n_pairs}, "
                f"sep={self.separation_nm:.1f}nm, pixel_size={self.pixel_size})")


# --------------------------------------------------------------------------- #
# 目录级工具
# --------------------------------------------------------------------------- #
def list_datasets(kind: str | None = None, root=None) -> list[str]:
    """列出 ``<root>`` 下所有数据集 stem(可只保留某一 kind)。"""
    root = DATASET_ROOT if root is None else root
    if not os.path.isdir(root):
        return []
    kinds = [kind] if kind else [
        d for d in sorted(os.listdir(root))
        if os.path.isdir(os.path.join(root, d))
    ]
    out = []
    for k in kinds:
        d = os.path.join(root, k)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".json"):
                out.append(os.path.join(d, fn[: -len(".json")]))
    return out


def load_all(kind: str | None = None, root=None) -> list[TwoEmitterDataset]:
    """加载并按键 ``separation_nm`` 排序(单发射体 kind 的间距为 0)。"""
    ds = [TwoEmitterDataset.load(p) for p in list_datasets(kind, root)]
    ds.sort(key=lambda d: (d.separation_nm, d.path or ""))
    return ds


def write_manifest(entries, config=None, root=None) -> str:
    """写 ``<root>/manifest.json``: 全局配置 + 每个数据集一行摘要。"""
    root = DATASET_ROOT if root is None else root
    os.makedirs(root, exist_ok=True)
    manifest = {
        "schema": SCHEMA,
        "config": dict(config or {}),
        "datasets": [dict(e) for e in entries],
    }
    path = os.path.join(root, "manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False, default=_json_default)
    return path


def load_manifest(root=None) -> dict:
    """读取 ``<root>/manifest.json``。"""
    root = DATASET_ROOT if root is None else root
    with open(os.path.join(root, "manifest.json"), "r", encoding="utf-8") as f:
        return json.load(f)
