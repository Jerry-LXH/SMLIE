# smlm/multiplex/pair.py
"""通道之间的一对一配对, 以及多文件配对集合的合并。

配对三条约束(缺一不可):
  1. 预对齐: 用中位位移把 mov 通道平移到 ref 附近 (初值可外部给);
  2. 一对一: 一个 ref 点只允许被一个 mov 点占用, 按残差从小到大抢占;
  3. 由粗到细: 容差逐级收紧, 每轮用当前匹配的中位位移更新平移量(类 ICP)。

多文件合并的坑: 亮文件能配出几千对, 暗文件只有几百对, 直接 concat 会让亮的
那个完全主导拟合。所以 ``merge_pairs(cap_per_file=...)`` 按文件配额抽样。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
from scipy.spatial import cKDTree

DEFAULT_TOLS = (15.0, 6.0, 3.0, 1.5)


@dataclass
class PairSet:
    """一批 (mov, ref) 配对点, 附带每个点的来源文件以便分文件诊断。"""
    mov: np.ndarray                       # (N,2) mov 通道 local 坐标
    ref: np.ndarray                       # (N,2) ref 通道 local 坐标
    file: np.ndarray = None               # (N,) 每对来自哪个文件
    shift: np.ndarray = None              # (2,) 配对时用的整体平移 (dy?, 见下)

    def __post_init__(self):
        self.mov = np.asarray(self.mov, float).reshape(-1, 2)
        self.ref = np.asarray(self.ref, float).reshape(-1, 2)
        if len(self.mov) != len(self.ref):
            raise ValueError(f'mov/ref 长度不一致: {len(self.mov)} vs {len(self.ref)}')
        if self.file is None:
            self.file = np.full(len(self.mov), '', dtype=object)
        else:
            self.file = np.asarray(self.file, dtype=object).reshape(-1)
            if len(self.file) != len(self.mov):
                raise ValueError('file 长度与配对点数不一致')

    def __len__(self):
        return len(self.mov)

    def __add__(self, other):
        if len(other) == 0:
            return self
        if len(self) == 0:
            return other
        return PairSet(mov=np.vstack([self.mov, other.mov]),
                       ref=np.vstack([self.ref, other.ref]),
                       file=np.concatenate([self.file, other.file]))

    def subset(self, mask):
        mask = np.asarray(mask)
        return PairSet(self.mov[mask], self.ref[mask], self.file[mask])

    def residuals(self):
        """每对的残差向量 mov - ref 与模长。"""
        d = self.mov - self.ref
        return d, np.hypot(d[:, 0], d[:, 1])

    def files(self):
        """出现过哪些文件(保序去重)。"""
        seen, out = set(), []
        for f in self.file:
            if f not in seen:
                seen.add(f)
                out.append(f)
        return out

    def __repr__(self):
        return f'PairSet(n={len(self)}, files={self.files()})'


# ---------------------------------------------------------------------------
#  坐标提取
# ---------------------------------------------------------------------------
def xy_of(obj):
    """Localizations / 数组 / (x,y) 列表 都统一成 (N,2)。"""
    if obj is None:
        return np.zeros((0, 2))
    if hasattr(obj, 'x') and hasattr(obj, 'y'):
        return np.column_stack([np.asarray(obj.x, float),
                                np.asarray(obj.y, float)])
    return np.asarray(obj, float).reshape(-1, 2)


# ---------------------------------------------------------------------------
#  单轮匹配
# ---------------------------------------------------------------------------
def match_pairs(mov, ref, tol, shift):
    """把 mov 匹配到 ref: 先按 shift 预对齐, 取残差 < tol 的最近点, 保证一对一。

    返回 (i_mov, i_ref)。一个 ref 点只留残差最小的那个 mov 点。
    """
    mov = np.asarray(mov, float).reshape(-1, 2)
    ref = np.asarray(ref, float).reshape(-1, 2)
    d, j = cKDTree(ref).query(mov + np.asarray(shift, float), k=1)
    ok = d < tol
    if not ok.any():
        return np.empty(0, int), np.empty(0, int)
    i_mov, i_ref = np.flatnonzero(ok), j[ok]
    order = np.argsort(d[ok], kind='stable')      # 残差小的优先占用目标点
    i_mov, i_ref = i_mov[order], i_ref[order]
    _, first = np.unique(i_ref, return_index=True)
    return i_mov[first], i_ref[first]


def pair_channels(mov, ref, tols=DEFAULT_TOLS, init_shift=None,
                  min_pairs=10, verbose=False, tag=''):
    """由粗到细地配对两个通道的定位点。

    Parameters
    ----------
    mov, ref : Localizations 或 (N,2) 数组
        要变换的通道 / 参考通道的定位点, 坐标是各自的 local 坐标。
    tols : 逐级收紧的容差(px)。第一级要能覆盖粗偏移残差。
    init_shift : 初始平移(px), 形如 (dx, dy)。默认用最近邻位移的中位数估。
    min_pairs : 某级配到的对数少于这个值就不再往下细化。

    Returns
    -------
    (PairSet, shift)
        PairSet 里 ``mov``/``ref`` 都是**原始坐标**(没有预先平移),
        方便直接拿去拟合 —— ``shift`` 只是搜索用的中间量。
        拟合时多项式要拟合的是 ``ref - mov`` 这个相对位移, 它自己就是粗偏移。
    """
    mov = xy_of(mov)
    ref = xy_of(ref)
    if len(mov) == 0 or len(ref) == 0:
        return PairSet(np.zeros((0, 2)), np.zeros((0, 2))), np.zeros(2)

    if init_shift is None:
        d, j = cKDTree(ref).query(mov, k=1)
        shift = np.median(ref[j] - mov, axis=0)
    else:
        shift = np.asarray(init_shift, float).reshape(2)

    i_mov, i_ref = np.empty(0, int), np.empty(0, int)
    for tol in tols:
        i_mov, i_ref = match_pairs(mov, ref, tol, shift)
        if len(i_mov) < min_pairs:
            if verbose:
                print(f'  [{tag}] tol={tol:.1f} px 只配到 {len(i_mov)} 对, 停止细化')
            break
        shift = np.median(ref[i_ref] - mov[i_mov], axis=0)
        if verbose:
            print(f'  [{tag}] tol={tol:5.1f} px -> {len(i_mov):5d} pairs, '
                  f'shift=({shift[0]:+.2f}, {shift[1]:+.2f}) px')
    if len(i_mov) < min_pairs:
        return PairSet(np.zeros((0, 2)), np.zeros((0, 2))), shift
    return PairSet(mov[i_mov], ref[i_ref], shift=shift), shift


# ---------------------------------------------------------------------------
#  多文件合并
# ---------------------------------------------------------------------------
def merge_pairs(sets, cap_per_file=None, rng=None, shuffle=True):
    """把多个文件的 PairSet 合并成一个, 可选按文件配额抽样。

    cap_per_file : 每个文件最多保留多少对。None = 不抽。
                   作用是避免"亮文件几千对 vs 暗文件几百对"导致拟合被单文件主导。
    """
    sets = [s for s in sets if s is not None and len(s) > 0]
    if not sets:
        return PairSet(np.zeros((0, 2)), np.zeros((0, 2)))
    if rng is None:
        rng = np.random.default_rng(0)

    kept = []
    for s in sets:
        if cap_per_file is not None and len(s) > int(cap_per_file):
            idx = rng.choice(len(s), int(cap_per_file), replace=False)
            idx.sort()
            s = s.subset(idx)
        kept.append(s)

    out = kept[0]
    for s in kept[1:]:
        out = out + s

    if shuffle and len(out) > 1:
        idx = rng.permutation(len(out))
        out = out.subset(idx)
    return out
