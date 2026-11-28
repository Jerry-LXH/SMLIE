"""Trace 的逐字段填充框架。

为什么不再是"一个函数做完所有事"
--------------------------------
取一条 trace 牵涉若干**彼此独立**的决定: 信号用多大核/什么形状、背景用环还是平面、
要不要算 SNR、邻居要不要按 PSF 模型扣掉。把它们编进一个函数, 任何一个想换都得改
那个函数; 编进一个 dataclass, 每加一个量都要动 dataclass。

这里改成: ``Trace`` 是**开放字段的 dict** (见 :class:`smlm.datatypes.Trace`),
每个阶段 (:class:`TraceFiller`) 只声明自己写哪些键 (``provides``)、需要哪些键已经
在 (``requires``), 由 :class:`TraceFillPipeline` 串起来跑。加一个新量 = 加一个类。

三个阶段动作
------------
``plan(ctx) -> int``
    不看图像。可以在黑板 ``ctx`` 上写几何量 (典型是 ``ctx.radius``), 返回自己需要
    多大的裁切块。pipeline 取所有阶段的最大值, **只裁一次**。
``fill(trace, ctx, box)``
    原地写 ``trace`` 的字段。各阶段的 ``provides`` 互不重叠, 所以原地写是安全的。

``plan`` 存在的理由: 背景环的框大小 = 信号半径 + gap + width, 而信号半径由信号阶段
算出。若阶段只能"读 trace 已有字段", 就变成"先算半径 → 再裁切 → 才能填"的三明治,
排不出一条线性顺序。让几何量走 ``ctx`` 而不是走 ``trace``, 依赖就仍然是单向的。
"""
from __future__ import annotations

import warnings
from collections.abc import Sequence

import numpy as np
from scipy.spatial import cKDTree

from ..datatypes import Trace
from .aperture import crop_centered
from .background import valid_mask

# 字段登记表: 名字 -> 'array' (T,) | 'scalar' | 'any'。
# 它替代 dataclass 的字段检查: pipeline 在开跑前校验每个阶段声明的键都在这里。
TRACE_FIELDS = {
    # --- 骨架 (pipeline 自己写) ---
    'emitter_id': 'scalar',
    'frames': 'array',
    'x': 'scalar',
    'y': 'scalar',
    'channel': 'scalar',
    # --- 信号 ---
    'sigma': 'scalar',
    'radius': 'scalar',
    'mode': 'scalar',
    'sw': 'scalar',
    'sw2': 'scalar',
    'area': 'scalar',
    'raw': 'array',
    # --- 背景 ---
    'bg': 'array',
    'bg_std': 'array',
    'bg_flux': 'array',
    'n_bg_pixels': 'scalar',
    'neighbours': 'scalar',
    'bg_ref': 'scalar',
    'bg_smooth': 'scalar',
    'bg_model': 'scalar',
    'bg_sub': 'scalar',
    'bg_sub_rel': 'scalar',
    'bg_drift': 'any',
    # --- 由 raw/bg 推出的量 ---
    'intensity': 'array',
    'snr': 'array',
    # --- 状态判定 ---
    'state': 'array',
    'fit': 'array',
    'n_states': 'scalar',
    'bleach_frame': 'scalar',
    'alive': 'array',
    'n_on': 'scalar',
    'n_switches': 'scalar',
    'n_cp': 'scalar',
    'on_frames': 'scalar',
    'on_fraction': 'scalar',
    'noise': 'scalar',
}

SKELETON_KEYS = ('emitter_id', 'frames', 'x', 'y', 'channel')


def register_field(name, kind='array'):
    """登记一个自定义字段, 供自定义 :class:`TraceFiller` 的 ``provides`` 使用。

    Parameters
    ----------
    name : str
        字段名。
    kind : 'array' | 'scalar' | 'any'
        ``array`` = 逐帧数组 (长度 = 帧数); ``scalar`` = 单个数。
        只用于文档与校验, 不强制。

    Returns
    -------
    str
        ``name``, 便于写在类体里: ``provides = (register_field('duty'),)``。
    """
    if kind not in ('array', 'scalar', 'any'):
        raise ValueError(f"kind 必须是 'array' / 'scalar' / 'any', 得到 {kind!r}")
    TRACE_FIELDS[str(name)] = kind
    return str(name)


def as_centers_yx(obj):
    """归一化 ``Emitter`` 序列或 ``(N, 2)`` 数组到 ``[row(y), col(x)]``。"""
    arr = np.asarray(obj)
    if arr.ndim == 2 and arr.shape[1] == 2 and arr.dtype != object:
        return arr.astype(float)
    return np.array([[e.centroid[1], e.centroid[0]] for e in obj], dtype=float)


# ---------------------------------------------------------------------------
#  裁切块
# ---------------------------------------------------------------------------
class TraceBox:
    """一个分子的裁切块 ``(T, n, n)``, 中心像素在 ``[half, half]``。

    ``r0`` / ``c0`` 是块的左上角在整图上的位置(可能为负, 越界处补 0)。
    """

    __slots__ = ('data', 'half', 'r0', 'c0', 'image_shape')

    def __init__(self, data, half, r0, c0, image_shape):
        self.data = data
        self.half = int(half)
        self.r0 = int(r0)
        self.c0 = int(c0)
        self.image_shape = (int(image_shape[0]), int(image_shape[1]))

    @property
    def n(self):
        return 2 * self.half + 1

    def signal(self, radius):
        """中心 ``radius`` 半径的信号子块 (view, 不复制)。"""
        r = int(radius)
        sl = slice(self.half - r, self.half + r + 1)
        return self.data[:, sl, sl]

    def valid(self):
        """块内哪些像素真的落在图里。"""
        return valid_mask(self.half, self.r0, self.c0, self.image_shape)

    def __repr__(self):
        return f'TraceBox(T={self.data.shape[0]}, n={self.n}, r0={self.r0}, c0={self.c0})'


# ---------------------------------------------------------------------------
#  逐分子黑板
# ---------------------------------------------------------------------------
class FillContext:
    """一个分子的共享输入 + 各阶段之间传递几何量的黑板。

    每分子新建一次, 不跨分子保留状态。
    """

    __slots__ = ('emitter_id', 'cy', 'cx', 'sigma', 'image_shape', 'n_frames',
                 'tag', 'radius', 'mask_yx', 'self_mask', 'mask_tree', 'notes')

    def __init__(self, emitter_id, cy, cx, sigma, image_shape, n_frames, tag='',
                 mask_yx=None, self_mask=True, mask_tree=None):
        self.emitter_id = int(emitter_id)
        self.cy = float(cy)
        self.cx = float(cx)
        self.sigma = float(sigma)
        self.image_shape = (int(image_shape[0]), int(image_shape[1]))
        self.n_frames = int(n_frames)
        self.tag = str(tag)
        self.radius = 0                    # 由信号阶段的 plan() 填
        self.mask_yx = mask_yx
        self.self_mask = bool(self_mask)
        self.mask_tree = mask_tree
        self.notes = {}

    @property
    def center_frac(self):
        """块坐标相对整数中心的小数偏移 (fy, fx)。"""
        return (self.cy - round(self.cy), self.cx - round(self.cx))

    def neighbours(self, half, excl_radius, self_radius=1.5):
        """掩膜中心里落在环附近的那些, 换算到**块坐标**。

        Returns
        -------
        (nrel, nidx) : ``(K,2)`` 块坐标 + ``(K,)`` 原索引; 没有邻居时为 ``(None, None)``。
        """
        if self.mask_tree is None or self.mask_yx is None or not len(self.mask_yx):
            return None, None
        qr = float(half) + float(excl_radius) + 1.0
        idx = np.asarray(self.mask_tree.query_ball_point([self.cy, self.cx], qr),
                         int)
        near = self.mask_yx[idx].astype(float) if idx.size else np.zeros((0, 2))
        if len(near):
            # 去掉"自己"(最近的那个中心), 否则会把内圈环挖空
            d_self = np.hypot(near[:, 0] - self.cy, near[:, 1] - self.cx)
            j = int(np.argmin(d_self))
            if self.self_mask or d_self[j] <= float(self_radius):
                near = np.delete(near, j, axis=0)
                idx = np.delete(idx, j)
        if not len(near):
            return None, None
        ir, ic = int(round(self.cy)), int(round(self.cx))
        return near + np.array([half - ir, half - ic], float), idx

    def __repr__(self):
        return (f'FillContext(#{self.emitter_id}, y={self.cy:.2f}, '
                f'x={self.cx:.2f}, σ={self.sigma:.2f}, r={self.radius})')


# ---------------------------------------------------------------------------
#  阶段基类
# ---------------------------------------------------------------------------
class TraceFiller:
    """往 ``Trace`` 里填一组字段的阶段。

    子类声明 ``provides`` (要写的键) 与 ``requires`` (需要已存在的键), 实现
    ``plan`` / ``fill``。
    """

    name = 'filler'
    provides: tuple[str, ...] = ()
    requires: tuple[str, ...] = ()

    def plan(self, ctx) -> int:
        """不看图像, 返回需要的 half; 可写 ``ctx.radius`` 等几何量。"""
        return 0

    def fill(self, trace, ctx, box) -> None:
        raise NotImplementedError

    def describe(self) -> str:
        return self.name

    def __repr__(self):
        return f'{type(self).__name__}({self.describe()})'


# ---------------------------------------------------------------------------
#  pipeline
# ---------------------------------------------------------------------------
class TraceFillPipeline:
    """在一个通道上, 对每个 emitter 跑同一组阶段。

    每个分子只裁一次块 (取所有阶段 ``plan()`` 的最大 half), 然后依次 ``fill``。

    两个中心表是**分开**的
    ----------------------
    ``centers_for_raw``
        要提取信号的分子 (一维结果与它对齐)。决定在哪里取 `raw`。
    ``centers_for_bg``
        背景环要**挖掉**的分子中心。默认等于 ``centers_for_raw`` (挖自己);
        高密度或多通道情形要传一张更完整的表 (例如本通道所有分子的中心),
        否则邻居的 PSF 会漏进背景环、把 `bg` 抬高。
        这张表最好**包含自己的中心** —— 剔除"自己"用的是"最近的那个"规则。

    σ 必须显式给
    ------------
    PSF 宽度只接受显式输入 (``sigma``), 不再从 ``locs`` 逐分子估计: 高密度时
    "中心 2px 内的定位点"会混进邻居, σ 偏大 -> 半径偏大 -> 背景环收到更多邻居光,
    而 ``matched`` 的权重核同时也被拉宽, 一处错两处错。
    典型取值来自光学参数 ``σ ≈ 0.2165·λ/NA / pixel_size``。
    """

    def __init__(self, image, centers_for_raw, *, centers_for_bg=None,
                 sigma=None, sigma_default=1.5, n_frames=None, tag=''):
        image = np.asarray(image)
        if image.ndim != 3:
            raise ValueError(f'image 必须是 (T, H, W), 得到 {image.shape}')
        self.image = image
        F, H, W = image.shape
        self.n_frames = int(F if n_frames is None else min(int(n_frames), F))
        self.image = image[:self.n_frames]
        self.image_shape = (H, W)
        self.tag = str(tag)
        self.frames = np.arange(self.n_frames)

        self.centers_for_raw = as_centers_yx(centers_for_raw).reshape(-1, 2)
        n = len(self.centers_for_raw)

        if sigma is None:
            warnings.warn(
                f'没有指定 PSF 宽度 sigma, 全部按 {float(sigma_default):g} px 处理; '
                '建议用光学参数算出 σ (≈ 0.2165·λ/NA / pixel_size) 再传入',
                stacklevel=2)
            sig = np.full(n, float(sigma_default))
        else:
            sig = np.asarray(sigma, float).reshape(-1)
            if sig.size == 1:
                sig = np.full(n, float(sig[0]))
        if sig.size != n:
            raise ValueError(f'sigma 长度 {sig.size} 与中心数 {n} 不符')
        self.sigmas = sig

        self.self_mask = centers_for_bg is None
        mc = (self.centers_for_raw if self.self_mask
              else np.asarray(centers_for_bg, float).reshape(-1, 2))
        self.centers_for_bg = mc
        self.mask_tree = cKDTree(mc) if len(mc) else None

    # ---------- 校验 ----------
    @staticmethod
    def check(fillers, base=None):
        """静态校验阶段链。

        * ``requires`` 必须已被提供 —— 否则报错;
        * ``provides`` 必须是已登记字段 —— 否则报错;
        * ``provides`` 与已有字段重复 —— **警告并覆盖** (后面的阶段写新值)。
          注意: 若前面的阶段已经用旧值算过别的量, 那些量不会自动重算。

        ``base`` 是**已经存在**的字段集合; 默认骨架字段 (新建 trace 的情形)。
        给已有 trace 追加字段时传 ``set(trace)``。
        """
        have = set(SKELETON_KEYS if base is None else base)
        for f in fillers:
            label = type(f).__name__
            missing = [k for k in f.requires if k not in have]
            if missing:
                raise ValueError(
                    f'{label}.requires {missing} 还没有任何阶段提供; '
                    f'可用: {sorted(have)}')
            unknown = [k for k in f.provides if k not in TRACE_FIELDS]
            if unknown:
                raise ValueError(f'{label}.provides 含未登记字段 {unknown}; '
                                 '请先加进 TRACE_FIELDS')
            dup = [k for k in f.provides if k in have]
            if dup:
                warnings.warn(
                    f'{label} 覆盖已存在的字段 {dup}; 若前面的阶段已经用旧值算过'
                    f'结果, 那些结果不会自动重算', UserWarning, stacklevel=3)
            have |= set(f.provides)
        return have

    # ---------- 执行 ----------
    def run(self, fillers, *, verbose=False):
        fillers = list(fillers)
        if not fillers:
            raise ValueError('至少需要一个 TraceFiller')
        available = self.check(fillers)

        n = len(self.centers_for_raw)
        if verbose:
            print(f'  [{self.tag}] {n} 条 trace, 阶段: '
                  f'{" -> ".join(type(f).__name__ for f in fillers)}')
            for f in fillers:
                print(f'  [{self.tag}]   {type(f).__name__}: {f.describe()}')

        traces = []
        for i in range(n):
            ctx = FillContext(
                i, self.centers_for_raw[i, 0], self.centers_for_raw[i, 1],
                self.sigmas[i],
                self.image_shape, self.n_frames, tag=self.tag,
                mask_yx=self.centers_for_bg, self_mask=self.self_mask,
                mask_tree=self.mask_tree)

            half = max(int(f.plan(ctx)) for f in fillers)
            box = None
            if half > 0:
                data, r0, c0 = crop_centered(self.image, ctx.cy, ctx.cx, half)
                box = TraceBox(data, half, r0, c0, self.image_shape)

            trace = Trace(emitter_id=i, frames=self.frames,
                          x=float(ctx.cx), y=float(ctx.cy), channel=self.tag)
            for f in fillers:
                f.fill(trace, ctx, box)
            traces.append(trace)

        self._check_homogeneous(traces, available)
        if verbose:
            self._report(traces)
        return traces

    def add_fields(self, traces, fillers, *, verbose=False):
        """给**已有的** traces 追加字段 (原地修改), 不重裁 ROI、不新建列表。

        ``ctx`` 从 trace 自己的 ``x / y / sigma / radius`` 重建, 图像与掩膜中心表
        沿用本 pipeline 的, 所以对同一次取 trace 追加字段不需要再传任何东西。

        Parameters
        ----------
        traces : sequence of Trace
            同一个 pipeline 产出的、字段集一致的那些 trace。
        fillers : sequence of TraceFiller
            新阶段。``requires`` 必须已在 trace 里, ``provides`` 不能与已有字段重复。

        Returns
        -------
        list[Trace]
            包含**同一批** trace 对象的列表 (对象已被原地修改)。

        Examples
        --------
        >>> fill_traces = pipe.add_fields(fill_traces, [DutyCycle()])
        """
        traces = list(traces)
        fillers = list(fillers)
        if not fillers:
            raise ValueError('至少需要一个 TraceFiller')
        if not traces:
            return traces

        keys = set(traces[0])
        for i, t in enumerate(traces):
            if set(t) != keys:
                raise ValueError(f'trace #{i} 的字段集与 #0 不一致: '
                                 f'多 {sorted(set(t) - keys)}, '
                                 f'缺 {sorted(keys - set(t))}')
        available = self.check(fillers, base=keys)

        if verbose:
            print(f'  [{self.tag}] 追加字段: '
                  f'{" -> ".join(type(f).__name__ for f in fillers)}')

        for i, trace in enumerate(traces):
            ctx = FillContext(
                trace.get('emitter_id', i), trace.get('y', 0.0),
                trace.get('x', 0.0), float(trace.get('sigma', 1.5)),
                self.image_shape,
                int(np.size(trace['frames'])) if 'frames' in trace
                else self.n_frames,
                tag=trace.get('channel', self.tag),
                mask_yx=self.centers_for_bg, self_mask=self.self_mask,
                mask_tree=self.mask_tree)
            # 已有 radius 的阶段 (如 AnnulusRing) 要能从 trace 里接着用
            ctx.radius = int(trace.get('radius') or 0)

            half = max(int(f.plan(ctx)) for f in fillers)
            box = None
            if half > 0:
                data, r0, c0 = crop_centered(self.image, ctx.cy, ctx.cx, half)
                box = TraceBox(data, half, r0, c0, self.image_shape)
            for f in fillers:
                f.fill(trace, ctx, box)

        self._check_homogeneous(traces, available)
        return traces

    @staticmethod
    def _check_homogeneous(traces, expected):
        """所有 trace 的键集合必须一致, 且覆盖所有声明过的键。"""
        if not traces:
            return
        keys = set(traces[0])
        for i, t in enumerate(traces):
            if set(t) != keys:
                raise RuntimeError(f'trace #{i} 字段集与 #0 不一致: '
                                   f'多 {sorted(set(t) - keys)}, '
                                   f'缺 {sorted(keys - set(t))}')
        missing = set(expected) - keys
        if missing:
            raise RuntimeError(f'阶段声明了但没有写入的字段: {sorted(missing)}')

    def _report(self, traces):
        if not traces:
            return
        radii = np.array([t['radius'] for t in traces], int)
        sig = np.array([t['sigma'] for t in traces], float)
        print(f'  [{self.tag}]   半径 {radii.min()}~{radii.max()} px '
              f'(中位 {int(np.median(radii))}), σ 中位 {np.median(sig):.2f} px')
        if 'n_bg_pixels' in traces[0]:
            nb = np.array([t['neighbours'] for t in traces], int)
            npx = np.array([t['n_bg_pixels'] for t in traces], int)
            print(f'  [{self.tag}]   背景中心 {len(self.centers_for_bg)} 个'
                  f'{"(外部屏蔽)" if not self.self_mask else "(自己)"}, '
                  f'环像素中位 {int(np.median(npx))}, '
                  f'邻居中位 {int(np.median(nb))}')
        if 'bg_model' in traces[0] and any(t['bg_model'] for t in traces):
            print(f'  [{self.tag}]   + 环内邻居模型扣除')


def apply_fillers(traces, fillers, *, verbose=False):
    """在**已有** trace 上跑一组"不看影像"的阶段 (原地修改)。

    与 :meth:`TraceFillPipeline.add_fields` 的区别: 这里没有 ``image``, 所以只接受
    ``plan()`` 返回 ``0`` 的阶段 —— 即 ``fill()`` 只读 trace 自己字段的那些
    (``ThresholdStates`` / ``ChangepointStates`` 之类, ``ctx`` / ``box`` 完全没用)。
    这样从磁盘读回来的 trace 不必重新载入影像就能做状态判定。

    ``intensity`` 里**不能有 NaN**: 池化文件补在尾部的 NaN 一旦进 ``np.cumsum``,
    整条 trace 会被判成一个暗段, 而且不报任何错 —— 所以这里直接拦下来。

    Parameters
    ----------
    traces : sequence of Trace
        字段集必须一致。
    fillers : sequence of TraceFiller
        要追加的阶段。
    verbose : bool
        打印每个阶段的 ``describe()``。

    Returns
    -------
    list[Trace]
        同一批对象 (已被原地修改)。

    Examples
    --------
    >>> traces = load_pooled('TRACE_pooled/SPARK_6mW/traces_pooled.npz')[0]
    >>> apply_fillers(traces, [ChangepointStates(bleach_frames=500)])
    """
    traces = list(traces)
    fillers = list(fillers)
    if not fillers:
        raise ValueError('至少需要一个 TraceFiller')
    if not traces:
        return traces

    keys = set(traces[0])
    for i, t in enumerate(traces):
        if set(t) != keys:
            raise ValueError(f'trace #{i} 的字段集与 #0 不一致: '
                             f'多 {sorted(set(t) - keys)}, '
                             f'缺 {sorted(keys - set(t))}')
    for i, t in enumerate(traces):
        y = t.get('intensity')
        if y is not None and not np.isfinite(np.asarray(y, float)).all():
            raise ValueError(
                f'trace #{i} 的 intensity 含 NaN/inf; 池化数据要先按 n_frames '
                '截断 (用 smlm.trace.io.load_pooled 读就不会踩到)')
    available = TraceFillPipeline.check(fillers, base=keys)

    if verbose:
        for f in fillers:
            print(f'  {type(f).__name__}: {f.describe()}')

    for i, trace in enumerate(traces):
        n_frames = (int(np.size(trace['frames'])) if 'frames' in trace
                    else int(np.size(trace.get('intensity', ()))))
        ctx = FillContext(
            int(trace.get('emitter_id', i)), trace.get('y', 0.0),
            trace.get('x', 0.0), float(trace.get('sigma', 1.5)), (0, 0),
            n_frames, tag=trace.get('channel', ''))
        ctx.radius = int(trace.get('radius') or 0)      # 已有几何量的阶段接着用
        needed = max(int(f.plan(ctx)) for f in fillers)
        if needed > 0:
            culprit = next(f for f in fillers if int(f.plan(ctx)) > 0)
            raise ValueError(
                f'{type(culprit).__name__}.plan() 需要 {needed} px 的影像块, '
                '但 apply_fillers 没有影像; 需要在影像上取数的阶段请用 '
                'TraceFillPipeline.run / add_fields')
        for f in fillers:
            f.fill(trace, ctx, None)

    TraceFillPipeline._check_homogeneous(traces, available)
    return traces


def to_table(traces, fields=None):
    """把 ``list[Trace]`` 堆成 ``{字段: (N,) 或 (N,T) 数组}``, 便于批量统计。"""
    traces = list(traces)
    if not traces:
        return {}
    names = list(fields) if fields is not None else list(traces[0])
    names = [k for k in names if all(k in t for t in traces)]
    out = {}
    for k in names:
        vals = [t[k] for t in traces]
        if all(np.ndim(v) == 0 for v in vals):
            out[k] = np.asarray(vals)
        else:
            out[k] = np.vstack([np.atleast_1d(v) for v in vals])
    return out
