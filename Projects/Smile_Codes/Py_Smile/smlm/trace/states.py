"""状态判定: 把 ``intensity`` 变成逐帧**亮度等级**, 再挑出漂白帧。

两条路线只差"怎么二值化成等级"
------------------------------
:class:`ThresholdStates`
    使用者给一个阈值, ``intensity > high`` 就是亮的 -> 只有 0/1 两个等级。
:class:`ChangepointStates`
    不给阈值, 让变点检测自己找出亮度相近的段 -> 等级数由数据决定 (0..k-1)。

**其余全部共用** :func:`summarize_states` —— 漂白判定、ON 段计数、切换计数、
拟合亮度都是同一个实现。这样两条路线的 ``n_on`` / ``bleach_frame`` 语义永远一致,
换方法不会让下游筛选条件变味。

漂白不是"状态码"
----------------
``state`` 现在表示**亮度等级** (0 = 最暗), 多级时 1..k-1 都可能出现, 所以漂白
不能再借用某个码。漂白是**时间判定**: 最后一个 ON 之后连续 OFF 至少
``bleach_frames`` 帧。它单独存在 ``bleach_frame`` (帧索引, ``-1`` = 未漂白) 和
派生掩码 ``alive`` (``False`` 段 = 已漂白) 里。

为什么放在 ``Trace`` 的字段里
-----------------------------
状态序列和它的计数都是"逐 trace"的量, 而 ``Trace`` 已经是开放字段的 dict。
写回去之后::

    keep = [i for i, t in enumerate(traces) if t['n_switches'] <= 5]

不需要再维护一条 ``results[i]`` 与 ``traces[i]`` 平行的索引; 批量统计用
:func:`smlm.trace.fill.to_table` 直接拿到 ``(N,)`` / ``(N, T)`` 数组。
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import changepoint as cp
from .fill import TraceFiller

OFF, ON = 0, 1                      # 等级 0 = 最暗(暗态), 1 = 第一亮级
STATE_FIELDS = ('state', 'fit', 'n_states', 'bleach_frame', 'alive',
                'n_on', 'n_switches', 'n_cp', 'on_frames', 'on_fraction')


def runs_of(mask):
    """连续 ``True`` 段: ``[(start, stop_exclusive), ...]``。"""
    mask = np.asarray(mask, bool).reshape(-1)
    if not mask.any():
        return []
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def _clean(on, min_on, min_off):
    """去掉太短的 ON 段、填掉太短的 OFF 段 (两轮, 直到稳定)。

    短 OFF 段通常是闪烁噪声: 分子其实一直亮着, 只是某一帧暗了一下; 把它填掉
    才不会把一次闪烁数成三次。
    """
    on = np.array(on, bool)
    for _ in range(2):
        changed = False
        for s, e in runs_of(on):
            if e - s < int(min_on):
                on[s:e] = False
                changed = True
        for s, e in runs_of(~on):
            if e - s < int(min_off):
                on[s:e] = True
                changed = True
        if not changed:
            break
    return on


def summarize_states(levels, bleach_frames, y=None):
    """把逐帧亮度等级整理成一套统一的状态字段。

    Parameters
    ----------
    levels : array_like
        ``(T,)`` 整数等级, ``0`` = 最暗, 数值越大越亮。两条路线都产出这个。
    bleach_frames : int
        **时间判定**: 最后一个 ON 之后连续 OFF 至少这么多帧 = 已漂白。
    y : array_like, optional
        原始强度; 给了就一并算出 ``fit`` (每个等级取该级均值)。

    Returns
    -------
    dict
        ``state`` / ``fit`` / ``n_states`` / ``bleach_frame`` / ``alive`` /
        ``n_on`` / ``n_switches`` / ``n_cp`` / ``on_frames`` / ``on_fraction``。

    Notes
    -----
    计数都在 ``alive`` 区间内统计 —— 漂白之后的不该算成闪烁。
    ``n_cp`` 是**等级**的变化次数, ``n_switches`` 是**二值 ON/OFF** 的变化次数;
    只有两个等级时二者相等。
    """
    levels = np.asarray(levels, np.int8).reshape(-1)
    n = levels.size
    on = levels > OFF

    idx_on = np.flatnonzero(on)
    last_on = int(idx_on[-1]) if idx_on.size else -1
    tail = n - 1 - last_on
    bleach_frame = (last_on + 1
                    if idx_on.size and tail >= int(bleach_frames) else -1)

    alive_end = bleach_frame if bleach_frame > 0 else n
    alive = np.zeros(n, bool)
    alive[:alive_end] = True

    live_levels = levels[:alive_end]
    live_on = live_levels > OFF
    n_on = len(runs_of(live_on))
    n_switches = (int(np.count_nonzero(np.diff(live_on.astype(np.int8))))
                  if alive_end > 1 else 0)
    n_cp = int(np.count_nonzero(np.diff(live_levels))) if alive_end > 1 else 0

    out = dict(
        state=levels,
        n_states=int(np.unique(live_levels).size) if alive_end else 0,
        bleach_frame=int(bleach_frame),
        alive=alive,
        n_on=int(n_on),
        n_switches=int(n_switches),
        n_cp=int(n_cp),
        on_frames=int(live_on.sum()),
        on_fraction=float(live_on.mean()) if alive_end else 0.0,
    )
    if y is not None:
        y = np.asarray(y, float).reshape(-1)
        fit = np.zeros(n)
        for lv in np.unique(levels):
            m = levels == lv
            fit[m] = y[m].mean()
        out['fit'] = fit
    return out


# ---------------------------------------------------------------------------
#  阈值法
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ThresholdStates(TraceFiller):
    """把 ``intensity`` 按使用者给的阈值二值化, 再交给 :func:`summarize_states`。

    Parameters
    ----------
    high : float
        ON 阈值 (与 ``intensity`` 同单位), 由使用者给定, 不做自适应。
    bleach_frames : int
        **时间判定**: 最后一个 ON 之后连续 OFF 至少这么多帧, 就认为分子已经漂白。
        典型值取总帧数的 1/4 左右 (例如 2000 帧数据用 500)。
    min_on : int
        短于这么多帧的 ON 段直接丢掉 (默认 1 = 不清理)。
    min_off : int
        短于这么多帧的 OFF 段填成 ON (默认 1 = 不清理)。

    Notes
    -----
    阈值法只有两个等级, 所以 ``n_states`` 恒为 1 或 2, ``n_cp == n_switches``。
    """

    name = 'threshold_states'
    provides = STATE_FIELDS
    requires = ('intensity',)

    high: float = 100.0
    bleach_frames: int = 500
    min_on: int = 1
    min_off: int = 1

    def __post_init__(self):
        if int(self.bleach_frames) < 1:
            raise ValueError(f'bleach_frames 至少为 1, 得到 {self.bleach_frames}')
        if int(self.min_on) < 1 or int(self.min_off) < 1:
            raise ValueError('min_on / min_off 至少为 1')

    def fill(self, trace, ctx, box) -> None:
        y = np.asarray(trace['intensity'], float).reshape(-1)
        on = _clean(y > float(self.high), self.min_on, self.min_off)
        trace.update(summarize_states(on.astype(np.int8), self.bleach_frames, y))

    def describe(self) -> str:
        extra = ''
        if int(self.min_on) > 1 or int(self.min_off) > 1:
            extra = f', clean on<{self.min_on}/off<{self.min_off}'
        return (f'intensity > {self.high:g} = ON, '
                f'尾部连续 OFF >= {self.bleach_frames} 帧 = 漂白{extra}')


# ---------------------------------------------------------------------------
#  变点法
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ChangepointStates(TraceFiller):
    """用变点检测找亮度段, 再按段均值聚成等级, 交给 :func:`summarize_states`。

    与 :class:`ThresholdStates` 的唯一区别是**不要求使用者给阈值**: 等级边界由
    数据定 (``penalty`` 控制"多愿意切")。所以它能处理亮度不均一的情况
    (例如同一分子在不同时间亮度略有差别), 也能识别"一个窗口里两个分子"这类
    多级异常 —— 用 ``n_states`` 筛::

        keep = select_traces(traces, n_states=(1, 2))

    Parameters
    ----------
    bleach_frames : int
        与阈值法**同名同义**: 最后一个 ON 之后连续 OFF 至少这么多帧 = 已漂白。
    penalty : float
        分割惩罚, 越大 -> 切点越少 -> 等级越少 (默认 3.0)。
    min_seg : int
        最短段长 (帧), 默认 5。天然抑制单帧 glitch。
    merge_thr : float
        段均值合并阈值 (× 稳健噪声), 默认 2.0。
    refine_iter : int
        切点精修轮数, 默认 5。
    vst : None | 'sqrt' | 'anscombe'
        方差稳定变换, 默认不做。暗分子 (亮/暗 < 5x) 建议开 ``'anscombe'``,
        亮分子开了反而会多切几刀 —— 见 :mod:`smlm.trace.changepoint`。
    """

    name = 'changepoint_states'
    provides = STATE_FIELDS + ('noise',)
    requires = ('intensity',)

    bleach_frames: int = 500
    penalty: float = 3.0
    min_seg: int = 5
    merge_thr: float = 2.0
    refine_iter: int = 5
    vst: str | None = None

    def __post_init__(self):
        if int(self.bleach_frames) < 1:
            raise ValueError(f'bleach_frames 至少为 1, 得到 {self.bleach_frames}')
        if int(self.min_seg) < 1:
            raise ValueError(f'min_seg 至少为 1, 得到 {self.min_seg}')
        cp.forward(np.zeros(2), self.vst)          # 非法 vst 在构造时就报错

    def fill(self, trace, ctx, box) -> None:
        y = np.asarray(trace['intensity'], float).reshape(-1)
        levels, _ = cp.segment(
            y, penalty=self.penalty, min_seg=self.min_seg,
            refine_iter=self.refine_iter, merge_thr=self.merge_thr,
            vst=self.vst)
        fields = summarize_states(levels, self.bleach_frames, y)
        fields['noise'] = cp.robust_noise(y)
        trace.update(fields)

    def describe(self) -> str:
        extra = f', vst={self.vst}' if self.vst else ''
        return (f'变点检测 (penalty={self.penalty:g}, min_seg={self.min_seg}, '
                f'merge={self.merge_thr:g}σ), 尾部 OFF >= '
                f'{self.bleach_frames} 帧 = 漂白{extra}')


# ---------------------------------------------------------------------------
#  筛选
# ---------------------------------------------------------------------------
def select_traces(traces, **criteria):
    """按字段筛 trace, 返回**索引列表** (可直接用于绘图/视频)。

    每条规则可以是:

    * ``(lo, hi)`` —— 闭区间, 任一端 ``None`` 表示不限;
    * 单个数 —— 要求相等;
    * 可调用对象 —— ``func(值) -> bool``。

    Examples
    --------
    >>> keep = select_traces(fill_traces, n_on=(2, 8), n_switches=(2, None))
    >>> fig = plot_trace_fields(fill_traces, indices=keep[:6], bands='state')
    """
    traces = list(traces)
    out = []
    for i, trace in enumerate(traces):
        ok = True
        for name, rule in criteria.items():
            if name not in trace:
                raise KeyError(f'trace #{i} 没有字段 {name!r}; '
                               f'可用: {sorted(trace)}')
            value = trace[name]
            if callable(rule):
                keep = bool(rule(value))
            elif isinstance(rule, tuple):
                lo, hi = rule
                keep = ((lo is None or value >= lo)
                        and (hi is None or value <= hi))
            else:
                keep = bool(value == rule)
            if not keep:
                ok = False
                break
        if ok:
            out.append(i)
    return out
