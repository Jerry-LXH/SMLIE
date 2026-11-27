# smlm/fret/pairs.py
"""挑选合格的 FRET trace 对。

判据分两类
----------
**污染判据 (不循环)**: 只问"这条 trace 的时间/空间结构像不像一个孤立的分子对",
与 E 的大小无关 ——

``red``      整条红光 ``intensity`` 的最大值 < ``red_min``: 这个孔径里根本没受体。
``event``    :func:`~smlm.fret.events.find_step` 找不到"红降 + 绿升"台阶。
``flags``    ``EventSpec`` 的物理硬门限: 台阶前受体在发光、台阶后真的灭掉、红降显著、
             绿升显著且持续、去淬灭后供体在发光。
``spread``   **pre 窗的绿光离散度必须与光子噪声一致**。被 FRET 淬灭的供体是一个
             *恒定暗态*; 邻居闯入/闪烁会把离散度顶上去。这是**相对**判据
             (观测离散度 / 期望噪声), 所以低 E 但干净的分子不会被误杀。
``attr``     pre 窗内落在孔径里的绿光分子数 ≤ ``max_n_pre``。这是"绿光通道密、
             亮度不可靠"的直接检验; 按 ``attr_dup`` 做空间去重, 避免同一个物理
             分子被定位器拆成两条。
``neigh``    绿光背景环里被挖掉的核心数 ≤ ``neighbours_max`` (环背景可信度)。
``after``    红光灭掉以后**不再常亮** (灭后 ``post_win`` 帧里超过阈值的占空比
              ≤ ``after_duty``)。用占空比而不是 p95: p95 是"灭后有没有再亮过",
              **会闪的受体在 200 帧里几乎必然再亮一次**, 于是被误判成"孔径里有
              ̄第二个受体"; 占空比才真正区分"另一个常亮的受体"和"同一个分子在闪"。

**刻意没有** ``I_D^pre ≤ f·I_D^post`` 这一条: 那个比值就是 ``1 − E``, 拿它当质量
门限等于"按 E 筛分子", 会把低 E 的好分子全砍掉, 而且筛完就不能再拿这批数据报告
无偏的 E 分布。真需要按样品先验收紧, 在调用方自己对算出的 E 做截断, 并写明。

用法
----
>>> recs, stats = select_pairs(traces_r, traces_g, emitters_d=multiEmitters['g'])
>>> events = [r.event for r in recs if r.keep]
>>> cal = calibrate(events, traces_a=traces_r, traces_d=traces_g)
>>> res, _ = fret_efficiency(events, traces_r, traces_g, cal)
>>> attach_efficiency(recs, res)          # 把 E 写回 record, 供画廊显示
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import numpy as np

from .datatypes import CRITERIA, CRITERIA_STRICT, passed
from .events import EventSpec, find_step

__all__ = ['PairCriteria', 'PairRecord', 'select_pairs', 'noise_per_frame',
           'pre_spread', 'donor_stats', 'attach_efficiency', 'FUNNEL']

# 漏斗顺序 (也是 reasons 的排列顺序)
FUNNEL = (('event', '没找到红降+绿升台阶'),
          ('red', '漂白前窗口不够亮'),
          ('flags', '物理门限不合格'),
          ('spread', 'pre 窗绿光不平稳'),
          ('attr', '孔径内绿光分子 > 1'),
          ('neigh', '环内邻居过多'),
          ('after', '红光灭后又亮'))

# 高斯分布下 (p84 − p16) / σ
_P16P84 = 2.0 * 0.994457883209753


@dataclass(frozen=True)
class PairCriteria:
    """挑选合格 trace 对的全部阈值。

    ``red_min``
        台阶前 ``red_win`` 帧的 ``intensity`` **中位**亮度下限。受体在"漂白"之前
        必须真的在发亮; 用整条 trace 的最大值不行 —— 一个热像素或一次短暂爆闪
        就能骗过去。
    ``red_win``
        取平均(中位)的帧数, 从事件帧往前数。
    ``lag_max``
        绿光去淬灭允许的相对红光漂白的滞后帧数, 同时传给 ``EventSpec``。
    ``pre_win``
        "洁净度"检查用的台阶前窗口 (帧)。这是**诊断窗**, 与 ``EfficiencySpec.pre_win``
        (算 E 用的窗) 无关。
    ``k_spread``
        pre 窗归一化离散度上限 (1 = 纯光子噪声)。大 = 该窗里有闪烁的邻居。
    ``attr_scale`` / ``attr_px``
        归属半径 = ``attr_scale × 孔径半径``, 孔径半径缺失时退回 ``attr_px``。
    ``attr_min_locs``
        pre 窗内至少这么多定位点才算"这个绿光分子当时在发光"。
    ``attr_dup``
        两个绿光中心相距小于它就当同一个物理分子 (定位器会把一个分子拆成两条)。
    ``max_n_pre``
        pre 窗内允许的绿光分子数上限。
    ``neighbours_max``
        绿光环内邻居数上限; ``None`` = 不设限 (只报数)。
    ``post_win``
        "红光灭后又亮"检查窗 (帧)。
    ``after_duty``
        灭后 ``post_win`` 帧里亮度超过阈值的**占空比**上限。``None`` = 退回旧口径
        (灭后 p95 低于阈值), 只对完全不闪的样品合适。
    ``require``
        ``find_step`` 的物理门限里必须全部通过的那些。
    """
    red_min: float = 150.0
    red_win: int = 10
    lag_max: int = 10
    pre_win: int = 60
    k_spread: float = 2.5
    attr_scale: float = 0.8
    attr_px: float = 3.0
    attr_min_locs: int = 5
    attr_dup: float = 0.9
    max_n_pre: int = 1
    neighbours_max: int | None = None
    post_win: int = 200
    after_duty: float | None = 0.20
    require: tuple = CRITERIA_STRICT

    def report(self):
        nb = '不限' if self.neighbours_max is None else f'≤{self.neighbours_max}'
        af = ('p95 口径' if self.after_duty is None
              else f'占空比≤{self.after_duty:.0%}')
        return (f'PairCriteria(漂白前 {self.red_win} 帧中位 ≥{self.red_min:.0f}, '
                f'lag≤{self.lag_max} 帧, pre 窗 {self.pre_win} 帧 '
                f'离散度≤{self.k_spread:g}σ, '
                f'孔径 {self.attr_scale:g}×r 内 ≤{self.max_n_pre} 个, '
                f'环内邻居 {nb}, 灭后 {af})')


@dataclass
class PairRecord:
    """一条候选 trace 对的全部诊断。``keep`` 为假时 ``reasons`` 说明原因。"""
    emitter_id: int
    keep: bool = False
    reasons: list = field(default_factory=list)
    failed_flags: list = field(default_factory=list)   # require 里没通过的子判据
    event: object = None                 # StepEvent | None
    # --- 红光 ---
    IA_max: float = np.nan               # 全程 max(intensity), 仅诊断
    IA_pre_win: float = np.nan           # 漂白前 red_win 帧的中位亮度 (判据用)
    red_after: float = np.nan            # 灭后 p95(intensity), 仅诊断
    red_after_thr: float = np.nan
    red_after_duty: float = np.nan        # 灭后超过阈值的帧占比 (判据用)
    # --- 绿光洁净度 ---
    spread: float = np.nan               # pre 窗离散度 / 实测噪声 (1 = 纯噪声)
    pre_level: float = np.nan            # pre 窗中位亮度 / σ (单位 σ)
    noise_ratio: float = np.nan          # σ_光子预算 / σ_实测 (≫1 = 门限偏松)
    n_d_pre: int = -1                    # pre 窗内孔径里的绿光分子数 (去重后)
    n_d_post: int = -1                   # post 窗同上 (仅诊断)
    n_in_ap: int = -1                    # 孔径内绿光分子总数 (未按窗口筛, 仅诊断)
    d_nearest: float = np.nan            # 最近的绿光中心到孔径中心的距离
    neighbours: int = -1                 # 绿光环内被挖洞的邻居数
    radius: float = np.nan               # 归属半径
    # --- 结果 ---
    E: float = np.nan                    # attach_efficiency 填

    @property
    def j(self):
        return -1 if self.event is None else int(self.event.frame)

    def why(self):
        if self.keep:
            return 'keep'
        return ', '.join(self.reasons) or 'keep'

    def __repr__(self):
        e = '' if not np.isfinite(self.E) else f', E={self.E:.2f}'
        return (f'PairRecord(#{self.emitter_id}, {self.why()}{e}, '
                f'I_A<={self.IA_max:.0f}, n_d^pre={self.n_d_pre})')


# ---------------------------------------------------------------------------
#  逐帧噪声
# ---------------------------------------------------------------------------
def noise_per_frame(trace):
    """由 ``intensity/bg_std/sw/sw2`` 推出的每帧强度噪声 σ (photon/frame 同单位)。

    与 :class:`~smlm.trace.calculate.PoissonSNR` 用同一个噪声模型::

        total = intensity · sw2
        noise = sqrt(total + area_eff · bg_std²),  area_eff = sw²/sw2
        σ_I   = noise / sw2

    这条式子只用光子预算, 不碰 trace 的实际涨落 —— 正好可以拿它当"应该有多少
    噪声"的参照来检验 pre 窗。
    """
    I = np.asarray(trace['intensity'], float)
    sw = float(trace.get('sw') or 0.0)
    sw2 = float(trace.get('sw2') or 0.0)
    if sw2 <= 0:
        return np.full(I.shape, np.nan)
    bg_std = np.asarray(trace.get('bg_std', np.zeros(I.shape)), float)
    area_eff = sw * sw / sw2
    noise = np.sqrt(np.maximum(I, 0.0) * sw2 + area_eff * bg_std ** 2)
    return noise / sw2


def pre_spread(trace_d, j, pre_win):
    """台阶前绿光的归一化离散度与相对电平。

    参照噪声用 **pre 窗自己的实测噪声** (``robust_noise``, 逐帧差分的 MAD), 而不是
    光子预算公式: 实测里公式常常把噪声高估一个量级 (环形背景的涨落被算了两遍),
    用它当分母会让 ``k_spread`` 实际上松掉 10 倍、一条都筛不掉。

    Returns
    -------
    (spread, level, noise_ratio) : (float, float, float)
        ``spread`` = ``(p84 − p16) / (1.989·σ_实测)``, 纯噪声下 ≈ 1;
        ``level``  = ``median(intensity) / σ`` (单位 σ, 看这个窗有多暗);
        ``noise_ratio`` = ``σ_公式 / σ_实测``, ≫1 说明光子预算高估了噪声。
        窗太短时都是 ``nan``。
    """
    from ..trace.changepoint import robust_noise

    I = np.asarray(trace_d['intensity'], float)
    a, b = max(0, int(j) - int(pre_win)), int(j)
    if b - a < 5:
        return np.nan, np.nan, np.nan
    seg = I[a:b]
    ok = np.isfinite(seg)
    if ok.sum() < 5:
        return np.nan, np.nan, np.nan
    seg = seg[ok]
    sig = np.asarray(noise_per_frame(trace_d), float)[a:b][ok]
    s_bud = float(np.nanmedian(sig)) if np.isfinite(sig).any() else np.nan

    s = float(robust_noise(seg))
    # 窗内几乎没有起伏时给 σ 一个下限 (取较大值, **不**替换实测值): 否则 0 除。
    # 下限只取预算的 2%, 所以预算被高估一个量级时实测值依然生效。
    if np.isfinite(s_bud) and s_bud > 0:
        s = max(s, 0.02 * s_bud)
    if not np.isfinite(s) or s <= 0:
        return np.nan, np.nan, np.nan
    p16, p84 = np.percentile(seg, [16, 84])
    ratio = float(s_bud / s) if (np.isfinite(s_bud) and s > 0) else np.nan
    return (float((p84 - p16) / (_P16P84 * s)),
            float(np.median(seg) / s), ratio)


# ---------------------------------------------------------------------------
#  绿光归属: 孔径里 pre 窗还有几个分子在发光
# ---------------------------------------------------------------------------
def pre_bleach_level(trace_a, j, red_win):
    """台阶前 ``red_win`` 帧的**亮态**亮度 (受体在"漂白"之前到底亮不亮)。

    先取窗内中位作阈值, 再对高于阈值的帧取中位。会闪的受体在 10 帧里可能只有
    5~6 帧亮着, 普通中位会被暗帧拉低一半, 从而把亮分子判成"没受体"; 取亮帧中位
    后, 干净 trace 与普通中位几乎相同 (只差 ~0.2σ), 闪烁 trace 则不被低估。

    窗口为空 (``j <= 0``) 时返回 ``nan``。
    """
    I = np.asarray(trace_a['intensity'], float)
    a, b = max(0, int(j) - int(red_win)), int(j)
    seg = I[a:b]
    seg = seg[np.isfinite(seg)]
    if seg.size < 2:
        return np.nan
    med = float(np.median(seg))
    on = seg[seg > med]
    return float(np.median(on)) if on.size >= 3 else med


def _n_clusters(pts, r):
    """单链聚类: 返回半径 ``r`` 内的连通簇数 (小数组, O(n²) 足够)。"""
    n = len(pts)
    if n == 0:
        return 0
    left = np.ones(n, bool)
    k = 0
    for i in range(n):
        if not left[i]:
            continue
        k += 1
        left[i] = False
        stack = [i]
        while stack:
            a = stack.pop()
            d = np.hypot(pts[:, 0] - pts[a, 0], pts[:, 1] - pts[a, 1])
            nb = np.flatnonzero(left & (d <= r))
            for b in nb:
                left[b] = False
                stack.append(int(b))
    return k


def donor_stats(tree, emitters, centre_xy, j, lag, *, pre_win, post_win,
                r_attr, min_locs, dup):
    """孔径内的绿光分子归属诊断。

    Parameters
    ----------
    tree : cKDTree
        绿光 emitter 中心 (x, y) 的 KD 树。
    emitters : sequence of Emitter
        与 ``tree`` 同序的绿光 emitter。
    centre_xy : (float, float)
        孔径中心, **该通道自己的坐标** (``trace_d['x'], trace_d['y']`` 就是)。

    Returns
    -------
    dict
        ``n_pre``    pre 窗内有 ≥``min_locs`` 个定位点、且经 ``dup`` 去重后的分子数
        ``n_post``   post 窗同上 (仅诊断)
        ``n_in_ap``  孔径内绿光分子总数 (不分窗口, 未经去重, 仅诊断)
        ``d_nearest`` 最近中心到孔径中心的距离
    """
    out = {'n_pre': -1, 'n_post': -1, 'n_in_ap': -1, 'd_nearest': np.nan}
    if tree is None or not len(emitters):
        return out                       # 没给绿光表: 未知
    c = np.asarray(centre_xy, float).reshape(2)
    idx = np.asarray(tree.query_ball_point(c, float(r_attr)), int)
    out['n_in_ap'] = int(idx.size)
    if idx.size == 0:
        # 孔径里一个绿光中心都没有 = 确定没有 pre 窗在发光的分子, 不是"未知"
        out['n_pre'] = out['n_post'] = 0
        return out
    xy = np.asarray([emitters[i].centroid for i in idx], float).reshape(-1, 2)
    d = np.hypot(xy[:, 0] - c[0], xy[:, 1] - c[1])
    out['d_nearest'] = float(d.min())

    a0, a1 = max(0, int(j) - int(pre_win)), int(j)
    b0, b1 = int(j) + int(lag), int(j) + int(lag) + int(post_win)

    def _active(lo, hi):
        keep = []
        for k, i in enumerate(idx):
            fr = np.asarray(emitters[i].frame, int)
            if int(((fr >= lo) & (fr < hi)).sum()) >= int(min_locs):
                keep.append(k)
        return _n_clusters(xy[keep], float(dup)) if keep else 0

    out['n_pre'] = _active(a0, a1)
    out['n_post'] = _active(b0, b1)
    return out


# ---------------------------------------------------------------------------
#  批量
# ---------------------------------------------------------------------------
def _p95(x, min_n=10):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, 95)) if x.size >= min_n else np.nan


def _trace_radius(trace_d, crit):
    r = float(trace_d.get('radius') or 0.0)
    return crit.attr_scale * r if r > 0 else float(crit.attr_px)


def _diagnose(trace_a, trace_d, emitter_id, *, tree, emitters_d, crit, spec):
    rec = PairRecord(emitter_id=int(emitter_id))
    ia = np.asarray(trace_a['intensity'], float)
    idd = np.asarray(trace_d['intensity'], float)

    fin = ia[np.isfinite(ia)]
    rec.IA_max = float(fin.max()) if fin.size else np.nan
    rec.neighbours = int(trace_d.get('neighbours') or 0)
    rec.radius = _trace_radius(trace_d, crit)

    ev = find_step(ia, idd, None, None, spec)
    rec.event = ev
    if ev is not None:
        # find_step 只填自己知道的东西, emitter_id 由调用方补 (scan_steps 也是这么做的);
        # 漏了它 fret_efficiency 会去读 traces[-1], 静默算错。
        ev.emitter_id = int(emitter_id)
        j, lag = int(ev.frame), int(ev.lag)
        rec.IA_pre_win = pre_bleach_level(trace_a, j, crit.red_win)
        rec.spread, rec.pre_level, rec.noise_ratio = pre_spread(
            trace_d, j, crit.pre_win)
        st = donor_stats(tree, emitters_d, (trace_d['x'], trace_d['y']), j, lag,
                         pre_win=crit.pre_win, post_win=crit.post_win,
                         r_attr=rec.radius, min_locs=crit.attr_min_locs,
                         dup=crit.attr_dup)
        rec.n_d_pre, rec.n_d_post = st['n_pre'], st['n_post']
        rec.n_in_ap, rec.d_nearest = st['n_in_ap'], st['d_nearest']
        rec.red_after_thr = max(3.0 * float(ev.sig_a),
                                float(spec.off_ratio) * max(ev.IA_pre, 0.0))
        tail = ia[j + lag: j + lag + crit.post_win]
        rec.red_after = _p95(tail)
        rec.red_after_duty = (float(np.mean(tail > rec.red_after_thr))
                              if tail.size else np.nan)

    # ---- 判据 (顺序 = FUNNEL) ----
    if ev is None:
        rec.reasons.append('event')
    else:
        # 受体在"漂白"之前必须真的在发亮: 台阶前 red_win 帧的中位亮度
        if not np.isfinite(rec.IA_pre_win) or rec.IA_pre_win <= crit.red_min:
            rec.reasons.append('red')
        # 记录**全部** 7 条子判据的结果 (含没纳入 require 的 drop_d) 便于诊断
        rec.failed_flags = [k for k in CRITERIA if not ev.flags.get(k, False)]
        if any(k in set(crit.require) for k in rec.failed_flags):
            rec.reasons.append('flags')
        if np.isfinite(rec.spread) and rec.spread > crit.k_spread:
            rec.reasons.append('spread')
        if emitters_d is not None and rec.n_d_pre > crit.max_n_pre:
            rec.reasons.append('attr')
        if crit.neighbours_max is not None and rec.neighbours > crit.neighbours_max:
            rec.reasons.append('neigh')
        if crit.after_duty is None:
            if np.isfinite(rec.red_after) and rec.red_after > rec.red_after_thr:
                rec.reasons.append('after')
        elif (np.isfinite(rec.red_after_duty)
              and rec.red_after_duty > float(crit.after_duty)):
            rec.reasons.append('after')
    rec.keep = not rec.reasons
    return rec


def select_pairs(traces_a, traces_d, *, emitters_d=None, crit=None, spec=None,
                 verbose=True):
    """逐对诊断双通道 trace, 返回 ``(records, stats)``。

    Parameters
    ----------
    traces_a, traces_d : sequence of Trace
        受体(红) / 供体(绿)。**下标必须一一对应** (同一批孔径中心)。
    emitters_d : sequence of Emitter, optional
        绿光通道**自己的** emitter (``multiEmitters['g']``), 用于归属判定。
        它们的 ``centroid`` 必须与 ``trace_d['x'] / trace_d['y']`` 同坐标系
        (即同一个分析窗口)。不给就跳过 ``attr`` 判据 (会打印提示)。
    crit : PairCriteria
    spec : EventSpec, optional
        台阶扫描参数; 不给就用 ``EventSpec(lag_max=crit.lag_max)``。
        给了也只覆盖 ``sw`` / ``sw_long`` / 各 σ 门限 —— ``lag_max`` 永远以
        ``crit`` 为准, 免得两个旋钮打架。

    Returns
    -------
    records : list[PairRecord]
        与输入**等长同序**, 通过与否看 ``.keep``; 未通过的带 ``.reasons``。
    stats : dict
        ``funnel`` (每级剩余条数)、``n``、以及几个分布的中位数。
    """
    crit = crit or PairCriteria()
    # lag_max 永远以 crit 为准: 它是"事件时间窗"这个对外旋钮, 不该因为传了 spec 就失效
    spec = replace(spec or EventSpec(), lag_max=int(crit.lag_max))
    if len(traces_a) != len(traces_d):
        raise ValueError(f'两通道 trace 数不一致: {len(traces_a)} vs {len(traces_d)}')

    tree, emitters_d = None, list(emitters_d or [])
    if emitters_d:
        from scipy.spatial import cKDTree
        xy = np.asarray([e.centroid for e in emitters_d], float).reshape(-1, 2)
        tree = cKDTree(xy)

    recs = [_diagnose(traces_a[i], traces_d[i], i, tree=tree,
                      emitters_d=emitters_d or None, crit=crit, spec=spec)
            for i in range(len(traces_a))]

    ok = np.ones(len(recs), bool)
    funnel = []
    for key, label in FUNNEL:
        if key == 'attr' and not emitters_d:
            continue
        if key == 'neigh' and crit.neighbours_max is None:
            continue
        ok &= np.array([key not in r.reasons for r in recs], bool)
        funnel.append((key, label, int(ok.sum())))

    kept = [r for r in recs if r.keep]
    stats = {'n': len(recs), 'n_keep': len(kept), 'funnel': funnel,
             'criteria': crit, 'spec': spec}
    if verbose:
        _report(recs, stats, crit, emitters_d)
    return recs, stats


def _report_flags(recs, spec, crit):
    """把"物理门限不合格"拆成每一条子判据的命中数 —— 不然看不出卡在哪。"""
    if spec is None:
        return
    bad = [r for r in recs if 'flags' in r.reasons]
    if not bad:
        return
    from collections import Counter
    cnt = Counter(k for r in bad for k in r.failed_flags)
    req = set(crit.require)
    names = {'pre_a': '台阶前受体在发光', 'off_a': '台阶后真的灭掉',
             'red': '红降显著', 'grn': '绿升显著', 'hold': '绿升持续',
             'post_d': '去淬灭后供体在发光', 'drop_d': 'ΔI_D 绝对幅度'}
    print(f'  物理门限明细 (EventSpec: sw={spec.sw}, 长窗 {spec.sw_long}, '
          f'lag≤{spec.lag_max}; 红降>{spec.min_snr_a:g}σ_A, '
          f'绿升>{spec.min_snr_d:g}σ_D, 持续>{spec.min_snr_d:g}σ_D, '
          f'去淬灭后>{spec.min_snr_on_d:g}σ_D, ΔI_D>{spec.min_drop_d:.0f} photon):')
    for k in ('pre_a', 'off_a', 'red', 'grn', 'hold', 'post_d', 'drop_d'):
        if cnt.get(k):
            req_mark = '' if k in req else '   (未纳入 require, 仅统计)'
            print(f'    {k:<7} {names[k]:<16} {cnt[k]:>5} 条{req_mark}')


def _report(recs, stats, crit, emitters_d):
    print(crit.report())
    if not emitters_d:
        print('  [!] 没给 emitters_d: 跳过"孔径内绿光分子数"判据'
              '(绿光信噪比不可靠时这是最关键的一条)')
    print(f'候选 trace 对: {stats["n"]}')
    prev = None
    for k, (key, label, n) in enumerate(stats['funnel']):
        cut = '' if prev is None else f'  (-{prev - n})'
        tail = f'   <-- 合格 ({n})' if k == len(stats['funnel']) - 1 else ''
        print(f'  {k + 1}. {label:<16}: {n:>6}{cut}{tail}')
        prev = n
    hits = [(label, sum(1 for r in recs if key in r.reasons))
            for key, label in FUNNEL]
    hits = [(lab, c) for lab, c in hits if c]
    if hits:
        print('  原因分布 (可重叠): ' +
              ', '.join(f'{lab} {c}' for lab, c in hits))
    _report_flags(recs, stats.get('spec'), crit)
    kept = [r for r in recs if r.keep]
    if not kept:
        print('[!] 一条都没通过。先把 crit 放松 / 看画廊再定阈值。')
        return
    ev = [r.event for r in kept]
    print('\n合格对的关键量 (中位):')
    print(f'  漂白前 {crit.red_win} 帧 {np.nanmedian([r.IA_pre_win for r in kept]):.0f}   '
          f'I_A^max {np.median([r.IA_max for r in kept]):.0f}   '
          f'pre 离散度 {np.nanmedian([r.spread for r in kept]):.2f}σ   '
          f'pre 电平 {np.nanmedian([r.pre_level for r in kept]):.1f}σ   '
          f'噪声模型/实测 {np.nanmedian([r.noise_ratio for r in kept]):.1f}×')
    print(f'  I_A^pre {np.median([e.IA_pre for e in ev]):.0f}  '
          f'I_D^pre {np.median([e.ID_pre for e in ev]):.0f}  '
          f'I_D^post {np.median([e.ID_post for e in ev]):.0f}  '
          f'lag {np.median([e.lag for e in ev]):.0f} 帧')
    if emitters_d:
        for key, lab in (('n_in_ap', '孔径内绿光分子总能找到'),
                         ('n_d_pre', '其中 pre 窗在发光')):
            v = [max(getattr(r, key), 0) for r in kept]
            print(f'  {lab}: {np.bincount(v).tolist()}   '
                  f'(下标 = 个数, 值 = 条数)')


def attach_efficiency(records, results):
    """把 :func:`~smlm.fret.efficiency.fret_efficiency` 的 ``E`` 写回 record。

    按 ``emitter_id`` 对齐, 原地修改并返回 ``records``, 方便画廊里显示 E。
    """
    by_id = {int(x.emitter_id): float(x.E) for x in results}
    for r in records:
        r.E = by_id.get(int(r.emitter_id), np.nan)
    return records
