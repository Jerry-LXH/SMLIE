# smlm/fret/events.py
"""「受体漂白 + 供体去淬灭」事件扫描。

不做整条 trace 的台阶分解, 只在一个窗口内找
``红降 ΔI_A(j) > 0`` 且 ``绿升 ΔI_D(j) > 0`` 同时最显著的那一帧 ``j``;
供体一侧允许滞后 0..lag_max 帧再取窗(去淬灭常滞后 1~2 s)。

为什么要有"物理前提"硬门限
--------------------------
只靠"红降绿升"四个字会被背景过扣 / 邻居闯入骗过。所以额外要求:
台阶前受体真的在亮、台阶后真的灭掉、台阶后供体真的亮着。

供体在去淬灭后常常**闪** (亮/灭交替)。所以供体的"亮态电平"一律用
**长窗内亮帧的中位数** (:func:`post_level`: 先按 ``I_D^pre + 3σ_D`` 挑出亮帧,
亮帧不足 ``sw`` 个就算没有可信亮态), 而不是某 4 帧的中位数 —— 后者会被
"恰好灭着"的窗口坑死: 实测里有占空比 0.8 的分子被判成"绿升不持续"。
"持续"本身用 ``hold`` 表达 (长窗内亮帧占比 ≥ ``hold_frac``), 不再要求连续。

注意 **不要**用 ``I_D^pre < r·I_D^post`` 当质量门限: 这个比值本身就是
``1 − E``, 拿它筛等于"按 E 筛分子", 会把低 E 的好分子全砍掉。
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ..trace.changepoint import robust_noise
from .datatypes import StepEvent, CRITERIA, CRITERIA_STRICT, passed


@dataclass(frozen=True)
class EventSpec:
    """事件扫描的全部阈值。"""
    sw: int = 4                 # 台阶前后短窗(帧)
    sw_long: int = 20           # 长窗: 在这么长的窗里统计供体的亮帧
    lag_max: int = 16           # 供体抬升允许的最大滞后(帧)
    lag_step: int = 2

    min_snr_on_a: float = 3.0   # 台阶前受体 > 3σ_A(真的在发光)
    off_ratio: float = 0.25     # 台阶后受体 < 25%(或 < 3σ_A), 真的灭掉
    min_snr_a: float = 4.0      # 红降 > 4σ_A
    min_snr_d: float = 2.0      # 绿升 > 2σ_D
    min_snr_on_d: float = 3.0   # 台阶后供体 > 3σ_D(去淬灭成功)
    min_drop_d: float = 150.0   # ΔI_D 绝对下限(photon); 只卡信噪比会让 γ 的分母是噪声
    hold_frac: float = 0.4      # 长窗内至少这么多比例的帧处于亮态 (持续 ≠ 连续)

    max_events: int = 2000      # 超过就随机抽样(固定种子)

    def report(self):
        return (f'EventSpec(sw={self.sw}, long={self.sw_long}, lag≤{self.lag_max}; '
                f'红降>{self.min_snr_a}σ, 绿升>{self.min_snr_d}σ, '
                f'亮帧占比>{self.hold_frac:.0%}, '
                f'ΔI_D>{self.min_drop_d:.0f})')


# ---------------------------------------------------------------------------
#  工具
# ---------------------------------------------------------------------------
def _roll_mean(x, w):
    """滑动均值: 第 k 项 = mean(x[k : k+w]), 长度 len(x)-w+1。"""
    c = np.concatenate(([0.0], np.cumsum(np.nan_to_num(np.asarray(x, float)))))
    return (c[w:] - c[:-w]) / w


def _med(y, a, b):
    s = np.asarray(y, float)[max(0, int(a)):max(0, int(b))]
    s = s[np.isfinite(s)]
    return float(np.median(s)) if s.size else np.nan


def post_level(y, a, b, thr, min_n):
    """``[a,b)`` 内**亮帧**的中位数, 以及亮帧个数。

    供体在受体漂白后常常闪, 直接取这 4 帧的中位数会被"恰好灭着"的窗口坑死
    (实测占空比 0.8 的分子也会被判成"绿升不持续")。所以先按 ``thr`` 挑出亮帧,
    再对亮帧取中位数; 亮帧不足 ``min_n`` 个就返回 ``nan`` (= 没有可信的亮态)。

    Returns
    -------
    (level, n_on) : (float, int)
    """
    s = np.asarray(y, float)[max(0, int(a)):max(0, int(b))]
    s = s[np.isfinite(s)]
    if s.size == 0:
        return np.nan, 0
    on = s[s > float(thr)]
    if on.size < int(min_n):
        return np.nan, int(on.size)
    return float(np.median(on)), int(on.size)


# ---------------------------------------------------------------------------
#  单条 trace
# ---------------------------------------------------------------------------
def find_step(ya, yd, sa=None, sd=None, spec: EventSpec = None):
    """在一条双通道 trace 上找最显著的红降+绿升台阶。

    Parameters
    ----------
    ya, yd : (F,) ndarray
        受体(红) / 供体(绿) 的强度 trace。
    sa, sd : float, optional
        各自的噪声 σ; 不给就用 ``robust_noise`` 现估。
    """
    spec = spec or EventSpec()
    ya = np.asarray(ya, float)
    yd = np.asarray(yd, float)
    n = ya.size
    if n != yd.size or n < spec.sw + 6:
        return None
    sa = float(robust_noise(ya) if sa is None else sa)
    sd = float(robust_noise(yd) if sd is None else sd)
    sa, sd = max(sa, 1e-6), max(sd, 1e-6)

    lo, hi = spec.sw + 2, n - spec.sw_long - spec.lag_max - 2
    if hi <= lo:
        return None

    ra, rd = _roll_mean(ya, spec.sw), _roll_mean(yd, spec.sw)
    js = np.arange(lo, hi + 1)
    bA = ra[js - spec.sw]                     # 台阶前的红
    bD = rd[js - spec.sw]                     # 台阶前的绿
    dA = bA - ra[js]                          # 红降(正 = 掉了)
    dD = np.full(js.size, -np.inf)
    for L in range(0, spec.lag_max + 1, spec.lag_step):
        dD = np.maximum(dD, rd[js + L] - bD)  # 绿升(允许滞后)
    score = np.minimum(dA / sa, dD / sd)
    score = np.where((dA > 0) & (dD > 0), score, -np.inf)
    k = int(np.argmax(score))
    if not np.isfinite(score[k]):
        return None

    # 中位数精修(±3 帧)。供体一侧用**亮帧**统计量: 先按 d_b+3σ_D 挑亮帧, 再取其中位
    thr_d = None          # 依赖 d_b, 每个候选 j 单独算
    best = None
    for j in range(max(lo, js[k] - 3), min(hi, js[k] + 3) + 1):
        a_b = _med(ya, j - spec.sw, j)
        a_a = _med(ya, j, j + spec.sw)
        d_b = _med(yd, j - spec.sw, j)
        thr_d = d_b + 3.0 * sd                      # 亮帧阈值
        d_best, l_best, n_best = -np.inf, 0, 0
        for L in range(0, spec.lag_max + 1, spec.lag_step):
            lvl, n_on = post_level(yd, j + L, j + L + spec.sw_long,
                                    thr_d, spec.sw)
            if not np.isfinite(lvl):
                continue
            v = lvl - d_b
            if v > d_best:
                d_best, l_best, n_best = v, L, n_on
        if not np.isfinite(d_best):
            # 找不到可信亮态: 退回短窗中位数, **只为定位**(判据会把这个候选否掉),
            # 这样 find_step 仍然"总能给出一个候选", 和以前的约定一致。
            for L in range(0, spec.lag_max + 1, spec.lag_step):
                v = _med(yd, j + L, j + L + spec.sw) - d_b
                if v > d_best:
                    d_best, l_best, n_best = v, L, 0
        sc = min((a_b - a_a) / sa, d_best / sd)
        if best is None or sc > best['score']:
            best = dict(score=sc, j=int(j), lag=int(l_best), IA_pre=a_b,
                        IA_post=a_a, ID_pre=d_b, ID_post=d_b + d_best,
                        n_on=int(n_best))
    if best is None:
        return None

    j, lag = best['j'], best['lag']
    # 亮态电平 (抗闪烁) 与它的纯中位数版本 (诊断: 闪烁会让后者明显偏低)
    id_blink, n_on_d = post_level(yd, j + lag, j + lag + spec.sw_long,
                                   best['ID_pre'] + 3.0 * sd, spec.sw)
    id_long = _med(yd, j + lag, j + lag + spec.sw_long)

    dIA = best['IA_pre'] - best['IA_post']
    dID = id_blink - best['ID_pre']          # nan = 没有可信亮态 -> 相关判据全 False
    dID_long = id_long - best['ID_pre']
    need_on = int(np.ceil(spec.hold_frac * spec.sw_long))
    off_a = max(3.0 * sa, spec.off_ratio * max(best['IA_pre'], 0.0))

    flags = {
        'pre_a':  bool(best['IA_pre'] > spec.min_snr_on_a * sa),
        'off_a':  bool(best['IA_post'] < off_a),
        'red':    bool(dIA > spec.min_snr_a * sa),
        'grn':    bool(np.isfinite(dID) and dID > spec.min_snr_d * sd),
        'hold':   bool(n_on_d >= max(need_on, spec.sw)),
        'post_d': bool(np.isfinite(id_blink)
                       and id_blink > spec.min_snr_on_d * sd),
        'drop_d': bool(np.isfinite(dID) and dID > spec.min_drop_d),
    }
    ev = StepEvent(emitter_id=-1, frame=j, lag=lag,
                   IA_pre=best['IA_pre'], IA_post=best['IA_post'],
                   ID_pre=best['ID_pre'], ID_post=id_blink,
                   dIA=dIA, dID=dID, dID_long=dID_long,
                   sig_a=sa, sig_d=sd, score=float(best['score']),
                   flags=flags)
    # 用 setattr 而不是构造参数: 交互式会话里 %autoreload 可能只重载了本模块,
    # 此时旧 StepEvent 不认 n_on_d, 传参会让整条流程崩掉。
    ev.n_on_d = int(n_on_d)
    return ev


# ---------------------------------------------------------------------------
#  批量
# ---------------------------------------------------------------------------
def scan_steps(traces_a, traces_d, spec: EventSpec = None, *,
               ids=None, verbose=True):
    """对一批双通道 trace 扫事件。返回 ``list[StepEvent]``(未筛选)。

    ``traces_a`` = 受体(红), ``traces_d`` = 供体(绿), 下标必须一一对应。
    """
    spec = spec or EventSpec()
    n = min(len(traces_a), len(traces_d))
    if len(traces_a) != len(traces_d):
        raise ValueError(f'两通道 trace 数不一致: {len(traces_a)} vs {len(traces_d)}')

    ids = list(range(n)) if ids is None else sorted(int(i) for i in ids)
    if len(ids) > spec.max_events:
        rng = np.random.default_rng(0)
        ids = sorted(int(v) for v in rng.choice(ids, spec.max_events, replace=False))
        if verbose:
            print(f'  (emitter 太多: 随机抽样 {spec.max_events} / {n} 条)')

    out = []
    for i in ids:
        ev = find_step(np.asarray(traces_a[i].intensity, float),
                       np.asarray(traces_d[i].intensity, float), None, None, spec)
        if ev is None:
            continue
        ev.emitter_id = int(i)
        out.append(ev)

    if verbose:
        print(f'事件扫描: {len(ids)} 条 trace  ->  {len(out)} 个候选台阶')
        _report(out, spec)
    return out


def _report(events, spec: EventSpec):
    """逐级判据的漏斗报表。"""
    if not events:
        return
    total = len(events)
    print(f'  有候选极值点          : {total}')
    acc = []
    names = {'pre_a': '台阶前受体在发光',
             'off_a': '台阶后受体灭掉',
             'red': '红降显著',
             'grn': '绿升显著',
             'hold': '绿升持续(亮帧占比)',
             'post_d': '去淬灭后供体在发光',
             'drop_d': f'绿升幅度 > {spec.min_drop_d:.0f} photon'}
    for k in CRITERIA:
        acc.append(k)
        cnt = sum(1 for e in events if passed(e, tuple(acc)))
        mark = '   <-- 合格事件' if tuple(acc) == CRITERIA else ''
        print(f'  + {names[k]:<18}: {cnt:>6}{mark}')
