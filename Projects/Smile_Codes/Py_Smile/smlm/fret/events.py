# smlm/fret/events.py
"""「受体漂白 + 供体去淬灭」事件扫描。

受体一侧: **只把"亮段的结束"当候选台阶**
------------------------------------------
真数据里受体经常**闪** (亮/灭交替)。早期实现是"用 4 帧滑动均值算红降, 取最显著的
那一帧", 结果在会闪的样品上大面积失效: 暗隙本身也是一个"红降", 而且往往比真正的
漂白更陡, 于是 ``j`` 落在暗隙上, ``I_A^pre`` 只有半亮、``off_a`` 判失败、
``t_bleach`` 也错 (实测 348 条里有 250 条因此被判成"灭后又亮")。

现在改成: 先按 p5/p95 定出暗/亮电平, 用**迟滞**二值化得到逐帧 ON 掩码
(:func:`_on_mask`), 取每段 ON 的**结束帧**作为候选, ``I_A^pre`` 取该亮段内的中位。
候选之间按两条规则挑:

1. **优先"真的掉到背景"的候选** (短窗 ``sw`` 后已低于 ``off_ratio·I_A^pre`` 或 3σ)——
   闪烁的暗隙掉不下去, 漂白才掉得下去;
2. 其中取 ``min(持续落差/σ_A, 供体升幅/σ_D)`` 最大者, 持续落差用 ``sw_post``
   (默认 100 帧) 的长窗中位算 —— 暗隙在这么长的窗里会被"拉回来", 落差自然小。

早期实现只用台阶后 4 帧的落差打分, 于是**暗隙和真漂白在打分上没区别**
(两者 4 帧后都接近 0), 谁赢取决于供体那边谁碰巧有跳变: 实测 fov6 的 #156
真漂白在 57.7 s, 却被判到 29.6 s 的一次部分下跌上。

暗隙仍是合法的候选, 但会被 ``after`` 判据 (灭后占空比) 正确否掉 —— 各判据
只负责自己那件事, 不再靠窗口统计量互相纠缠。

供体一侧: 亮态电平用**亮帧中位**
--------------------------------
供体去淬灭允许滞后 0..lag_max 帧, 电平用 :func:`post_level` (长窗内亮帧的中位数;
亮帧不足 ``sw`` 个就算没有可信亮态), 而不是某 4 帧的中位数 —— 后者会被"恰好灭着"
的窗口坑死。"持续"用 ``hold`` 表达 (长窗内亮帧占比 ≥ ``hold_frac``), 不要求连续。

为什么要有"物理前提"硬门限
--------------------------
只靠"红降绿升"四个字会被背景过扣 / 邻居闯入骗过。所以额外要求:
台阶前受体真的在亮、台阶后真的灭掉、台阶后供体真的亮着。

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
    sw_post: int = 100          # "再也不回来"的判断窗: 用这么长窗的中位算持续落差
    pre_duty_win: int = 150     # 红光状态的统计窗: 台阶前这么多帧里 ON 的占比
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


def _on_mask(y, sd):
    """受体的迟滞二值化: 返回 (逐帧 ON 布尔, 亮态电平)。

    闪烁的受体在 4 帧中位里会被平均成"半亮", 窗口法就找不到真正的漂白点。
    所以先按 p5/p95 定出暗态与亮态电平, 再用**迟滞**二值化: 高于上阈值才转入
    亮态, 低于下阈值才转回暗态 —— 段内的单帧抖动不会切断亮段。

    亮暗分不开 (p95-p5 < 4σ) 时返回 ``None`` = 这条 trace 没有可辨的受体亮态。
    """
    p5, p95 = np.percentile(y, [5.0, 95.0])
    if not np.isfinite(p95 - p5) or p95 - p5 < 4.0 * sd:
        return None
    thr_hi = p5 + 0.55 * (p95 - p5)
    thr_lo = p5 + 0.25 * (p95 - p5)
    on = np.zeros(y.size, bool)
    state = False
    for i, v in enumerate(y):
        if state:
            state = v >= thr_lo
        else:
            state = v > thr_hi
        on[i] = state
    return on


def _episodes(on, min_len):
    """ON 段列表 ``[(start, stop)]`` (stop 开区间), 只留长度 ≥ ``min_len`` 的。"""
    edges = np.flatnonzero(np.diff(on.astype(np.int8))) + 1
    bounds = np.concatenate(([0], edges, [on.size]))
    return [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:])
            if on[a] and b - a >= int(min_len)]


# ---------------------------------------------------------------------------
#  红光状态: 平稳 / 闪烁 / 闪光的结束 / 没真灭
# ---------------------------------------------------------------------------
@dataclass
class RedState:
    """漂白事件附近受体通道的状态刻画 (全部由迟滞 ON 掩码得到)。

    ``duty``      台阶前 ``window`` 帧里 ON 的占比 (1 = 一直亮着)
    ``n_off``     同一窗口内 ≥2 帧的暗段个数 (闪了几次)
    ``level``     窗口内**亮帧**的中位亮度
    ``post_ratio`` 灭后长窗电平 / ``level`` (≈0 = 真的灭了)
    ``label``     ``'steady' | 'blinking' | 'dark' | 'unbleached' | 'unknown'``

    为什么需要它: 受体闪不闪决定了**供体的判据能有多严**。受体平稳亮着时, 供体
    在 pre 窗里也应当只有光子噪声; 受体自己在闪时, 供体跟着闪是正常的, 用同一把
    尺子量会把它们全判成"孔径里有别人"。反过来 ``dark`` 说明这个"台阶"只是一次
    闪光的结束, 根本不是漂白。
    """
    duty: float = np.nan
    n_off: int = -1
    level: float = np.nan
    post_ratio: float = np.nan
    label: str = 'unknown'

    def __repr__(self):
        return (f'RedState({self.label}, duty={self.duty:.2f}, '
                f'off={self.n_off}, post={self.post_ratio:.2f})')


def red_state(ya, j, sa=None, *, window=150, post_win=100, dark_duty=0.5,
              steady_off=2, off_ratio=0.25):
    """估计受体 trace 在台阶 ``j`` 附近的状态 (:class:`RedState`)。

    亮暗电平由**[j-window, j+post_win) 这一段自己**定 (迟滞二值化), 不用整条 trace
    的 p5/p95: 只短暂亮过一次的分子, 全程分位数全被背景占住, 全局阈值会把噪声帧
    也判成亮态, 于是占空比虚高、灭后比值虚低。

    判类顺序 (先否后正):

    1. ``duty < dark_duty`` → ``'dark'``: 台阶前大半个窗口是暗的, 这是"一次闪光
       的结束", 不是漂白;
    2. ``post_ratio > off_ratio`` → ``'unbleached'``: 灭后没掉到背景, 台阶是暗隙;
    3. ``n_off <= steady_off`` → ``'steady'``: 一直亮着, 供体也该只有噪声;
    4. 其余 → ``'blinking'``: 受体自己在闪, 供体的闪不能算污染。

    窗口内亮暗分不开 (p95-p5 < 4σ) 时返回 ``label='unknown'``, 所有数值为 ``nan``
    —— 没有证据时调用方**不应**据此否决。
    """
    ya = np.asarray(ya, float)
    j = int(j)
    out = RedState()
    if sa is None:
        sa = float(robust_noise(ya))
    sa = max(float(sa), 1e-6)
    a, b = max(0, j - int(window)), min(ya.size, j + int(post_win))
    k = j - a                                   # 台阶在局部窗口里的下标
    if k < 2 or b - k < 2:
        return out
    mask = _on_mask(ya[a:b], sa)
    if mask is None:
        return out
    pre = mask[:k]
    out.duty = float(pre.mean())
    offs = _episodes(~pre, 2)
    out.n_off = int(sum(1 for s, e in offs if e <= k))
    lit = ya[a:a + k][pre]
    lit = lit[np.isfinite(lit)]
    out.level = float(np.median(lit)) if lit.size else np.nan
    post = _med(ya, a + k, b)
    if np.isfinite(post) and np.isfinite(out.level) and out.level > 0:
        out.post_ratio = float(post / out.level)
    if out.duty < float(dark_duty):
        out.label = 'dark'
    elif np.isfinite(out.post_ratio) and out.post_ratio > float(off_ratio):
        out.label = 'unbleached'
    elif out.n_off <= int(steady_off):
        out.label = 'steady'
    else:
        out.label = 'blinking'
    return out


# ---------------------------------------------------------------------------
#  单条 trace
# ---------------------------------------------------------------------------
def find_step(ya, yd, sa=None, sd=None, spec: EventSpec = None):
    """在一条双通道 trace 上找最显著的红降+绿升台阶。

    受体(红)一侧**只把"亮段结束"当候选**, 不在一刀切的滑动窗里找最大值:
    闪烁的受体有多个暗隙, 窗口均值法的 argmax 很容易落在暗隙上, 于是
    ``I_A^pre`` / ``t_bleach`` / ``off_a`` 全跟着错(实测 250/348 的分子因此
    被判成"灭后又亮")。改成亮段结束之后, 真正的漂白会进入候选, 暗隙虽然也在
    候选里, 但会被 ``after`` 判据正确否掉。

    供体(绿)一侧仍按 ``lag`` 搜索, 电平用**亮帧中位** (:func:`post_level`)。

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

    lo = spec.sw + 2
    hi = min(n - spec.sw_long - spec.lag_max - 2, n - spec.sw_post)
    if hi <= lo:
        return None

    on = _on_mask(ya, sa)
    if on is None:
        return None
    episodes = [e for e in _episodes(on, spec.sw) if lo <= e[1] <= hi]
    if not episodes:
        return None

    # 每个亮段的结束都是一次"亮 -> 暗"转变: 候选台阶
    cands = []
    for start, j in episodes:
        seg = ya[max(start, j - spec.sw_long):j]      # 亮段内的帧, 全是 ON
        if seg.size < 2:
            continue
        a_b = float(np.median(seg))
        a_a = _med(ya, j, j + spec.sw)                # 台阶后短窗 (给 off_a 用)
        # 打分用**持续落差**: 闪烁的暗隙在 sw_post 长窗里会被"拉回来",
        # 所以它这一项很小; 真漂白才会给出大的持续落差。
        dA_sus = a_b - _med(ya, j, j + spec.sw_post)
        d_b = _med(yd, j - spec.sw, j)
        thr_d = d_b + 3.0 * sd                        # 供体的亮帧阈值
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
            # 找不到可信亮态: 退回短窗中位数, **只为定位**(判据会否掉这个候选),
            # 这样 find_step 仍然"总能给出一个候选", 和以前的约定一致。
            for L in range(0, spec.lag_max + 1, spec.lag_step):
                v = _med(yd, j + L, j + L + spec.sw) - d_b
                if v > d_best:
                    d_best, l_best, n_best = v, L, 0
        cands.append(dict(score=min(dA_sus / sa, d_best / sd), j=int(j),
                          lag=int(l_best), IA_pre=a_b, IA_post=a_a,
                          ID_pre=d_b, ID_post=d_b + d_best, n_on=int(n_best)))
    if not cands:
        return None
    # "闪"和"漂白"的差别就是**掉到不到背景**: 优先在"短窗后已接近背景"的候选里挑。
    # 全都不满足时退回全部候选 (此时判据会把这个事件否掉, 但 tag 信息留得住)。
    clean = [c for c in cands
             if c['IA_post'] < max(3.0 * sa, spec.off_ratio * c['IA_pre'])]
    best = max(clean or cands, key=lambda c: c['score'])

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
    # 受体状态 (平稳/闪烁/闪光结束/没真灭): 供体侧的判据要按它调整严格度
    ev.red_state = red_state(ya, j, sa, window=spec.pre_duty_win,
                             post_win=spec.sw_post, off_ratio=spec.off_ratio)
    # 受体的逐帧 ON 掩码: pre 窗里的供体电平/离散度**只在受体亮着的帧上**统计 ——
    # 受体熄灭的间隙里供体本来就应该去淬灭变亮, 那不是污染
    ev.red_on = on
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
