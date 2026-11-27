"""供体去淬灭后闪烁时, 事件判据必须仍然认得出来。

真实数据里踩到过: 占空比 0.8 的分子被判成"绿升不持续", 因为固定取
``[j+lag, j+lag+4)`` 这 4 帧的中位数, 恰好落在熄灭的几帧上。
"""
import numpy as np
import pytest

from smlm.fret.calibrate import local_step
from smlm.fret.events import EventSpec, find_step, post_level

F, J, LAG = 400, 200, 2
SD_N = 20.0          # 红通道白噪声
SD_G = 8.0           # 绿通道白噪声


def _pair(*, ia_pre=3000.0, ia_post=0.0, id_pre=250.0, post_on=2500.0,
          blink=(None, None), start_off=True, seed=0):
    """blink=(on 帧数, off 帧数) 在 J+LAG 之后交替; ``(None, None)`` = 常亮。

    ``start_off=True`` 让"灭"的阶段排在最前 —— 这正是真数据里的失败模式:
    台阶后那几帧恰好是灭的。
    """
    rng = np.random.default_rng(seed)
    ia = np.where(np.arange(F) < J, ia_pre, ia_post) + rng.normal(0, SD_N, F)
    idd = np.full(F, id_pre, float)
    t0 = J + LAG
    on_n, off_n = blink
    if on_n is None:
        idd[t0:] = post_on
    else:
        k = np.arange(F - t0) % (on_n + off_n)
        mask = (k >= off_n) if start_off else (k < on_n)
        idd[t0:] = np.where(mask, post_on, 0.0)
    return ia, idd + rng.normal(0, SD_G, F)


def _step(ia, idd, spec=None):
    ev = find_step(ia, idd, None, None, spec or EventSpec(lag_max=10))
    assert ev is not None
    return ev


# ---------------------------------------------------------------------------
#  post_level 本身
# ---------------------------------------------------------------------------
def test_post_level_takes_the_median_of_bright_frames_only():
    y = np.array([0, 0, 100, 110, 0, 120, 0, 90, 0, 0], float)
    lvl, n = post_level(y, 0, y.size, thr=50.0, min_n=3)
    assert n == 4
    assert lvl == pytest.approx(105.0)         # median([90,100,110,120])


def test_post_level_needs_enough_bright_frames():
    y = np.array([0, 0, 500, 0, 0, 0, 0, 0, 0, 0], float)
    lvl, n = post_level(y, 0, y.size, thr=50.0, min_n=4)
    assert n == 1 and np.isnan(lvl)


# ---------------------------------------------------------------------------
#  find_step 在闪烁下仍然合格
# ---------------------------------------------------------------------------
def test_steady_donor_passes_and_reports_the_level():
    ia, idd = _pair()
    ev = _step(ia, idd)
    assert all(ev.flags[k] for k in ('pre_a', 'off_a', 'red', 'grn',
                                     'hold', 'post_d', 'drop_d'))
    assert ev.ID_post == pytest.approx(2500.0, rel=0.1)
    assert ev.dID == pytest.approx(2250.0, rel=0.1)
    assert ev.n_on_d >= 10


def test_blinking_donor_is_still_accepted():
    """台阶后先灭 5 帧、再亮 5 帧 (占空比 50%): 亮态电平不能被灭帧拉低。"""
    ia, idd = _pair(blink=(5, 5), start_off=True)
    ev = _step(ia, idd)
    for k in ('grn', 'hold', 'post_d', 'drop_d'):
        assert ev.flags[k], f'{k} 不该失败'
    # 抗闪烁口径 ≈ 亮态; 纯中位数口径被灭帧拉下来 —— 两者差得多正是闪烁的证据
    assert ev.dID == pytest.approx(2250.0, rel=0.15)
    assert ev.dID_long < 0.5 * ev.dID
    assert ev.n_on_d >= 8


def test_transient_spike_fails_hold_but_not_grn():
    """只亮 4 帧就灭 -> 有"绿升"但没有"持续", hold 负责挡。"""
    ia, idd = _pair(blink=(4, 20), start_off=False)
    ev = _step(ia, idd)
    assert ev.flags['grn'] and ev.flags['post_d']
    assert not ev.flags['hold']
    assert ev.n_on_d < int(np.ceil(0.4 * EventSpec().sw_long))


def test_no_green_rise_has_no_bright_state():
    """供体根本没亮起来: find_step 仍给候选, 但绿光判据全 False。"""
    ia, idd = _pair(post_on=250.0)
    ev = _step(ia, idd)
    assert not (ev.flags['grn'] or ev.flags['hold'] or ev.flags['post_d'])
    assert not np.isfinite(ev.dID)


# ---------------------------------------------------------------------------
#  受体在闪: 台阶必须落在"最后一次亮 -> 暗", 不能落在暗隙上
# ---------------------------------------------------------------------------
def test_blink_gap_in_the_acceptor_is_not_mistaken_for_the_bleach():
    """受体: 亮 30 帧 -> 暗 20 帧(闪) -> 亮 30 帧 -> 灭; 供体在真漂白处才亮。

    老实现 (滑动窗均值取最大红降) 会把 j 放在暗隙的起点上, 于是
    ``I_A^pre`` / ``off_a`` / ``t_bleach`` 全错, 并被 ``after`` 判成"灭后又亮"。
    """
    rng = np.random.default_rng(0)
    lo, hi, blink_a, blink_b = 100, 130, 150, 180       # 亮段1 / 暗隙 / 亮段2
    ia = np.full(F, 5.0)
    ia[lo:hi] = 3000.0
    ia[blink_a:blink_b] = 3000.0
    ia += rng.normal(0, 20.0, F)

    idd = np.full(F, 250.0)
    idd[blink_b:] = 2500.0
    idd += rng.normal(0, 8.0, F)

    ev = _step(ia, idd)
    assert ev.frame == pytest.approx(blink_b, abs=1)     # 真漂白, 不是暗隙起点
    assert ev.IA_pre == pytest.approx(3000.0, rel=0.05)  # 亮段电平, 不是半亮
    assert ev.IA_post < 0.1 * ev.IA_pre
    assert ev.flags['off_a'] and ev.flags['red'] and ev.flags['pre_a']


def test_blink_gap_is_scored_below_the_real_bleach():
    """受体闪 40 帧后复位, 再真灭: 选中的必须是真灭那一刻。

    暗隙的 4 帧后电平也接近 0, 只看短窗落差会被它骗到; 持续落差 (sw_post)
    与"是否掉到背景"两条规则一起才分得开。
    """
    rng = np.random.default_rng(3)
    ia = np.full(F, 5.0)
    ia[40:120] = 2000.0            # 亮
    ia[120:170] = 5.0              # 闪 50 帧
    ia[170:260] = 2000.0           # 又亮
    ia[260:] = 5.0                 # 真灭
    ia += rng.normal(0, 15.0, F)
    idd = np.full(F, 5.0)
    idd[260:] = 1200.0             # 供体只在真灭处亮
    idd += rng.normal(0, 6.0, F)

    ev = find_step(ia, idd, None, None, EventSpec(lag_max=10))
    assert ev is not None
    assert ev.frame == pytest.approx(260, abs=2)
    assert ev.IA_pre > 1500 and ev.IA_post < 100
    assert ev.flags['off_a'] and ev.flags['grn'] and ev.flags['red']


def test_acceptor_without_any_bright_episode_has_no_event():
    """受体全程没有可辨亮态 (只有噪声) -> 连候选都不给。"""
    rng = np.random.default_rng(1)
    ia = 30.0 + rng.normal(0, 20.0, F)
    idd = 250.0 + rng.normal(0, 8.0, F)
    assert find_step(ia, idd, None, None, EventSpec(lag_max=10)) is None


# ---------------------------------------------------------------------------
#  local_step (γ 的分母) 用同一口径
# ---------------------------------------------------------------------------
def test_local_step_is_blink_robust():
    ia, idd = _pair(blink=(5, 5), start_off=True)
    ev = _step(ia, idd)
    loc = local_step(ia, idd, ev.frame, ev.lag, 4, sd=SD_G)
    assert loc is not None
    assert loc['dID'] == pytest.approx(2250.0, rel=0.15)
    assert loc['n_on'] >= 2                # >= max(2, win//2)


def test_local_step_plain_window_would_have_been_low():
    """对照: 老口径 (固定 4 帧中位数) 在闪烁下会明显偏低。"""
    ia, idd = _pair(blink=(5, 5), start_off=True)
    ev = _step(ia, idd)
    j, lag = int(ev.frame), int(ev.lag)
    seg = idd[j + lag: j + lag + 4]
    assert float(np.median(seg)) < 0.3 * ev.ID_post
