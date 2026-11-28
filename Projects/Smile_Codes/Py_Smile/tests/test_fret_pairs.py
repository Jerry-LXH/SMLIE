"""挑选合格 FRET trace 对的判据: 每条门限各一个合成样例。

合成 trace 的噪声用固定 σ, 不追求与 Poisson 模型完全一致 —— 判据本身是
"观测离散度 / 期望噪声" 的比值, 只要比值落在门限两侧就足以验证行为。
"""
import numpy as np
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from smlm.datatypes import Emitter, Trace                     # noqa: E402
from smlm.fret import pairs as P                              # noqa: E402
from smlm.fret.events import EventSpec                        # noqa: E402

F = 400
J = 200                     # 红光漂白帧
LAG = 2                     # 绿光去淬灭滞后
SW, SW2, AREA = 5.0, 8.0, 50.0
BG_STD = 2.0
RADIUS = 4.0                # -> 归属半径 0.8*4 = 3.2 px


def _trace(intensity, *, x=10.0, y=10.0, neighbours=0, radius=RADIUS):
    v = np.asarray(intensity, float)
    n = v.size
    return Trace(emitter_id=0, frames=np.arange(n), x=float(x), y=float(y),
                 channel='t', intensity=v, raw=v.copy(),
                 bg_flux=np.zeros(n), bg=np.zeros(n),
                 bg_std=np.full(n, BG_STD), sw=SW, sw2=SW2, area=AREA,
                 radius=float(radius), neighbours=int(neighbours),
                 sigma=1.2)


def _pair(*, ia_pre=3000.0, ia_post=0.0, id_pre=250.0, id_post=2500.0,
          seed=0, sd_a=20.0, sd_d=8.0, n=F, j=J, lag=LAG):
    rng = np.random.default_rng(seed)
    ia = np.full(n, ia_pre)
    ia[j:] = ia_post
    idd = np.full(n, id_pre)
    if lag is not None:
        idd[j + lag:] = id_post
    return (_trace(ia + rng.normal(0, sd_a, n)),
            _trace(idd + rng.normal(0, sd_d, n)))


def _donor(x=10.0, y=10.0, frames=None, n_photons=2000.0):
    fr = np.arange(F) if frames is None else np.asarray(frames, int)
    return Emitter(frame=fr, x=np.full(fr.size, float(x)),
                   y=np.full(fr.size, float(y)), loc_idx=fr.copy())


ON_N, OFF_N = 30, 10        # 受体: 亮 30 帧 / 灭 10 帧
UNQ = 2500.0                # 供体未被淬灭时的亮度


def _anti_correlated(*, id_pre=250.0, sd_a=20.0, sd_d=8.0):
    """受体在闪, 供体在受体熄灭的间隙里**去淬灭变亮** —— 物理行为, 不是污染。

    受体亮着时供体被淬灭到 ``id_pre``, 受体灭时升到 ``UNQ``; 台阶之后供体常亮。
    """
    rng = np.random.default_rng(0)
    red_on = (np.arange(F) % (ON_N + OFF_N)) < ON_N
    ia = np.where(red_on, 3000.0, 0.0)
    ia[J:] = 0.0
    idd = np.where(red_on, id_pre, UNQ)
    idd[J + LAG:] = UNQ
    return (_trace(ia + rng.normal(0, sd_a, F)),
            _trace(idd + rng.normal(0, sd_d, F))), red_on


def _run(pa, pd, emitters=None, **kw):
    crit = P.PairCriteria(**kw)
    recs, stats = P.select_pairs([pa], [pd], emitters_d=emitters,
                                 crit=crit, verbose=False)
    return recs[0], stats


# ---------------------------------------------------------------------------
#  通过
# ---------------------------------------------------------------------------
def test_clean_pair_is_kept():
    pa, pd = _pair()
    rec, stats = _run(pa, pd, [_donor()])
    assert rec.keep, rec.reasons
    assert rec.reasons == []
    assert rec.n_d_pre == 1          # 供体自己
    assert rec.spread < 2.5
    # find_step 的 j 是"台阶边界帧", ±1 会落在真实帧上; lag 只要不超上限
    assert abs(rec.j - J) <= 3
    assert 0 <= rec.event.lag <= 10
    assert stats['n_keep'] == 1


def test_split_fragments_of_one_donor_count_as_one():
    """定位器把同一分子拆成两条 (相距 < attr_dup) 时不该被误杀。"""
    pa, pd = _pair()
    em = [_donor(), _donor(x=10.4, y=10.3)]      # 0.5 px 内
    rec, _ = _run(pa, pd, em)
    assert rec.n_in_ap == 2
    assert rec.n_d_pre == 1
    assert rec.keep, rec.reasons


# ---------------------------------------------------------------------------
#  各级否决
# ---------------------------------------------------------------------------
def test_dark_before_bleach_is_rejected():
    """漂白前 10 帧亮度 ≤150 -> 'red'。"""
    pa, pd = _pair(ia_pre=40.0, ia_post=0.0, sd_a=5.0, sd_d=8.0)
    rec, stats = _run(pa, pd, [_donor()])
    assert not rec.keep
    assert rec.reasons[0] == 'red'
    assert rec.IA_pre_win <= 150
    assert stats['funnel'][1][2] == 0        # 第 2 级 (red) 就全砍


def test_red_gate_uses_pre_bleach_window_not_the_global_max():
    """一个热像素也能让 max(intensity) 很高, 但漂白前其实没在发光。

    旧口径 ``max(intensity) < red_min`` 会把这个放进来; 新口径看台阶前 10 帧,
    会把它挡掉。
    """
    rng = np.random.default_rng(4)
    ia = np.full(F, 40.0) + rng.normal(0, 5.0, F)
    ia[50] = 3000.0                            # 单帧热像素/爆闪
    ia[J:] = 0.0                               # 200 帧之后的"漂白"
    pa = _trace(ia + rng.normal(0, 5.0, F))
    pd = _pair()[1]

    rec, _ = _run(pa, pd, [_donor()])
    assert rec.IA_max > 150                    # 全程最大值很亮
    assert rec.IA_pre_win <= 150               # 但漂白前是暗的
    assert not rec.keep and 'red' in rec.reasons


def test_pre_bleach_level_is_a_median_of_the_window():
    ia = np.zeros(F)
    ia[J - 10:J] = [10.0, 20.0, 30.0, 40.0, 50.0,
                    60.0, 70.0, 80.0, 90.0, 1000.0]      # 一帧是尖峰
    # 高于窗内中位的那些帧的中位: 60..1000 -> 80 (干净 trace 与普通中位只差一点)
    assert P.pre_bleach_level(_trace(ia), J, 10) == pytest.approx(80.0)
    assert np.isnan(P.pre_bleach_level(_trace(ia), 0, 10))


def test_pre_bleach_level_ignores_dark_gaps_from_blinking():
    """受体在漂白前闪了几下: 普通中位会被暗帧拉低, 亮帧中位不会。"""
    ia = np.zeros(F)
    ia[J - 10:J] = [0.0, 300.0, 0.0, 300.0, 0.0, 300.0, 0.0, 300.0, 0.0, 300.0]
    assert P.pre_bleach_level(_trace(ia), J, 10) == pytest.approx(300.0)
    assert float(np.median(ia[J - 10:J])) == pytest.approx(150.0)


def test_no_red_drop_is_rejected():
    """红光不漂白 -> 受体没有"亮段结束", 连候选都没有。"""
    pa, pd = _pair(ia_pre=3000.0, ia_post=3000.0)
    rec, _ = _run(pa, pd, [_donor()])
    assert not rec.keep
    assert rec.event is None
    assert rec.reasons == ['event']


def test_too_short_trace_has_no_event():
    """帧数不够做台阶搜索时 event 为 None, 落在 'event' 级。"""
    pa, pd = _pair(n=20, j=10)
    rec, _ = _run(pa, pd, [_donor()])
    assert rec.event is None
    assert rec.reasons == ['event']


def test_blinking_neighbour_window_is_rejected():
    """pre 窗里绿光忽高忽低 = 邻居在闪, 不是被淬灭的供体。"""
    pa, pd = _pair()
    idd = np.array(pd['intensity'])
    blink = (np.arange(F) // 8) % 2 == 0
    idd[:J] += blink[:J] * 1600.0            # 暗/亮交替
    rec, _ = _run(pa, _trace(idd), [_donor()])
    assert not rec.keep
    assert 'spread' in rec.reasons
    assert rec.spread > 2.5


def test_two_green_molecules_in_aperture_are_rejected():
    """pre 窗内孔径里有第二个绿光分子 -> 绿光亮度不可信。"""
    pa, pd = _pair()
    em = [_donor(), _donor(x=11.5, y=10.0)]       # 1.5 px, 超过去重半径
    rec, _ = _run(pa, pd, em)
    assert rec.n_d_pre == 2
    assert not rec.keep and 'attr' in rec.reasons


def test_emitter_only_after_bleach_does_not_count():
    """只在漂白之后才出现的绿光分子不算 pre 窗污染。"""
    pa, pd = _pair()
    late = np.arange(J + LAG, F)
    em = [_donor(), _donor(x=11.5, y=10.0, frames=late)]
    rec, _ = _run(pa, pd, em)
    assert rec.n_in_ap == 2 and rec.n_d_pre == 1
    assert rec.keep, rec.reasons


def test_red_recovering_after_bleach_is_rejected():
    """孔径里还有第二个受体: 红光灭了又亮 -> I_A^pre / ΔI_A 不可信。"""
    pa, pd = _pair()
    ia = np.array(pa['intensity'])
    ia[J + 100:] += 3000.0
    rec, _ = _run(_trace(ia), pd, [_donor()])
    assert not rec.keep and 'after' in rec.reasons
    assert rec.red_after_duty > 0.2


def test_after_gate_tolerates_a_short_blink_after_the_bleach():
    """灭后只是短暂闪一下 (占空比 10%) -> 不算"又亮", 该合格。

    这条是给会闪的样品用的: 用 p95 判"有没有再亮过"会把它们全砍掉。
    """
    pa, pd = _pair()
    ia = np.array(pa['intensity'])
    ia[J + 50:J + 70] += 2500.0          # 200 帧里亮 20 帧
    rec, _ = _run(_trace(ia), pd, [_donor()])
    assert rec.red_after_duty == pytest.approx(0.10, abs=0.02)
    assert 'after' not in rec.reasons
    assert rec.keep, rec.reasons


def test_after_gate_switches_back_to_p95_when_duty_is_disabled():
    pa, pd = _pair()
    ia = np.array(pa['intensity'])
    ia[J + 50:J + 70] += 2500.0
    rec, _ = _run(_trace(ia), pd, [_donor()], after_duty=None)
    assert 'after' in rec.reasons       # 老口径: p95 摸到那 20 帧就否


def test_flash_ending_is_rejected_by_the_duty_gate():
    """受体在"漂白"之前大部分时间是暗的 -> 这是一次闪光的结束, 不是漂白。"""
    rng = np.random.default_rng(0)
    ia = np.full(F, 5.0)
    ia[J - 60:J] = 3000.0                    # 只亮 60 帧
    pa = _trace(ia + rng.normal(0, 20.0, F))
    pd = _pair()[1]

    rec, _ = _run(pa, pd, [_donor()])
    assert rec.red_label == 'dark'
    assert rec.red_duty < 0.5
    assert 'duty' in rec.reasons and not rec.keep
    # 关掉这条判据时事件本身还在, 说明拦它的是占空比而不是别的
    off, _ = _run(pa, pd, [_donor()], red_duty_min=None)
    assert off.event is not None and 'duty' not in off.reasons


def test_pre_spread_ignores_the_donor_lighting_up_in_acceptor_gaps():
    """受体闪时供体跟着闪是物理行为; 只在受体亮着的帧上量离散度就不会误杀。"""
    (pa, pd), red_on = _anti_correlated()
    rec, _ = _run(pa, pd, [_donor()])
    assert rec.red_label == 'blinking'
    assert rec.spread < 2.5, f'spread={rec.spread}'
    assert rec.n_pre_gate == pytest.approx(60 * ON_N / (ON_N + OFF_N), abs=8)
    assert rec.keep, rec.reasons
    # 对照: 不做门控时, 这段"受体灭-供体亮"的交替会被当成邻居污染
    assert P.pre_spread(pd, rec.j, 60)[0] > 5.0


def test_efficiency_measures_id_pre_only_while_the_acceptor_is_on():
    """I_D^pre 混进"受体灭、供体已去淬灭"的帧, E 会掉一半。"""
    import smlm.fret as fr

    (pa, pd), _ = _anti_correlated()
    rec, _ = _run(pa, pd, [_donor()])
    cal = fr.FretCalibration(gamma=1.1)

    gated, _ = fr.fret_efficiency([rec.event], [pa], [pd], cal, verbose=False)
    plain, _ = fr.fret_efficiency([rec.event], [pa], [pd], cal,
                                  fr.EfficiencySpec(gate_pre_d=False),
                                  verbose=False)
    assert gated[0].ID_pre == pytest.approx(250.0, abs=60)
    assert gated[0].E == pytest.approx(0.92, abs=0.04)
    assert plain[0].ID_pre > 2.5 * gated[0].ID_pre      # 未转移态被平均进来
    assert plain[0].E < gated[0].E - 0.10


def test_neighbour_gate_is_optional():
    """neighbours 取的是**供体**通道的环内邻居数 (污染 I_D^pre 的那一侧)。"""
    pa, pd = _pair()
    pd = _trace(pd['intensity'], neighbours=45)
    assert _run(pa, pd, [_donor()])[0].keep          # 默认不限
    rec, _ = _run(pa, pd, [_donor()], neighbours_max=10)
    assert not rec.keep and 'neigh' in rec.reasons


def test_without_emitters_the_attribution_stage_is_skipped():
    pa, pd = _pair()
    rec, stats = _run(pa, pd, None)
    assert rec.keep
    assert 'attr' not in [k for k, _, _ in stats['funnel']]


# ---------------------------------------------------------------------------
#  漏斗 / 批量
# ---------------------------------------------------------------------------
def test_funnel_is_monotonic_and_matches_keep():
    batch = [_pair()[0:2], _pair(ia_pre=30.0)[0:2], _pair(ia_post=3000.0)[0:2]]
    ta = [p[0] for p in batch]
    td = [p[1] for p in batch]
    recs, stats = P.select_pairs(ta, td, emitters_d=[_donor()], verbose=False)
    counts = [n for _, _, n in stats['funnel']]
    assert counts == sorted(counts, reverse=True)
    assert counts[-1] == stats['n_keep'] == sum(r.keep for r in recs)
    for r in recs:
        assert r.keep == (not r.reasons)
    assert len(recs) == len(ta)


def test_spec_lag_max_is_overridden_by_the_criteria():
    """``spec`` 只覆盖窗宽/门限, ``lag_max`` 永远听 ``crit`` 的。

    两个旋钮指同一件事时, 曾经是"传了 spec 就静默丢掉 crit.lag_max",
    真数据上表现为漏斗数字对不上、却没人报错。
    """
    pa, pd = _pair()
    crit = P.PairCriteria(lag_max=LAG)
    narrow, _ = P.select_pairs([pa], [pd], crit=crit,
                               spec=EventSpec(lag_max=0), verbose=False)
    wide, _ = P.select_pairs([pa], [pd], crit=crit,
                             spec=EventSpec(lag_max=16), verbose=False)
    assert narrow[0].event is not None
    assert narrow[0].event.lag == wide[0].event.lag == LAG


def test_trace_count_mismatch_raises():
    with pytest.raises(ValueError):
        P.select_pairs([_pair()[0]], [], verbose=False)


def test_events_carry_their_emitter_id_for_efficiency():
    """回归: ``find_step`` 返回的 emitter_id 是 -1, select_pairs 必须补上。

    漏了的话 ``fret_efficiency`` 会按 -1 去读 ``traces[-1]``, 结果全部静默丢掉
    (实测真数据里 56 条一条都没算出 E, γ 也估不出来)。
    """
    import smlm.fret as fr

    batch = [_pair(seed=s) for s in range(3)]
    ta = [p[0] for p in batch]
    td = [p[1] for p in batch]
    em = []
    for k, tr in enumerate(td):
        tr['x'] = tr['y'] = 10.0 + 6.0 * k        # 各自独立的孔径位置
        em.append(_donor(x=10.0 + 6.0 * k))

    recs, _ = P.select_pairs(ta, td, emitters_d=em, verbose=False)
    events = [r.event for r in recs if r.keep]
    assert len(events) == 3
    assert [e.emitter_id for e in events] == [0, 1, 2]

    cal = fr.FretCalibration(gamma=1.1)
    res, stats = fr.fret_efficiency(events, ta, td, cal, verbose=False)
    assert len(res) == len(events)
    assert stats['n_drop'] == 0 and stats['n_noon'] == 0
    assert all(0.0 < r.E < 1.0 for r in res)


def test_efficiency_refuses_negative_emitter_id():
    """宁可响亮报错, 也不要用负下标静默读错 trace。"""
    import smlm.fret as fr

    pa, pd = _pair()
    ev = fr.find_step(pa['intensity'], pd['intensity'])
    assert ev.emitter_id == -1
    with pytest.raises(ValueError, match='emitter_id'):
        fr.fret_efficiency([ev], [pa], [pd], fr.FretCalibration(gamma=1.1),
                           verbose=False)


def test_no_green_emitter_in_aperture_reports_zero_not_unknown():
    """孔径里一个绿光中心都没有时 n_d_pre 是 0 (确定没有), 不是 -1 (未知)。"""
    pa, pd = _pair()
    pd['x'], pd['y'] = 10.0, 10.0
    rec, _ = _run(pa, pd, [_donor(x=40.0, y=40.0)])     # 绿光分子远在别处
    assert rec.n_in_ap == 0
    assert rec.n_d_pre == 0
    assert rec.keep, rec.reasons


def test_min_pre_gate_rejects_short_pre_window():
    """提前漂白的分子 pre 窗很短 (这里 j=20), min_pre 要能把它们挡掉。"""
    import smlm.fret as fr

    pa, pd = _pair(n=400, j=20)
    ev = P.select_pairs([pa], [pd], emitters_d=[_donor()], verbose=False)[0][0].event
    cal = fr.FretCalibration(gamma=1.1)

    loose = fr.EfficiencySpec()
    res, st = fr.fret_efficiency([ev], [pa], [pd], cal, loose, verbose=False)
    assert len(res) == 1 and st['n_pre'] == 0

    tight = fr.EfficiencySpec(min_pre=50)
    res2, st2 = fr.fret_efficiency([ev], [pa], [pd], cal, tight, verbose=False)
    assert res2 == [] and st2['n_pre'] == 1


def test_flags_rejection_records_which_sub_criterion_failed(capsys):
    """'物理门限不合格' 必须能拆到具体哪一条, 否则没法判断。"""
    from smlm.fret import CRITERIA

    # 受体只掉一半 (3000->1500): off_a 过不了, 其余都过
    pa, pd = _pair(ia_post=1500.0)
    crit = P.PairCriteria()
    recs, stats = P.select_pairs([pa], [pd], emitters_d=[_donor()],
                                 crit=crit, verbose=True)
    out = capsys.readouterr().out
    rec = recs[0]
    assert not rec.keep and 'flags' in rec.reasons
    assert rec.failed_flags == ['off_a']
    assert set(rec.failed_flags) <= set(CRITERIA)
    assert '物理门限明细' in out and 'off_a' in out


def test_kept_record_has_no_failed_flags():
    pa, pd = _pair()
    rec, _ = _run(pa, pd, [_donor()])
    assert rec.keep and rec.failed_flags == []


# ---------------------------------------------------------------------------
#  辅助量
# ---------------------------------------------------------------------------
def test_noise_per_frame_follows_poisson_budget():
    tr = _trace(np.full(10, 1000.0))
    # σ = sqrt(I*sw2 + (sw²/sw2)*bg_std²)/sw2
    want = np.sqrt(1000.0 * SW2 + (SW ** 2 / SW2) * BG_STD ** 2) / SW2
    np.testing.assert_allclose(P.noise_per_frame(tr), want)


def test_noise_per_frame_is_nan_without_sw2():
    tr = _trace(np.full(10, 100.0))
    tr['sw2'] = 0.0
    assert np.all(np.isnan(P.noise_per_frame(tr)))


def test_pre_spread_of_pure_noise_is_about_one():
    """噪声正好等于光子预算给出的 σ 时, 归一化离散度 ≈ 1。"""
    pd = _trace(np.full(F, 250.0))
    s = float(np.nanmedian(P.noise_per_frame(pd)))
    pd = _trace(250.0 + np.random.default_rng(3).normal(0, s, F))
    spread, level, ratio = P.pre_spread(pd, J, 60)
    assert spread == pytest.approx(1.0, abs=0.25)
    assert level > 0
    assert ratio == pytest.approx(1.0, abs=0.5)


def test_spread_gate_survives_a_wrong_photon_budget():
    """光子预算高估噪声时, 参照必须改用**实测**噪声, 否则门限形同虚设。

    真数据里就踩到了这个坑: ``spread`` 报 0.11σ —— 公式 σ 比实测大约 9 倍,
    ``k_spread=2.5`` 一条都筛不掉。这里直接对比两种参照口径。
    """
    pa, _ = _pair()
    rng = np.random.default_rng(11)
    idd = 250.0 + rng.normal(0, 8.0, F)               # 真实噪声 σ = 8
    blink = (np.arange(F) // 8) % 2 == 0
    idd[:J] += blink[:J] * 200.0                      # pre 窗里有弱邻居在闪
    pd = _trace(idd)
    pd['bg_std'] = np.full(F, 700.0)                  # 光预算被高估 ~20 倍

    spread, level, ratio = P.pre_spread(pd, J, 60)
    s_budget = float(np.nanmedian(P.noise_per_frame(pd)))
    seg = pd['intensity'][J - 60:J]
    old = (np.percentile(seg, 84) - np.percentile(seg, 16)) / (1.989 * s_budget)

    assert ratio > 5.0 and s_budget > 5 * 8.0         # 预算确实高估
    assert old < 2.5                                  # 旧口径 (按预算归一) 抓不住
    assert spread > 2.5                               # 新口径 (按实测归一) 抓得住
    rec, _ = _run(pa, pd, [_donor()])
    assert not rec.keep and 'spread' in rec.reasons


def test_pre_spread_needs_enough_frames():
    pa, pd = _pair()
    assert np.isnan(P.pre_spread(pd, 3, 60)[0])


def test_attach_efficiency_matches_by_emitter_id():
    class _R:
        def __init__(self, i, e):
            self.emitter_id, self.E = i, e

    recs = [P.PairRecord(emitter_id=0), P.PairRecord(emitter_id=1),
            P.PairRecord(emitter_id=2)]
    P.attach_efficiency(recs, [_R(0, 0.71), _R(2, 0.33)])
    assert recs[0].E == pytest.approx(0.71)
    assert np.isnan(recs[1].E)
    assert recs[2].E == pytest.approx(0.33)


# ---------------------------------------------------------------------------
#  画廊
# ---------------------------------------------------------------------------
def test_efficiency_plot_hides_e_don_by_default():
    """默认只画 E；E_don 要显式 show_don=True 才叠上去（它常被当成第二个分布）。"""
    import matplotlib.pyplot as plt
    import smlm.fret as fr
    from smlm.fret import viz

    res = [fr.FretResult(emitter_id=i, frame=0, E=0.5 + 0.05 * i,
                         E_don=0.4 + 0.05 * i) for i in range(5)]

    fig, ax = plt.subplots()
    viz.plot_efficiency(res, ax=ax)
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert len(labels) == 1 and 'don' not in labels[0].lower()

    fig2, ax2 = plt.subplots()
    viz.plot_efficiency(res, ax=ax2, show_don=True)
    labels2 = [t.get_text() for t in ax2.get_legend().get_texts()]
    assert any('don' in v.lower() for v in labels2)


def test_pair_gallery_renders_both_groups():
    from matplotlib.figure import Figure
    from smlm.fret import viz

    batch = [_pair()[0:2], _pair(ia_pre=30.0)[0:2], _pair(ia_post=3000.0)[0:2]]
    recs, _ = P.select_pairs([p[0] for p in batch], [p[1] for p in batch],
                             emitters_d=[_donor()], verbose=False)
    fig = viz.plot_pair_gallery(recs, [p[0] for p in batch],
                                [p[1] for p in batch], n=3, ex_time=0.1)
    assert isinstance(fig, Figure)
    titles = [ax.get_title() for ax in fig.axes]
    assert any('KEEP' in t for t in titles)
    assert any('REJECT' in t for t in titles)


def test_pair_gallery_wraps_into_rows():
    """条数超过 ncols 时折行, 面板总数要覆盖全部样例。"""
    from matplotlib.figure import Figure
    from smlm.fret import viz

    batch = [_pair(seed=s) for s in range(5)]
    for k, (pa, pd) in enumerate(batch):
        pa['x'] = pd['x'] = pa['y'] = pd['y'] = 10.0 + 6.0 * k
    recs, _ = P.select_pairs([p[0] for p in batch], [p[1] for p in batch],
                             emitters_d=[_donor(x=10.0 + 6.0 * k)
                                         for k in range(5)], verbose=False)
    fig = viz.plot_pair_gallery(recs, [p[0] for p in batch],
                                [p[1] for p in batch], which='kept',
                                n=5, ncols=2, ex_time=0.1)
    assert isinstance(fig, Figure)
    n_titles = sum(1 for ax in fig.axes if ax.get_title())
    assert n_titles == 5
    assert fig.axes[0].get_subplotspec().get_gridspec().nrows == 3   # ceil(5/2)
