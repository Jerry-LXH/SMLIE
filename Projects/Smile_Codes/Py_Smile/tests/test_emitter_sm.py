"""analyze_emitters_sm: 适用于闪烁的固定单分子吗?"""
import inspect

import numpy as np
import pytest

from smlm.datatypes import Emitter, Localizations
from smlm.emitter import (EccentricityFilter, EmitterMerger,
                          analyze_emitters_sm, ecc_pvalue)
from smlm.emitter.base import _TrackBuilder

KW = dict(min_frames=3, first_frame_range=None)


def locs_of(spots, jitter=0.15, seed=0):
    """spots: [((x, y), frames), ...] -> 带高斯定位噪声的 Localizations。"""
    rng = np.random.default_rng(seed)
    frames, xs, ys = [], [], []
    for (x, y), fr in spots:
        fr = list(fr)
        frames += fr
        xs += list(x + rng.normal(0, jitter, len(fr)))
        ys += list(y + rng.normal(0, jitter, len(fr)))
    return Localizations(frame=np.array(frames), x=np.array(xs), y=np.array(ys),
                         photons=np.full(len(xs), 500.0), bg=np.zeros(len(xs)),
                         sigma=np.full(len(xs), 1.5), lp=np.full(len(xs), 0.2))


def emitter_of(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    return Emitter(frame=np.arange(xs.size), x=xs, y=ys,
                   loc_idx=np.arange(xs.size))


# ---------------------------------------------------------------------------
#  连接
# ---------------------------------------------------------------------------
def test_default_motion_is_static():
    """固定分子是这条链的默认假设。"""
    p = inspect.signature(analyze_emitters_sm).parameters['motion']
    assert p.default == 'static'


def test_blinking_static_molecule_stays_one_emitter():
    spot = ((10.0, 10.0), list(range(0, 10)) + list(range(12, 22))
            + list(range(24, 34)))                  # 三次闪烁, 暗态 2 帧
    kept, all_em = analyze_emitters_sm(
        locs_of([spot]), 40, search_radius=1.0, max_gap=10, **KW)

    assert len(all_em) == 1
    assert len(kept) == 1
    assert kept[0].n_locs == 30
    assert kept[0].duty_cycle == pytest.approx(30 / 34)


def test_static_holds_across_a_long_dark_gap():
    spot = ((10.0, 10.0), list(range(0, 10)) + list(range(40, 50)))
    kept, all_em = analyze_emitters_sm(
        locs_of([spot], seed=1), 60, search_radius=1.0, max_gap=60, **KW)
    assert len(all_em) == 1
    assert len(kept) == 1


def test_linear_prediction_drifts_away_for_a_static_molecule():
    """linear 的速度是逐帧差分出来的噪声, 暗态越长外推越偏 —— 默认 static 的原因。"""
    rng = np.random.default_rng(0)
    tb = _TrackBuilder(0, 10.0, 10.0, 0)
    for f in range(1, 6):
        tb.append(f, 10.0 + rng.normal(0, 0.15), 10.0 + rng.normal(0, 0.15), f)

    far = 60
    d_static = np.hypot(*(tb.predict(far, 'static') - [10.0, 10.0]))
    d_linear = np.hypot(*(tb.predict(far, 'linear') - [10.0, 10.0]))
    assert d_static < 0.5
    assert d_linear > 3 * d_static


# ---------------------------------------------------------------------------
#  合并
# ---------------------------------------------------------------------------
def test_merge_reunites_segments_split_by_a_too_long_gap():
    spot = ((10.0, 10.0), list(range(0, 10)) + list(range(40, 50)))
    locs = locs_of([spot], seed=1)

    split = analyze_emitters_sm(locs, 60, search_radius=1.0, max_gap=5, **KW)[1]
    assert len(split) == 2

    kept, all_em = analyze_emitters_sm(
        locs, 60, search_radius=1.0, max_gap=5,
        do_merge=True, r_merge=1.5, **KW)
    assert len(all_em) == 1
    assert len(kept) == 1
    assert kept[0].n_locs == 20


def test_merge_keeps_two_molecules_that_shine_together():
    """真能分辨的近邻分子(同时亮)不该被合并, 也不该被距离阈值删掉。"""
    spots = [((10.0, 10.0), range(20)), ((11.5, 10.0), range(20))]
    kept, all_em = analyze_emitters_sm(
        locs_of(spots, seed=2), 20, search_radius=1.0, max_gap=2,
        do_merge=True, r_merge=1.5, **KW)
    assert len(all_em) == 2
    assert len(kept) == 2


def test_merger_needs_low_temporal_overlap():
    def em(x, frames):
        frames = np.asarray(frames, int)
        return Emitter(frame=frames, x=np.full(frames.size, float(x)),
                       y=np.full(frames.size, 10.0),
                       loc_idx=np.arange(frames.size))

    assert len(EmitterMerger(r_merge=1.5, max_overlap=0.3)(
        [em(10.0, range(0, 5)), em(10.2, range(20, 25))])) == 1
    assert len(EmitterMerger(r_merge=1.5, max_overlap=0.3)(
        [em(10.0, range(10)), em(11.0, range(10))])) == 2


# ---------------------------------------------------------------------------
#  偏心率: 零分布与判据
# ---------------------------------------------------------------------------
def _mc_tail(n, t, reps=20000, seed=0):
    """Monte Carlo: 各向同性 n 点云中 λmax/λmin > t 的比例。"""
    rng = np.random.default_rng(seed)
    p = rng.normal(0, 1, (reps, n, 2))
    p -= p.mean(axis=1, keepdims=True)
    s = np.einsum('nki,nkj->nij', p, p) / (n - 1)
    w = np.linalg.eigvalsh(s)                       # (reps, 2) 升序
    return float((w[:, 1] / w[:, 0] > t).mean())


@pytest.mark.parametrize('n,t', [(10, 6.0), (20, 4.9), (60, 2.2)])
def test_ecc_pvalue_matches_monte_carlo(n, t):
    """闭式零分布必须和模拟的尾概率一致 (这是判据正确性的根据)。"""
    emp = _mc_tail(n, t)
    pred = ecc_pvalue(t, n)
    se = np.sqrt(max(emp * (1 - emp), 1e-12) / 20000)
    assert abs(pred - emp) < 4 * se + 2e-3, (n, t, pred, emp)


def test_ecc_pvalue_bounds_and_monotonicity():
    assert ecc_pvalue(1.0, 20) == pytest.approx(1.0)
    vals = [ecc_pvalue(t, 20) for t in (2, 3, 4, 6, 10)]
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert all(0.0 <= v <= 1.0 for v in vals)
    assert ecc_pvalue(2.0, 2) == 1.0                # n 太小不判定


def test_eccentricity_false_positive_rate_is_p_crit():
    """n=10 的小轨迹: 旧判据 3*(1+2√(2/n)) 等价于 6.5% 的假阳性, 新判据应是 0.1%。"""
    rng = np.random.default_rng(11)
    ems = [emitter_of(10 + rng.normal(0, 0.15, 10),
                      10 + rng.normal(0, 0.15, 10)) for _ in range(1000)]

    kept = EccentricityFilter(p_crit=1e-3).keep(ems)
    assert kept.mean() >= 0.99                       # 1000 条里最多误杀 10 条

    # 旧口径在同一批数据上的误杀率明显更高 (n=10 时约 6.5%)
    legacy = np.mean([ecc_pvalue(e.eccentricity(), e.n_locs) < 0.065
                      for e in ems])
    assert legacy > 0.03


def test_eccentricity_catches_elongation_along_an_axis():
    """沿 x 排开的两个分子: 归一化会抹掉这种拉长, 协方差比值不会。"""
    rng = np.random.default_rng(12)
    n = 60
    em = emitter_of(10.0 + rng.normal(0, 0.50, n),      # 沿 x 拉长, 与 y 不相关
                    10.0 + rng.normal(0, 0.10, n))
    assert em.eccentricity() > 10
    assert not EccentricityFilter().keep([em])[0]
    # 各向同性对照: 同样的 n 与定位精度
    iso = emitter_of(10.0 + rng.normal(0, 0.30, n),
                     10.0 + rng.normal(0, 0.30, n))
    assert EccentricityFilter().keep([iso])[0]


def test_eccentricity_drops_a_two_molecule_centroid_track():
    """两个分子被拟合成一个质心 (亮度缓慢起伏) -> 点云沿连线拉长 -> 被剔除。"""
    rng = np.random.default_rng(3)
    n = 60
    frac = 0.5 + 0.4 * np.sin(np.arange(n) * 2 * np.pi / 40)
    em = emitter_of(10.0 + 1.5 * frac + rng.normal(0, 0.10, n),
                    10.0 + rng.normal(0, 0.10, n))
    assert not EccentricityFilter().keep([em])[0]


def test_isotropic_single_molecule_survives_the_whole_chain():
    rng = np.random.default_rng(4)
    n = 40
    locs = Localizations(frame=np.arange(n),
                         x=10.0 + rng.normal(0, 0.10, n),
                         y=10.0 + rng.normal(0, 0.10, n),
                         photons=np.full(n, 500.0), bg=np.zeros(n),
                         sigma=np.full(n, 1.5), lp=np.full(n, 0.2))
    kept, all_em = analyze_emitters_sm(locs, n, search_radius=1.0,
                                       max_gap=5, **KW)
    assert len(all_em) == 1
    assert len(kept) == 1


def test_long_traces_with_tiny_drift_survive_the_effect_size_floor():
    """长轨迹 + 0.3 px 残余漂移: 纯 p 值大批误杀, 效应量下限能救回来。

    n=828 时纯 p 值只要轴比 1.30 就判假, 而 0.3 px 漂移 (远小于分辨率)
    就能把轴比拉到 1.35 —— 绿通道每个 emitter 有 800+ 个定位点, 正是这种情况。
    """
    rng = np.random.default_rng(0)
    n, trials = 828, 60
    t = np.linspace(-0.5, 0.5, n)
    ems = [emitter_of(10.0 + 0.3 * t + rng.normal(0, 0.15, n),
                      10.0 + rng.normal(0, 0.15, n))
           for _ in range(trials)]

    kept_new = EccentricityFilter().keep(ems).mean()
    kept_old = EccentricityFilter(ecc_min=None).keep(ems).mean()

    assert kept_old < 0.45          # 纯 p 值: 大部分真轨迹被判成"两个分子"
    assert kept_new > 0.85          # 加效应量下限后绝大多数保留


def test_effect_size_floor_still_drops_a_real_double():
    """两个分子间距 0.5 px: 效应量足够, 加了下限也要剔除。"""
    rng = np.random.default_rng(22)
    n = 828
    w = rng.random(n)                      # 相对亮度起伏 -> 质心在两点间走
    em = emitter_of(10.0 + 0.5 * w + rng.normal(0, 0.15, n),
                    10.0 + rng.normal(0, 0.15, n))

    assert em.eccentricity() > 1.5
    assert not EccentricityFilter().keep([em])[0]


def test_effect_size_floor_protects_short_noisy_traces():
    """短轨迹的纯噪声轴比本来就能到 1.7+, 只要 p 不显著就不该剔除。"""
    rng = np.random.default_rng(23)
    ems = [emitter_of(10.0 + rng.normal(0, 0.15, 20),
                      10.0 + rng.normal(0, 0.15, 20)) for _ in range(200)]
    big = [e for e in ems if e.eccentricity() > 1.5]
    assert big                                  # 确实有一批轴比超过 1.5
    assert all(EccentricityFilter().keep([e])[0] for e in big)
