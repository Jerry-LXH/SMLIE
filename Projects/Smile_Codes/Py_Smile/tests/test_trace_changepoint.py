"""变点法状态判定: 分段准确性与和阈值法的一致性。"""
import numpy as np
import pytest

import smlm.trace as tr


def _intensity(blocks, n=400, on=100.0, off=5.0, noise=1.0, seed=0):
    """blocks: [(start, stop), ...] 里亮, 其余暗; 可叠加高斯噪声。"""
    y = np.full(n, float(off))
    for start, stop in blocks:
        y[start:stop] = on
    if noise:
        y = y + np.random.default_rng(seed).normal(0.0, noise, n)
    return y


def _trace(intensity, emitter_id=0):
    return {'emitter_id': emitter_id,
            'frames': np.arange(len(intensity)),
            'x': 1.0, 'y': 1.0,
            'intensity': np.asarray(intensity, float)}


def _pipe(n_frames):
    return tr.TraceFillPipeline(np.zeros((n_frames, 4, 4), np.float32),
                                np.array([[2.0, 2.0]]), sigma=1.5,
                                n_frames=n_frames)


def _fill(trace, filler):
    _pipe(len(trace['intensity'])).add_fields([trace], [filler])
    return trace


BLINKS = [(20, 120), (150, 260), (300, 340)]


def test_recovers_known_blinks_and_bleach():
    y = _intensity(BLINKS)
    trace = _fill(_trace(y), tr.ChangepointStates(bleach_frames=50))

    truth = np.zeros(y.size, bool)
    for start, stop in BLINKS:
        truth[start:stop] = True

    assert trace['bleach_frame'] == 340          # 末尾 60 帧 OFF >= 50
    assert trace['n_on'] == 3
    assert trace['n_switches'] == 5
    assert trace['n_states'] == 2
    # 漂白之前的 ON/OFF 逐帧判对率
    agree = ((trace['state'][:340] > tr.OFF) == truth[:340]).mean()
    assert agree > 0.98
    assert not trace['alive'][340:].any()
    assert trace['noise'] == pytest.approx(1.0, rel=0.2)


def test_agrees_with_threshold_on_clean_data():
    """低噪声下两条路线必须收敛到同一组闪烁计数。"""
    y = _intensity(BLINKS, noise=1.0)
    cp = _fill(_trace(y), tr.ChangepointStates(bleach_frames=50))
    th = _fill(_trace(y), tr.ThresholdStates(high=50.0, bleach_frames=50))

    for key in ('n_on', 'n_switches', 'n_cp', 'bleach_frame', 'on_frames'):
        assert cp[key] == th[key], key
    np.testing.assert_array_equal(cp['state'], th['state'])


def test_single_frame_glitch_is_not_a_segment():
    y = _intensity([(20, 120)], n=300)
    y[200] = 500.0                                # 孤立的 1 帧亮

    trace = _fill(_trace(y), tr.ChangepointStates(bleach_frames=500, min_seg=5))

    assert trace['n_on'] == 1
    assert trace['n_states'] == 2                 # 没有多出一个亮度等级
    assert int(np.asarray(trace['state']).max()) == 1


def test_three_levels_are_kept_as_three_states_and_binarizable():
    y = np.empty(400)
    y[:120] = 5.0
    y[120:260] = 50.0
    y[260:] = 100.0
    y = y + np.random.default_rng(3).normal(0.0, 1.0, y.size)

    trace = _fill(_trace(y), tr.ChangepointStates(bleach_frames=500))
    state = np.asarray(trace['state'])
    fit = np.asarray(trace['fit'])

    assert trace['n_states'] == 3
    assert sorted(np.unique(state).tolist()) == [0, 1, 2]
    # 等级按亮度升序: 0 = 最暗
    assert fit[10] < fit[150] < fit[300]
    # n_cp 数等级变化 (2 次); n_switches 只数二值 ON/OFF 翻转 (只有 1 次 OFF→ON)
    assert trace['n_cp'] == 2
    assert trace['n_switches'] == 1


def test_flat_and_dark_traces():
    flat = _fill(_trace(_intensity([(0, 400)], noise=1.0)),
                 tr.ChangepointStates(bleach_frames=50))
    assert flat['n_states'] == 1
    assert flat['bleach_frame'] == -1             # 一直亮 -> 没漂白
    assert flat['n_switches'] == 0

    dark = _fill(_trace(np.zeros(300)),
                 tr.ChangepointStates(bleach_frames=5))
    assert dark['n_on'] == 0
    assert (np.asarray(dark['state']) == tr.OFF).all()
    assert dark['bleach_frame'] == -1             # 从没亮过, 谈不上漂白


def test_vst_is_optional_and_validated():
    y = _intensity(BLINKS)
    both = _fill(_trace(y), tr.ChangepointStates(bleach_frames=50,
                                                 vst='anscombe'))
    assert both['n_on'] == 3

    with pytest.raises(ValueError, match='vst 只支持'):
        tr.ChangepointStates(vst='nope')
    with pytest.raises(ValueError, match='min_seg 至少为 1'):
        tr.ChangepointStates(min_seg=0)


def test_segment_returns_contiguous_bounds():
    y = _intensity(BLINKS)
    levels, bounds = tr.changepoint.segment(y)

    assert bounds[0] == 0 and bounds[-1] == y.size
    assert np.all(np.diff(bounds) > 0)
    assert levels.size == y.size
    # 每段内部等级恒定
    for a, b in zip(bounds[:-1], bounds[1:]):
        assert np.unique(levels[a:b]).size == 1
