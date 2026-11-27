"""阈值法状态判定与 trace 筛选。"""
import numpy as np
import pytest

import smlm.trace as tr


def _intensity(blocks, n=60, on=100.0, off=5.0):
    """blocks: [(start, stop), ...] 里亮, 其余暗。"""
    y = np.full(n, float(off))
    for start, stop in blocks:
        y[start:stop] = on
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


def _states(trace, **kw):
    _pipe(len(trace['intensity'])).add_fields(
        [trace], [tr.ThresholdStates(high=kw.pop('high', 50.0), **kw)])
    return trace


def test_three_blinks_then_bleach():
    y = _intensity([(5, 15), (20, 30), (35, 40)], n=200)
    trace = _states(_trace(y), bleach_frames=100)

    state = trace['state']
    assert state[0] == tr.OFF
    assert state[5] == state[10] == tr.ON
    assert state[15] == tr.OFF
    # 漂白不是状态码: 尾巴仍是等级 0, 由 bleach_frame / alive 标记
    assert state[40] == state[199] == tr.OFF
    assert not trace['alive'][40:].any()

    assert trace['n_on'] == 3                 # 闪烁次数
    assert trace['n_switches'] == 5           # OFF→ON 三次 + ON→OFF 两次
    assert trace['bleach_frame'] == 40        # 之后连续 OFF 160 帧 >= 100
    assert trace['on_frames'] == 25
    assert trace['on_fraction'] == pytest.approx(25 / 40)
    # 两个等级 -> state 里只有 0/1, n_cp 与 n_switches 相同
    assert sorted(set(state.tolist())) == [0, 1]
    assert trace['n_states'] == 2
    assert trace['n_cp'] == trace['n_switches']
    assert trace['fit'][5] == pytest.approx(100.0)
    assert trace['fit'][0] == pytest.approx(5.0)


def test_bleach_is_a_time_judgement():
    """同样的 OFF 尾巴: 阈值小就是漂白, 阈值大就不算。"""
    y = _intensity([(5, 40)], n=200)          # 最后一个 ON 在 39, 之后 160 帧 OFF
    short = _states(_trace(y), bleach_frames=100)
    long = _states(_trace(y), bleach_frames=500)

    assert short['bleach_frame'] == 40
    assert long['bleach_frame'] == -1
    assert not short['alive'][40:].any()
    assert long['alive'].all()                # 没判漂白 -> 整条都算活着
    # 没判漂白时, OFF 尾巴仍计入 on_fraction 的分母
    assert long['on_fraction'] == pytest.approx(35 / 200)


def test_no_bleach_when_trace_ends_on():
    y = _intensity([(5, 15), (40, 55)], n=55)
    trace = _states(_trace(y), bleach_frames=10)
    assert trace['bleach_frame'] == -1
    assert trace['state'][-1] == tr.ON
    assert trace['n_on'] == 2


def test_short_off_gap_is_filled():
    """中间只暗了 1 帧: min_off=2 时应并成一次闪烁。"""
    y = _intensity([(5, 15), (16, 25)])
    raw = _states(_trace(y), min_off=1)
    merged = _states(_trace(y), min_off=2)

    assert raw['n_on'] == 2
    assert merged['n_on'] == 1
    assert merged['state'][15] == tr.ON


def test_short_on_spike_is_dropped():
    y = _intensity([(5, 15)])
    y[30] = 100.0                             # 孤立的 1 帧亮
    raw = _states(_trace(y), min_on=1)
    clean = _states(_trace(y), min_on=2)

    assert raw['n_on'] == 2
    assert clean['n_on'] == 1
    assert clean['state'][30] == tr.OFF


def test_parameters_are_validated():
    with pytest.raises(ValueError, match='bleach_frames 至少为 1'):
        tr.ThresholdStates(bleach_frames=0)
    with pytest.raises(ValueError, match='至少为 1'):
        tr.ThresholdStates(min_on=0)
    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        _pipe(10).add_fields(
            [{'emitter_id': 0, 'frames': np.arange(10), 'x': 0.0, 'y': 0.0}],
            [tr.ThresholdStates()])


def test_select_traces_filters_by_field():
    traces = []
    for i, blocks in enumerate([[(5, 15)], [(5, 10), (20, 30)],
                                [(5, 8), (12, 20), (25, 40)]]):
        traces.append(_states(_trace(_intensity(blocks, n=200), emitter_id=i),
                              high=50, bleach_frames=100))

    assert tr.select_traces(traces, n_on=1) == [0]
    assert tr.select_traces(traces, n_on=(2, None)) == [1, 2]
    assert tr.select_traces(traces, n_on=(2, 2)) == [1]
    assert tr.select_traces(traces, bleach_frame=(0, None)) == [0, 1, 2]
    assert tr.select_traces(traces, bleach_frame=(-1, -1)) == []
    assert tr.select_traces(traces, n_switches=lambda v: v > 3) == [2]
    assert tr.select_traces(traces, n_on=(2, None),
                            n_switches=(4, None)) == [2]
    with pytest.raises(KeyError, match='没有字段'):
        tr.select_traces(traces, nope=(0, 1))


def test_bands_draw_state_spans():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    trace = _states(_trace(_intensity([(5, 15), (20, 30)], n=200)),
                    bleach_frames=100)
    fig = tr.plot_trace_fields([trace], indices=[0],
                               groups={'signal': ('intensity',)},
                               bands='state', scalars=('n_on', 'n_switches'))
    ax = fig.axes[0]
    # OFF / ON / OFF / ON 四段色带 + 漂白尾巴的斜纹
    assert len(ax.patches) >= 5
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert 'n_on = 2' in labels
    plt.close(fig)

    with pytest.raises(KeyError, match='没有字段'):
        tr.plot_trace_fields([trace], indices=[0],
                             groups={'s': ('intensity',)}, bands='nope')


def test_state_summary_reads_the_shared_fields():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    traces = [_states(_trace(_intensity([(5, 15)], n=200), emitter_id=i),
                      high=50, bleach_frames=100) for i in range(3)]
    fig = tr.plot_state_summary(traces)
    assert len(fig.axes) == 4
    # 三个面板是直方图, 一个是饼图
    assert len(fig.axes[3].patches) == 2
    plt.close(fig)

    bare = _trace(_intensity([(5, 15)], n=50))     # 没跑过状态判定
    with pytest.raises(KeyError, match='先跑 ThresholdStates'):
        tr.plot_state_summary([bare])


# ---------------------------------------------------------------------------
#  连到 emitter 小视频
# ---------------------------------------------------------------------------
def test_show_trace_movies_maps_trace_indices_to_emitters(monkeypatch):
    import smlm.viz as sv

    calls = {}

    def fake(stack, emitters, idx=None, **kw):
        calls['stack'] = stack
        calls['emitters'] = emitters
        calls['idx'] = list(idx)
        calls['kw'] = kw
        return 'IMAGE'

    monkeypatch.setattr(sv, 'show_emitter_movies', fake)

    traces = [_states(_trace(_intensity([(5, 30)], n=200), emitter_id=i),
                      high=50, bleach_frames=40) for i in range(4)]
    emitters = list(range(10))                   # 只用来核对索引映射

    out = tr.show_trace_movies('STACK', traces, emitters,
                               select={'n_on': (1, None)}, n_max=2,
                               half=7, frames=slice(0, 50))
    assert out == 'IMAGE'
    assert calls['idx'] == [0, 1]
    assert calls['stack'] == 'STACK'
    assert calls['emitters'] is emitters
    assert calls['kw'] == {'half': 7, 'frames': slice(0, 50)}


def test_show_trace_movies_validates_selection(monkeypatch):
    import smlm.viz as sv
    monkeypatch.setattr(sv, 'show_emitter_movies',
                        lambda *a, **k: pytest.fail('不该被调用'))

    traces = [_states(_trace(_intensity([(5, 30)], n=200)), high=50,
                      bleach_frames=40)]
    with pytest.raises(ValueError, match='没有 trace 可展示'):
        tr.show_trace_movies('S', traces, [0], select={'n_on': (5, None)})
    with pytest.raises(IndexError, match='一一对应'):
        tr.show_trace_movies('S', traces, [], indices=[0])


def test_plot_trace_fields_pins_the_bottom_to_zero():
    """下限统一为 0, 上限仍自适应 —— 否则多行零点不同, 看不出有没有信号。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    y = _intensity([(5, 15), (20, 30)], n=200) + 900.0     # 整体抬高, 远离 0
    trace = _states(_trace(y), high=50, bleach_frames=100)

    fig = tr.plot_trace_fields([trace], indices=[0],
                               groups={'s': ('intensity',)})
    bottom, top = fig.axes[0].get_ylim()
    assert bottom == 0.0
    assert top > 1000.0                      # 上限自适应, 没被压到 0
    plt.close(fig)

    fig = tr.plot_trace_fields([trace], indices=[0],
                               groups={'s': ('intensity',)}, ylim_bottom=None)
    assert fig.axes[0].get_ylim()[0] > 0.0   # None = 回到完全自适应
    plt.close(fig)
