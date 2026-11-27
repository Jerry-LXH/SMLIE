import numpy as np
import pytest

from smlm.datatypes import Emitter, Localizations, MultiImageStack
from smlm.multiplex import (DualViewCalibration, DualViewGeometry,
                            analyze_multichannel_emitters, emitter_centers_xy,
                            extract_multichannel_traces, map_channel_xy)
import smlm.trace as tr
from smlm.trace import ApertureRaw, TraceFillPipeline, as_centers_yx


def _calibration(moving, sx, sy):
    return DualViewCalibration(
        poly_a=np.array([0.0, sx, 0.0]),
        poly_b=np.array([0.0, 0.0, sy]),
        center=np.zeros(2),
        scale=1.0,
        affine=np.array([[sx, 0.0], [0.0, sy], [0.0, 0.0]]),
        poly_order=1,
        reference='g',
        moving=moving,
    )


def _multi(height=48, width=48, n_frames=12):
    calibration = _calibration('r', 2.0, 0.5)
    multi = MultiImageStack(
        data=np.zeros((2, n_frames, height, width), dtype=np.float32),
        channels=['g', 'r'], offsets=np.zeros((2, 2), dtype=int),
        pixel_size=100.0, reg={'r': calibration})
    return multi


def _geometry(multi):
    shape = multi.data.shape[-2:]
    windows = {channel: (0, 0, *shape) for channel in multi.channels}
    return DualViewGeometry(
        windows=windows, calib=dict(multi.reg), calib_windows=windows,
        reference='g', check_coverage=False)


def _blinking_locs(positions, n_frames=12, jitter=0.05, seed=0):
    rng = np.random.default_rng(seed)
    frames, xs, ys = [], [], []
    for x, y in positions:
        frames.append(np.arange(n_frames))
        xs.append(x + rng.normal(0, jitter, n_frames))
        ys.append(y + rng.normal(0, jitter, n_frames))
    frames = np.concatenate(frames).astype(int)
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    return Localizations(frame=frames, x=x, y=y,
                         photons=np.full(len(x), 1000.0),
                         bg=np.zeros(len(x)), sigma=np.full(len(x), 1.5),
                         lp=np.full(len(x), 0.2))


def _emitter(x, y, n_frames=12):
    return Emitter(frame=np.arange(n_frames), x=np.full(n_frames, float(x)),
                   y=np.full(n_frames, float(y)),
                   loc_idx=np.arange(n_frames))


LINK_KW = dict(search_radius=1.0, max_gap=2, motion='static',
               min_frames=5, first_frame_range=None, ecc_p=1e-3)


def test_emitter_centers_xy_uses_centroid_order():
    centers = emitter_centers_xy([_emitter(4.0, 8.0), _emitter(11.0, 13.0)])
    np.testing.assert_allclose(centers, [[4.0, 8.0], [11.0, 13.0]])
    assert centers.shape == (2, 2)


def test_analyze_multichannel_emitters_clusters_each_channel():
    multi = _multi()
    positions = [(10.0, 12.0), (20.0, 30.0)]
    locs = {'r': _blinking_locs(positions), 'g': _blinking_locs(positions)}

    emitters = analyze_multichannel_emitters(multi, locs, link_kw=LINK_KW)

    assert set(emitters) == {'g', 'r'}
    assert len(emitters['r']) == len(emitters['g']) == 2
    np.testing.assert_allclose(
        np.sort(emitter_centers_xy(emitters['r']), axis=0),
        np.array(sorted(positions)), atol=0.2)


def test_analyze_multichannel_emitters_reports_missing_channels():
    multi = _multi()
    with pytest.raises(ValueError, match='missing localizations'):
        analyze_multichannel_emitters(multi, {'r': _blinking_locs([(5.0, 5.0)])},
                                      link_kw=LINK_KW)
    with pytest.raises(ValueError, match='unknown channels'):
        analyze_multichannel_emitters(multi, {'b': _blinking_locs([(5.0, 5.0)])},
                                      channels=['b'], link_kw=LINK_KW)


def test_map_channel_xy_matches_geometry_mapping():
    multi = _multi()
    geom = _geometry(multi)
    xy = np.array([[10.0, 12.0], [20.0, 30.0]])
    expected = geom.to_channel(geom.to_reference(xy, 'r'), 'g')
    np.testing.assert_allclose(map_channel_xy(xy, geom, 'r', 'g'), expected)
    np.testing.assert_allclose(map_channel_xy(xy, geom, 'r', 'r'), xy)


def test_extract_multichannel_traces_maps_centres_and_skips_other_channels():
    multi = _multi()
    emitters = {'r': [_emitter(10.0, 12.0), _emitter(20.0, 30.0)],
                'g': [_emitter(20.0, 6.0), _emitter(40.0, 15.0)]}

    traces = extract_multichannel_traces(
        multi, emitters,
        centers_from={'g': 'r', 'r': 'r'},
        mask_from={'g': 'g', 'r': 'r'},
        sigma=1.5, mode='matched', k_sigma=2.0)

    assert set(traces) == {'g', 'r'}
    np.testing.assert_allclose([t.x for t in traces['r']], [10.0, 20.0])
    np.testing.assert_allclose([t.y for t in traces['r']], [12.0, 30.0])
    # r -> g calibration is (x * 2, y * 0.5).
    np.testing.assert_allclose([t.x for t in traces['g']], [20.0, 40.0])
    np.testing.assert_allclose([t.y for t in traces['g']], [6.0, 15.0])
    assert [t.emitter_id for t in traces['g']] == [0, 1]
    # 与单通道同一条阶段链 -> 同一套字段
    assert {'raw', 'bg', 'intensity', 'snr', 'radius', 'sigma'} <= set(traces['g'][0])

    # Channel omitted from centers_from is not extracted.
    only_r = extract_multichannel_traces(
        multi, emitters, centers_from={'r': 'r'},
        sigma=1.5, mode='matched', k_sigma=2.0)
    assert set(only_r) == {'r'}


def test_extract_multichannel_traces_uses_own_emitters_for_background():
    multi = _multi()
    emitters = {'r': [_emitter(10.0, 10.0)],
                'g': [_emitter(10.0, 10.0), _emitter(12.0, 10.0)]}
    kw = dict(mode='hard', radius=2, bg_gap=2, bg_width=3, excl_radius=2)

    self_masked = extract_multichannel_traces(
        multi, emitters, centers_from={'g': 'r'}, **kw)
    own_masked = extract_multichannel_traces(
        multi, emitters, centers_from={'g': 'r'}, mask_from={'g': 'g'}, **kw)

    # The second green emitter sits 2 px away and must be cut out of the ring.
    assert self_masked['g'][0]['neighbours'] == 0
    assert own_masked['g'][0]['neighbours'] >= 1


def test_extract_multichannel_traces_validates_sources():
    multi = _multi()
    emitters = {'r': [_emitter(10.0, 10.0)]}
    with pytest.raises(ValueError, match='unknown output channels'):
        extract_multichannel_traces(multi, emitters,
                                    centers_from={'b': 'r'})
    with pytest.raises(ValueError, match='no entry for source channel'):
        extract_multichannel_traces(multi, emitters, centers_from={'r': 'g'})
    with pytest.raises(ValueError, match='no entry for source channel'):
        extract_multichannel_traces(multi, emitters, centers_from={'r': 'r'},
                                    mask_from={'r': 'g'})
    with pytest.raises(TypeError, match='MultiImageStack'):
        extract_multichannel_traces((32, 32), emitters,
                                    centers_from={'r': 'r'})


def test_pipeline_accepts_centres_and_emitters_identically():
    rng = np.random.default_rng(1)
    movie = rng.random((6, 32, 32)).astype(np.float32)
    emitters = [_emitter(12.3, 9.7, n_frames=6), _emitter(20.1, 18.4, 6)]
    centers_yx = np.array([[9.7, 12.3], [18.4, 20.1]])
    fillers = [ApertureRaw(mode='hard', radius=2),
               tr.AnnulusRing(gap=2, width=3),
               tr.IntensityFromRawBg()]

    from_emitters = TraceFillPipeline(movie, as_centers_yx(emitters),
                                      sigma=1.5).run(fillers)
    from_centers = TraceFillPipeline(movie, centers_yx, sigma=1.5).run(fillers)

    assert len(from_emitters) == len(from_centers) == 2
    for left, right in zip(from_emitters, from_centers):
        np.testing.assert_allclose(left['intensity'], right['intensity'])
        np.testing.assert_allclose(left['raw'], right['raw'])
        assert (left['x'], left['y']) == pytest.approx((right['x'], right['y']))


def test_pipeline_writes_snr_and_honours_centers_for_bg():
    rng = np.random.default_rng(3)
    movie = (rng.random((8, 32, 32)) + 5.0).astype(np.float32)
    movie[:, 12, 12] += 200.0             # bright spot at (row, col) = (12, 12)
    centers_yx = np.array([[12.0, 12.0]])
    fillers = [ApertureRaw(mode='hard', radius=2),
               tr.AnnulusRing(gap=2, width=3, excl_radius=2),
               tr.IntensityFromRawBg(), tr.PoissonSNR()]

    own = TraceFillPipeline(movie, centers_yx, sigma=1.5).run(fillers)
    masked = TraceFillPipeline(movie, centers_yx, sigma=1.5,
                               centers_for_bg=np.array([[12.0, 14.0]])
                               ).run(fillers)

    assert own[0]['snr'] is not None
    assert np.all(own[0]['snr'] > 0)
    assert own[0]['neighbours'] == 0
    # The extra mask centre sits inside the ring, so it is cut out.
    assert masked[0]['neighbours'] >= 1


def test_per_channel_link_kw_overrides_the_common_ones():
    """绿通道只提供挖洞位置, 允许更短的轨迹; 红通道走公共的严格设置。"""
    multi = _multi(n_frames=12)
    locs = {'r': _blinking_locs([(10.0, 10.0)], n_frames=12),
            'g': _blinking_locs([(10.0, 10.0)], n_frames=2)}
    common = dict(search_radius=1.0, max_gap=2, motion='static',
                  min_frames=5, first_frame_range=None)

    both = analyze_multichannel_emitters(multi, locs, link_kw=common)
    assert len(both['r']) == 1
    assert len(both['g']) == 0                  # 2 个定位点 < min_frames=5

    tuned = analyze_multichannel_emitters(
        multi, locs, link_kw=common,
        link_kw_per_channel={'g': dict(min_frames=1)})
    assert len(tuned['r']) == 1
    assert len(tuned['g']) == 1                 # 逐通道放宽

    # 覆盖是"叠在公共参数上", 没写的键仍然继承
    assert tuned['r'][0].n_locs == both['r'][0].n_locs


def test_per_channel_link_kw_validates_channels():
    multi = _multi()
    locs = {'r': _blinking_locs([(10.0, 10.0)]),
            'g': _blinking_locs([(10.0, 10.0)])}
    with pytest.raises(ValueError, match='unknown channels'):
        analyze_multichannel_emitters(multi, locs,
                                      link_kw_per_channel={'b': dict()})
    with pytest.raises(TypeError, match='must be a dict'):
        analyze_multichannel_emitters(multi, locs,
                                      link_kw_per_channel={'g': 5})


def test_return_all_exposes_the_pre_filter_set():
    multi = _multi(n_frames=12)
    locs = {'r': _blinking_locs([(10.0, 10.0)], n_frames=12),
            'g': _blinking_locs([(10.0, 10.0)], n_frames=2)}
    common = dict(search_radius=1.0, max_gap=2, motion='static',
                  min_frames=5, first_frame_range=None)

    kept, every = analyze_multichannel_emitters(
        multi, locs, link_kw=common, return_all=True)
    assert set(kept) == set(every) == {'r', 'g'}
    assert len(every['g']) == 1 and len(kept['g']) == 0


def _fake_traces(n=3, n_frames=40, emitter_ids=True):
    """最小可用的 {通道: [Trace]}, 只要字段齐全。"""
    import smlm.trace as tr

    def one(ch, i):
        x = np.arange(n_frames)
        base = 100.0 + 10 * i + (500.0 if ch == 'r' else 300.0) * (x > 10)
        return tr.Trace(emitter_id=i if emitter_ids else 0,
                        frames=x, intensity=base, raw=base + 20.0,
                        bg=np.full(n_frames, 20.0 if ch == 'r' else 12.0),
                        radius=3, neighbours=i)

    return {'g': [one('g', i) for i in range(n)],
            'r': [one('r', i) for i in range(n)]}


def test_plot_channel_traces_overlays_channels_in_each_group():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from smlm.multiplex import plot_channel_traces

    traces = _fake_traces()
    fig = plot_channel_traces(
        traces, channels=('g', 'r'),
        colors={'g': 'tab:green', 'r': 'tab:red'},
        groups={'signal': dict(fields=('intensity',)),
                'raw + bg': dict(fields=('raw', 'bg'))},
        indices=range(3), scalars={'g': ('radius',), 'r': ('neighbours',)})

    assert len(fig.axes) == 6                     # 3 行 x 2 组
    # 组 'signal': 每图 2 条线 (两个通道)
    for k in (0, 2, 4):
        assert len(fig.axes[k].get_lines()) == 2
    first = fig.axes[0].get_lines()
    assert first[0].get_color() == 'tab:green' and first[1].get_color() == 'tab:red'

    # 组 'raw + bg': 每图 4 条线 (2 通道 x 2 字段), 同通道同色, bg 更暗
    for k in (1, 3, 5):
        assert len(fig.axes[k].get_lines()) == 4
    raw_bg = fig.axes[1].get_lines()
    alphas = [ln.get_alpha() for ln in raw_bg]
    assert alphas[0] == pytest.approx(1.0)
    assert alphas[1] < alphas[0]                   # bg 比 raw 暗
    assert raw_bg[0].get_color() == raw_bg[1].get_color()   # 同通道同色

    labels = [t.get_text() for t in fig.axes[1].get_legend().get_texts()]
    assert any(t.startswith('g · raw') for t in labels)
    assert any(t.startswith('g · bg') for t in labels)
    assert any(t.startswith('g radius') for t in labels)      # 标量也进了图例
    plt.close(fig)


def test_plot_channel_traces_validates_inputs():
    import matplotlib
    matplotlib.use('Agg')
    from smlm.multiplex import plot_channel_traces

    traces = _fake_traces()
    with pytest.raises(ValueError, match='未知通道'):
        plot_channel_traces(traces, channels=('b',))
    with pytest.raises(ValueError, match='索引必须对齐'):
        plot_channel_traces({'g': traces['g'], 'r': traces['r'][:2]})
    with pytest.raises(ValueError, match='未选中的通道'):
        plot_channel_traces(traces, colors={'b': 'k'})
    with pytest.raises(KeyError, match='没有字段'):
        plot_channel_traces(traces, groups={'x': ('nope',)})
    with pytest.raises(IndexError, match='超出'):
        plot_channel_traces(traces, indices=[99])


def test_plot_channel_traces_pins_the_bottom_to_zero():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from smlm.multiplex import plot_channel_traces

    traces = _fake_traces()
    for ch in traces:                        # 整体抬高, 远离 0
        for t in traces[ch]:
            t['intensity'] = t['intensity'] + 500.0

    fig = plot_channel_traces(traces, groups={'signal': ('intensity',)},
                              indices=range(2))
    for ax in fig.axes:
        assert ax.get_ylim()[0] == 0.0
        assert ax.get_ylim()[1] > 500.0
    plt.close(fig)

    fig = plot_channel_traces(traces, groups={'signal': ('intensity',)},
                              indices=range(2), ylim_bottom=None)
    assert fig.axes[0].get_ylim()[0] > 0.0
    plt.close(fig)
