"""TraceFillPipeline 的最小用法: 骨架 + ROI raw。"""
import numpy as np
import pytest

import smlm.trace as tr
from smlm.datatypes import Emitter

N_FRAMES = 10
POSITIONS = [(12.35, 20.62), (25.10, 14.90)]     # (y, x)


def _movie(size=40, n_frames=N_FRAMES, sigma=1.5, seed=0):
    rng = np.random.default_rng(seed)
    movie = (rng.normal(50.0, 2.0, (n_frames, size, size))).astype(np.float32)
    yy, xx = np.mgrid[0:size, 0:size]
    frames = np.arange(n_frames)
    for k, (cy, cx) in enumerate(POSITIONS):
        on = ((frames >= 3 * k) & (frames < 3 * k + 7)).astype(np.float32)
        movie += (800 * on[:, None, None]
                  * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2)
                           / (2 * sigma ** 2))).astype(np.float32)
    return movie


def _on_frames(index):
    """``_movie`` 里第 ``index`` 个分子的 ON 帧。"""
    frames = np.arange(N_FRAMES)
    return (frames >= 3 * index) & (frames < 3 * index + 7)


def _emitters(n_frames=N_FRAMES):
    return [Emitter(frame=np.arange(n_frames), x=np.full(n_frames, cx),
                    y=np.full(n_frames, cy), loc_idx=np.arange(n_frames))
            for cy, cx in POSITIONS]


def test_pipeline_builds_one_trace_per_emitter():
    centers = np.array(POSITIONS)
    traces = tr.TraceFillPipeline(_movie(), centers, sigma=1.5).run([tr.ApertureRaw()])

    assert len(traces) == len(POSITIONS)
    assert [t['emitter_id'] for t in traces] == [0, 1]
    assert [(t['y'], t['x']) for t in traces] == POSITIONS
    np.testing.assert_array_equal(traces[0]['frames'], np.arange(N_FRAMES))
    # raw 阶段只写信号相关的键, 不写背景/强度/SNR
    assert {'raw', 'sw', 'sw2', 'area', 'radius', 'sigma', 'mode'} <= set(traces[0])
    assert not {'bg', 'intensity', 'snr'} & set(traces[0])


def test_raw_equals_weighted_sum_over_the_roi():
    sigma = 1.5
    centers = np.array(POSITIONS)
    movie = _movie()
    traces = tr.TraceFillPipeline(movie, centers, sigma=[sigma, sigma]
                                  ).run([tr.ApertureRaw(k_sigma=2.0)])

    spec = tr.ApertureSpec(mode='matched', k_sigma=2.0)
    for trace, (cy, cx) in zip(traces, POSITIONS):
        half = int(round(cy)), int(round(cx))
        ap = tr.Aperture(spec, sigma, center_frac=(cy - half[0], cx - half[1]))
        box, _, _ = tr.crop_centered(movie, cy, cx, ap.radius)
        expected = np.einsum('tij,ij->t', box, ap.w) / ap.sw2
        np.testing.assert_array_equal(trace['raw'], expected)


def test_emitter_sequence_and_centre_array_agree():
    movie = _movie()
    from_emitters = tr.TraceFillPipeline(movie, _emitters(), sigma=1.5).run([tr.ApertureRaw()])
    from_centers = tr.TraceFillPipeline(
        movie, np.array(POSITIONS), sigma=1.5).run([tr.ApertureRaw()])
    for a, b in zip(from_emitters, from_centers):
        # Emitter.centroid 是逐帧坐标的均值, 可能比原坐标差 1 ULP
        np.testing.assert_allclose(a['raw'], b['raw'], rtol=1e-12)
        assert a['x'] == pytest.approx(b['x'])
        assert a['y'] == pytest.approx(b['y'])


def test_pipeline_rejects_inconsistent_stage_chains():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)

    class NeedsBackground(tr.TraceFiller):
        provides = ('snr',)
        requires = ('bg', 'raw')

        def fill(self, trace, ctx, box):
            pass

    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        pipe.run([tr.ApertureRaw(), NeedsBackground()])

    class UnknownField(tr.TraceFiller):
        provides = ('not_a_field',)

        def fill(self, trace, ctx, box):
            trace['not_a_field'] = 0

    with pytest.raises(ValueError, match='未登记字段'):
        pipe.run([UnknownField()])


def _gaussian_movie(amplitude=600.0, sigma=1.5, n_frames=6, size=48,
                    center=(24, 24)):
    yy, xx = np.mgrid[0:size, 0:size]
    frame = amplitude * np.exp(-((yy - center[0]) ** 2 + (xx - center[1]) ** 2)
                               / (2 * sigma ** 2))
    return np.broadcast_to(frame, (n_frames, size, size)).astype(
        np.float32).copy()


def test_aperture_raw_matched_equals_peak_times_2pi_sigma_squared():
    """matched 核与 PSF 同宽时, 归一化正好抵消截断: raw = A·2πσ²。"""
    sigma, amp = 1.5, 600.0
    movie = _gaussian_movie(amp, sigma)
    trace = tr.TraceFillPipeline(movie, np.array([[24.0, 24.0]]),
                                 sigma=sigma).run(
        [tr.ApertureRaw(mode='matched', k_sigma=2.0)])[0]

    assert trace['radius'] == 3                       # ceil(2 x 1.5)
    assert trace['raw'][0] == pytest.approx(amp * 2 * np.pi * sigma ** 2,
                                            rel=1e-4)


def test_aperture_raw_hard_equals_the_disc_sum():
    """hard 圆盘 sw2 = 1, 所以 raw 就是圆内像素直接求和。"""
    sigma, amp = 1.5, 600.0
    movie = _gaussian_movie(amp, sigma)
    trace = tr.TraceFillPipeline(movie, np.array([[24.0, 24.0]]),
                                 sigma=sigma).run(
        [tr.ApertureRaw(mode='hard', radius=3)])[0]

    yy, xx = np.mgrid[0:48, 0:48]
    disc = (yy - 24) ** 2 + (xx - 24) ** 2 <= 3 ** 2 + 1e-9
    assert trace['raw'][0] == pytest.approx(float(movie[0][disc].sum()),
                                            rel=1e-5)
    assert trace['sw2'] == 1.0
    assert trace['area'] == float(disc.sum())


def test_to_table_stacks_arrays_and_scalars():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw()])
    table = tr.to_table(traces)

    assert table['raw'].shape == (2, N_FRAMES)
    assert table['radius'].shape == (2,)
    np.testing.assert_array_equal(table['raw'][1], traces[1]['raw'])


# ---------------------------------------------------------------------------
#  背景阶段
# ---------------------------------------------------------------------------
def test_annulus_ring_recovers_the_flat_background():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw(), tr.AnnulusRing(gap=2, width=4)])

    for t in traces:
        assert t['n_bg_pixels'] > 50
        np.testing.assert_allclose(t['bg'], 50.0, atol=1.5)
        assert np.all(np.isfinite(t['bg_std']))
        assert np.all(t['bg_std'] > 0)
        # 环必须落在信号框外
        assert t['n_bg_pixels'] < (2 * t['radius'] + 2 * 7) ** 2


def test_annulus_ring_punches_out_neighbours():
    """环里另一个分子的核心应被挖掉, 否则背景被抬高。"""
    movie = _movie()
    target = np.array([POSITIONS[0]])
    ring = tr.AnnulusRing(gap=2, width=4, excl_radius=2)

    alone = tr.TraceFillPipeline(movie, target, sigma=1.5).run(
        [tr.ApertureRaw(), ring])[0]

    yy, xx = np.mgrid[0:movie.shape[-1], 0:movie.shape[-1]]
    neighbour = (19.35, 20.62)                   # 距 target 7 px, 落在环带里
    movie = movie + (600 * np.exp(-((yy - neighbour[0]) ** 2
                                    + (xx - neighbour[1]) ** 2) / (2 * 1.5 ** 2))
                     ).astype(np.float32)[None]
    masked = tr.TraceFillPipeline(movie, target, sigma=1.5,
                                  centers_for_bg=np.array([neighbour])).run(
        [tr.ApertureRaw(), ring])[0]
    unmasked = tr.TraceFillPipeline(movie, target, sigma=1.5).run(
        [tr.ApertureRaw(), ring])[0]

    assert alone['neighbours'] == 0
    assert masked['neighbours'] == 1
    assert masked['n_bg_pixels'] < unmasked['n_bg_pixels']
    assert abs(np.median(masked['bg']) - 50.0) < abs(np.median(unmasked['bg']) - 50.0)


def test_annulus_ring_reports_a_flat_background_exactly():
    """常数背景 -> bg 精确等于该常数, bg_std 精确为 0。"""
    level = 7.5
    movie = np.full((5, 40, 40), level, np.float32)
    trace = tr.TraceFillPipeline(movie, np.array([[20.0, 20.0]]), sigma=1.5).run(
        [tr.ApertureRaw(mode='hard', radius=2),
         tr.AnnulusRing(gap=2, width=3, excl_radius=2)])[0]

    np.testing.assert_allclose(trace['bg'], level)
    np.testing.assert_allclose(trace['bg_std'], 0.0, atol=1e-6)
    assert trace['neighbours'] == 0
    assert trace['n_bg_pixels'] > 0


def test_annulus_ring_must_follow_a_signal_stage():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        pipe.run([tr.AnnulusRing()])


def test_plot_trace_fields_selects_fields_and_validates():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw(), tr.AnnulusRing()])

    fig = tr.plot_trace_fields(traces, fields=('raw', 'bg'), indices=range(2))
    assert len(fig.axes) == 4                 # 2 条 trace x 2 个字段
    plt.close(fig)

    with pytest.raises(ValueError, match='不是逐帧序列'):
        tr.plot_trace_fields(traces, fields=('radius',))
    with pytest.raises(KeyError, match='没有字段'):
        tr.plot_trace_fields(traces, fields=('nope',))


def test_plot_trace_fields_groups_overlay_and_colours():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw(), tr.AnnulusRing()])

    # 两个分组 -> 两列; 组内叠加 -> 同一张图上两条线
    fig = tr.plot_trace_fields(traces, indices=range(2), groups={
        'signal': dict(fields=('raw',), color='tab:blue'),
        'background': dict(fields=('bg', 'bg_std'), color='tab:red')})
    assert len(fig.axes) == 4
    bg_ax = fig.axes[1]
    assert len(bg_ax.get_lines()) == 2
    assert {ln.get_color() for ln in bg_ax.get_lines()} == {'tab:red'}
    # 一个颜色给多条线 -> 用线型区分
    assert {ln.get_linestyle() for ln in bg_ax.get_lines()} == {'-', '--'}
    assert bg_ax.get_title() == '#0  background'
    plt.close(fig)

    # 逐字段上色
    fig = tr.plot_trace_fields(
        traces[:1], groups=[dict(fields=('bg', 'bg_std'),
                                 color=('tab:red', 'tab:orange'))])
    assert [ln.get_color() for ln in fig.axes[0].get_lines()] == \
        ['tab:red', 'tab:orange']
    plt.close(fig)

    # 颜色数与字段数不符
    with pytest.raises(ValueError, match='颜色'):
        tr.plot_trace_fields(traces[:1],
                             groups=[dict(fields=('bg', 'bg_std'),
                                          color=('r', 'g', 'b'))])


# ---------------------------------------------------------------------------
#  计算阶段: intensity / snr
# ---------------------------------------------------------------------------
def _full_chain():
    return [tr.ApertureRaw(k_sigma=2.0), tr.AnnulusRing(gap=2, width=4),
            tr.IntensityFromRawBg(), tr.PoissonSNR()]


def test_intensity_is_raw_minus_background_flux():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        _full_chain())

    for i, t in enumerate(traces):
        np.testing.assert_allclose(t['bg_flux'], t['bg'] * t['sw'] / t['sw2'],
                                   rtol=1e-15)
        np.testing.assert_allclose(t['intensity'],
                                   np.asarray(t['raw']) - t['bg_flux'],
                                   rtol=1e-15)
        on = _on_frames(i)
        assert np.median(t['intensity'][on]) > 100
        assert abs(np.median(t['intensity'][~on])) < 100


def test_snr_matches_its_definition():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        _full_chain())
    for t in traces:
        total = np.asarray(t['intensity'], float) * t['sw2']
        area_eff = t['sw'] ** 2 / t['sw2']
        noise = np.sqrt(np.maximum(total, 0.0) + area_eff * t['bg_std'] ** 2)
        expected = np.where((total > 0) & (t['bg_std'] > 0),
                            total / np.maximum(noise, 1e-12), 0.0)
        np.testing.assert_array_equal(t['snr'], expected)
        assert np.all(t['snr'] >= 0)


def test_intensity_and_bg_flux_follow_their_definition():
    """intensity = raw - bg*sw/sw2, bg_flux = bg*sw/sw2 (定义式, 不依赖旧实现)。"""
    movie = _movie()
    traces = tr.TraceFillPipeline(movie, np.array(POSITIONS),
                                  sigma=1.5).run(_full_chain())

    for t in traces:
        bg_flux = np.asarray(t['bg'], float) * (t['sw'] / t['sw2'])
        np.testing.assert_allclose(t['bg_flux'], bg_flux, rtol=1e-12)
        np.testing.assert_allclose(t['intensity'],
                                   np.asarray(t['raw'], float) - bg_flux,
                                   rtol=1e-12, atol=1e-9)


def test_calculate_stages_require_their_inputs():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        pipe.run([tr.ApertureRaw(), tr.PoissonSNR()])
    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        pipe.run([tr.ApertureRaw(), tr.AnnulusRing(), tr.PoissonSNR()])


def test_standard_fillers_builds_the_usual_chain():
    fillers = tr.standard_fillers(tr.ApertureSpec(mode='hard', radius=3))
    assert [type(f).__name__ for f in fillers] == [
        'ApertureRaw', 'AnnulusRing', 'IntensityFromRawBg', 'PoissonSNR']
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(fillers)
    assert {'raw', 'bg', 'intensity', 'snr'} <= set(traces[0])

    with pytest.raises(ValueError, match='只支持逐帧背景'):
        tr.standard_fillers(tr.ApertureSpec(bg_ref='tail'))


# ---------------------------------------------------------------------------
#  事后追加字段
# ---------------------------------------------------------------------------
def test_add_fields_in_two_steps_matches_one_step():
    movie = _movie()
    centers = np.array(POSITIONS)
    tail = [tr.AnnulusRing(gap=2, width=4), tr.IntensityFromRawBg(),
            tr.PoissonSNR()]

    one_step = tr.TraceFillPipeline(movie, centers, sigma=1.5).run(
        [tr.ApertureRaw(k_sigma=2.0)] + tail)

    pipe = tr.TraceFillPipeline(movie, centers, sigma=1.5)
    two_step = pipe.run([tr.ApertureRaw(k_sigma=2.0)])
    assert 'bg' not in two_step[0]
    returned = pipe.add_fields(two_step, tail)

    assert returned == two_step                  # 同一批对象, 原地修改
    assert all(a is b for a, b in zip(returned, two_step))
    for a, b in zip(one_step, two_step):
        for field in ('raw', 'bg', 'bg_std', 'intensity', 'snr', 'bg_flux'):
            np.testing.assert_array_equal(a[field], b[field])


def test_add_fields_can_add_a_registered_custom_field():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    traces = pipe.run([tr.ApertureRaw(), tr.AnnulusRing()])

    class DutyCycle(tr.TraceFiller):
        provides = ('duty',)
        requires = ('raw',)

        def fill(self, trace, ctx, box):
            trace['duty'] = float(np.mean(trace['raw'] > np.median(trace['raw'])))

    with pytest.raises(ValueError, match='未登记字段'):
        pipe.add_fields(traces, [DutyCycle()])

    tr.register_field('duty', 'scalar')
    pipe.add_fields(traces, [DutyCycle()])
    assert 0.0 <= traces[0]['duty'] <= 1.0
    assert all('duty' in t for t in traces)


def test_add_fields_rejects_unmet_requirements():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    traces = pipe.run([tr.ApertureRaw()])

    with pytest.raises(ValueError, match='还没有任何阶段提供'):
        pipe.add_fields(traces, [tr.PoissonSNR()])


def test_duplicate_provides_warns_and_overwrites():
    """重复提供字段: 警告 + 用后面的值覆盖。"""
    movie = _movie()
    centers = np.array(POSITIONS)

    with pytest.warns(UserWarning, match='覆盖已存在的字段'):
        traces = tr.TraceFillPipeline(movie, centers, sigma=1.5).run(
            [tr.ApertureRaw(mode='hard', radius=2), tr.ApertureRaw(k_sigma=2.0)])

    reference = tr.TraceFillPipeline(movie, centers, sigma=1.5).run(
        [tr.ApertureRaw(k_sigma=2.0)])
    for got, want in zip(traces, reference):
        assert got['mode'] == 'matched'
        assert got['radius'] == want['radius']          # 3, 而不是 hard 的 2
        np.testing.assert_array_equal(got['raw'], want['raw'])


def test_add_fields_overwrites_existing_field_with_warning():
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    traces = pipe.run([tr.ApertureRaw(mode='hard', radius=2)])
    before = np.asarray(traces[0]['raw'], float)

    with pytest.warns(UserWarning, match='覆盖已存在的字段'):
        pipe.add_fields(traces, [tr.ApertureRaw(k_sigma=2.0)])

    assert traces[0]['mode'] == 'matched'
    assert not np.array_equal(np.asarray(traces[0]['raw'], float), before)


# ---------------------------------------------------------------------------
#  环里有没有邻居: 挖洞 (excl_radius) 对背景的影响
# ---------------------------------------------------------------------------
def _movie_with_ring_neighbour():
    """在 POSITIONS[0] 的环带里再放一个亮分子。"""
    movie = _movie()
    yy, xx = np.mgrid[0:movie.shape[-1], 0:movie.shape[-1]]
    neighbour = (19.35, 20.62)                   # 距 target 7 px
    movie = movie + (600 * np.exp(-((yy - neighbour[0]) ** 2
                                    + (xx - neighbour[1]) ** 2) / (2 * 1.5 ** 2))
                     ).astype(np.float32)[None]
    return movie, np.array([neighbour])


def test_excl_radius_removes_neighbour_cores():
    movie, neighbour = _movie_with_ring_neighbour()
    target = np.array([POSITIONS[0]])

    def run(excl):
        pipe = tr.TraceFillPipeline(movie, target, sigma=1.5,
                                    centers_for_bg=neighbour)
        return pipe.run([tr.ApertureRaw(),
                         tr.AnnulusRing(excl_radius=excl)])[0]

    masked = run(2)
    plain = run(0)

    assert masked['neighbours'] == 1             # 邻居被数出来了
    assert plain['neighbours'] == 0              # excl_radius=0: 一个都没挖
    assert masked['n_bg_pixels'] < plain['n_bg_pixels']
    assert (abs(np.median(masked['bg']) - 50.0)
            < abs(np.median(plain['bg']) - 50.0))


def test_annulus_ring_reports_zero_neighbours_when_isolated():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw(), tr.AnnulusRing()])
    assert all(t['neighbours'] == 0 for t in traces)
    assert all(t['n_bg_pixels'] > 50 for t in traces)


# ---------------------------------------------------------------------------
#  可视化: 单个数放进图例
# ---------------------------------------------------------------------------
def test_plot_trace_fields_shows_scalars_in_legend():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5).run(
        [tr.ApertureRaw(), tr.AnnulusRing()])

    fig = tr.plot_trace_fields(traces, indices=[0],
                               groups={'background': ('bg', 'bg_std')},
                               scalars=('neighbours', 'n_bg_pixels'))
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert 'neighbours = 0' in labels
    assert any(text.startswith('n_bg_pixels = ') for text in labels)
    # 单字段那列也会因为标量而出现图例
    plt.close(fig)

    fig = tr.plot_trace_fields(traces, indices=[0], groups={'bg': ('bg',)},
                               legend=False, scalars='radius')
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels == ['radius = 3']              # 只列标量, 不列字段名
    plt.close(fig)

    # 单字段 + 不给 scalars: 字段名已经在 y 轴标签里, 不画图例
    fig = tr.plot_trace_fields(traces, indices=[0], groups={'bg': ('bg',)})
    assert fig.axes[0].get_legend() is None
    plt.close(fig)

    with pytest.raises(ValueError, match='不是单个数字'):
        tr.plot_trace_fields(traces, indices=[0], groups={'bg': ('bg',)},
                             scalars=('bg',))
    with pytest.raises(KeyError, match='没有字段'):
        tr.plot_trace_fields(traces, indices=[0], groups={'bg': ('bg',)},
                             scalars=('nope',))


# ---------------------------------------------------------------------------
#  σ 与两个中心表
# ---------------------------------------------------------------------------
def test_sigma_scalar_broadcasts_and_sets_radius():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS),
                                  sigma=2.0).run([tr.ApertureRaw(k_sigma=2.0)])
    assert {t['sigma'] for t in traces} == {2.0}
    assert {t['radius'] for t in traces} == {4}        # ceil(2 x 2.0)


def test_sigma_array_is_per_molecule():
    traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS),
                                  sigma=[1.0, 2.0]).run(
        [tr.ApertureRaw(k_sigma=2.0)])
    assert [t['sigma'] for t in traces] == [1.0, 2.0]
    assert [t['radius'] for t in traces] == [2, 4]     # ceil(k x sigma), 下界 2


def test_missing_sigma_warns_and_falls_back_to_default():
    with pytest.warns(UserWarning, match='没有指定 PSF 宽度'):
        traces = tr.TraceFillPipeline(_movie(), np.array(POSITIONS),
                                      sigma_default=1.8).run(
            [tr.ApertureRaw()])
    assert {t['sigma'] for t in traces} == {1.8}


def test_sigma_length_must_match_the_number_of_centres():
    with pytest.raises(ValueError, match='sigma 长度'):
        tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=[1.0, 2.0, 3.0])


def test_centers_for_bg_defaults_to_centers_for_raw():
    """不给 centers_for_bg 时, 背景屏蔽用自己那张表。"""
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5)
    assert pipe.self_mask is True
    np.testing.assert_array_equal(pipe.centers_for_bg, pipe.centers_for_raw)

    wide = np.vstack([np.array(POSITIONS), [[30.0, 30.0]]])
    pipe = tr.TraceFillPipeline(_movie(), np.array(POSITIONS), sigma=1.5,
                                centers_for_bg=wide)
    assert pipe.self_mask is False
    assert len(pipe.centers_for_bg) == 3
    assert len(pipe.centers_for_raw) == 2


def test_annulus_ring_smoothing_reduces_frame_to_frame_jitter():
    """bg_smooth=w 是按 √w 压帧间抖动; 默认 0 = 不平滑 (旧行为不变)。"""
    rng = np.random.default_rng(0)
    movie = (10.0 + rng.normal(0, 1.0, (200, 40, 40))).astype(np.float32)
    centers = np.array([[20.0, 20.0]])

    def run(smooth):
        return tr.TraceFillPipeline(movie, centers, sigma=1.5).run(
            [tr.ApertureRaw(mode='hard', radius=2),
             tr.AnnulusRing(gap=2, width=3, bg_smooth=smooth)])[0]

    plain, smoothed = run(0), run(10)
    jitter = lambda t: float(np.std(np.diff(t["bg"])))
    assert jitter(smoothed) < 0.5 * jitter(plain)
    assert np.median(smoothed["bg"]) == pytest.approx(
        np.median(plain["bg"]), abs=0.3)
