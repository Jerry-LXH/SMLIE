import numpy as np
import pytest

import smlm.drift as drift_module
from smlm.datatypes import DriftTrace, Localizations, MultiImageStack
from smlm.multiplex import (DualViewCalibration, DualViewGeometry,
                            correct_multi, drift_crop_window,
                            drift_requirements, estimate_drift_multi)
from smlm.multiplex.drift import CropWindow


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


def _multi_and_geometry(height=64, width=64):
    channels = ['g', 'r', 'b']
    calibrations = {
        'r': _calibration('r', 2.0, 0.5),
        'b': _calibration('b', 4.0, 2.0),
    }
    multi = MultiImageStack(
        data=np.ones((3, 3, height, width), dtype=np.float32),
        channels=channels,
        offsets=np.zeros((3, 2), dtype=int),
        pixel_size=100.0,
        reg=calibrations,
    )
    windows = {channel: (0, 0, height, width) for channel in channels}
    geom = DualViewGeometry(
        windows=windows,
        calib=calibrations,
        calib_windows=windows,
        reference='g',
        check_coverage=False,
    )
    return multi, geom


def _locs(x, y, lp=10.0):
    return Localizations(
        frame=np.array([0, 1, 2]),
        x=np.full(3, x, dtype=float),
        y=np.full(3, y, dtype=float),
        photons=np.ones(3),
        bg=np.zeros(3),
        lp=np.full(3, lp, dtype=float),
    )


def _traces():
    return {
        'g': DriftTrace(dx=np.array([0.0, 2.0, 4.0]),
                        dy=np.array([0.0, 1.0, 2.0])),
        'r': DriftTrace(dx=np.array([0.0, 1.0, 2.0]),
                        dy=np.array([0.0, 2.0, 4.0])),
        'b': DriftTrace(dx=np.array([0.0, 0.5, 1.0]),
                        dy=np.array([0.0, 0.5, 1.0])),
    }


def test_drift_requirements_are_per_side_and_per_channel():
    per_channel, need = drift_requirements(_traces(), (64, 64))

    assert per_channel['g'] == {'top': 0, 'bottom': 2, 'left': 0, 'right': 4}
    assert per_channel['r'] == {'top': 0, 'bottom': 4, 'left': 0, 'right': 2}
    assert per_channel['b'] == {'top': 0, 'bottom': 1, 'left': 0, 'right': 1}
    # Combined = per-side maximum across channels.
    assert need == {'top': 0, 'bottom': 4, 'left': 0, 'right': 4}


def test_drift_requirements_cover_negative_drift_asymmetrically():
    trace = DriftTrace(dx=np.array([0.0, -3.0, -1.0]),
                       dy=np.array([0.0, -2.0, 5.0]))
    per_channel, need = drift_requirements({'g': trace}, (20, 30))

    assert per_channel['g'] == {'top': 2, 'bottom': 5,
                                'left': 3, 'right': 0}
    assert need == per_channel['g']


def test_drift_crop_window_adds_pad_and_trim_to_every_side():
    window = drift_crop_window(_traces(), (64, 64), pad=2, extra_trim=3)

    assert isinstance(window, CropWindow)
    assert window.margins == {'top': 5, 'bottom': 9, 'left': 5, 'right': 9}
    assert window.y0 == 5 and window.x0 == 5
    assert window.shape == (50, 50)
    assert window.y == (5, 55) and window.x == (5, 55)


def test_drift_crop_window_rejects_impossible_request():
    with pytest.raises(ValueError, match='no common area remains'):
        drift_requirements(_traces(), (4, 64))
    trace = DriftTrace(dx=np.zeros(3), dy=np.zeros(3))
    with pytest.raises(ValueError, match='no valid area'):
        drift_crop_window({'g': trace}, (5, 5), pad=2, extra_trim=1)


def test_drift_crop_window_accepts_per_side_pad_and_trim():
    window = drift_crop_window(_traces(), (64, 64),
                               pad={'top': 0, 'bottom': 2},
                               extra_trim=(1, 2, 3, 4))

    # drift needs: top 0, bottom 4, left 0, right 4
    assert window.margins == {'top': 1, 'bottom': 8, 'left': 3, 'right': 8}
    assert window.shape == (55, 53)
    assert window.y0 == 1 and window.x0 == 3


def test_drift_crop_window_rejects_bad_per_side_values():
    with pytest.raises(ValueError, match='unknown sides'):
        drift_crop_window(_traces(), (64, 64), pad={'up': 1})
    with pytest.raises(ValueError, match='4 values'):
        drift_crop_window(_traces(), (64, 64), extra_trim=(1, 2, 3))
    with pytest.raises(ValueError, match='non-negative'):
        drift_crop_window(_traces(), (64, 64), pad={'top': -1})


def test_drift_crop_window_rejects_non_finite_trace():
    trace = DriftTrace(dx=np.array([0.0, np.nan, 1.0]),
                       dy=np.zeros(3))
    with pytest.raises(ValueError, match='non-finite'):
        drift_requirements({'g': trace}, (64, 64))


def test_estimate_returns_per_channel_dict_with_provenance(monkeypatch):
    multi, geom = _multi_and_geometry()
    locs = {'g': _locs(4, 6), 'r': _locs(3, 7), 'b': _locs(100, 100)}
    captured = {}

    class FakeEstimator:
        def __init__(self, **kwargs):
            captured['uncertainty'] = kwargs['uncertainty']

        def estimate(self, merged, n_frames):
            captured['merged'] = merged
            return DriftTrace(
                dx=np.arange(n_frames, dtype=float),
                dy=2 * np.arange(n_frames, dtype=float),
            )

    monkeypatch.setattr(drift_module, 'CloudDriftEstimator', FakeEstimator)
    drift = estimate_drift_multi(
        multi,
        [locs[channel] for channel in multi.channels],
        geom=geom,
        drift_reference='r',
        estimate_channels=('g', 'r'),
        verbose=False,
    )

    assert isinstance(drift, dict)
    assert list(drift) == multi.channels

    merged = captured['merged']
    np.testing.assert_allclose(merged.x[:3], locs['r'].x)
    np.testing.assert_allclose(merged.y[:3], locs['r'].y)
    np.testing.assert_allclose(merged.x[3:], locs['g'].x / 2)
    np.testing.assert_allclose(merged.y[3:], locs['g'].y * 2)

    np.testing.assert_allclose(drift['r'].dx, [0, 1, 2])
    np.testing.assert_allclose(drift['r'].dy, [0, 2, 4])
    np.testing.assert_allclose(drift['g'].dx, [0, 2, 4])
    np.testing.assert_allclose(drift['g'].dy, [0, 1, 2])
    np.testing.assert_allclose(drift['b'].dx, [0, 0.5, 1])
    np.testing.assert_allclose(drift['b'].dy, [0, 0.5, 1])

    meta = drift['b'].metadata
    assert meta['drift_reference'] == 'r'
    assert meta['coord_reference'] == 'g'
    assert meta['estimate_channels'] == ['g', 'r']
    assert meta['channel'] == 'b'

    # lp is a CRLB in pixels; only the channel->reference pixel scale applies.
    expected_uncertainty = np.median(
        [10.0] * 3 + [10.0 * np.sqrt((0.5 ** 2 + 2.0 ** 2) / 2)] * 3)
    assert captured['uncertainty'] == pytest.approx(expected_uncertainty)


def test_correct_multi_uses_explicit_asymmetric_window():
    multi, geom = _multi_and_geometry()
    locs = {'g': _locs(30, 30), 'r': _locs(30, 30), 'b': _locs(30, 30)}
    drift = _traces()
    window = drift_crop_window(drift, (64, 64), pad=2, extra_trim=3)

    corrected, corrected_locs = correct_multi(
        multi, locs, drift, window=window, geom=geom, edge=0, verbose=False)

    assert corrected.shape == (3, 3, 50, 50)
    np.testing.assert_array_equal(corrected.offsets, np.full((3, 2), 5))
    assert set(corrected_locs) == set(multi.channels)
    np.testing.assert_allclose(corrected_locs['b'].x, [25, 24.5, 24])
    np.testing.assert_allclose(corrected_locs['b'].y, [25, 24.5, 24])
    np.testing.assert_allclose(corrected_locs['g'].x, [25, 23, 21])
    np.testing.assert_allclose(corrected_locs['g'].y, [25, 24, 23])

    windows = corrected.geometry(check_coverage=False).windows
    assert windows['g'][:2] == (5, 5)
    assert windows['b'][2:] == (50, 50)
    assert corrected.metadata['drift_crop']['window'] == [5, 5, 50, 50]
    assert corrected.metadata['drift_crop']['drift_reference'] is None
    assert corrected.ref_name == 'g'


def test_correct_multi_rejects_window_smaller_than_drift_needs():
    multi, geom = _multi_and_geometry()
    locs = {'g': _locs(30, 30), 'r': _locs(30, 30), 'b': _locs(30, 30)}

    with pytest.raises(ValueError, match='too large for the drift'):
        correct_multi(multi, locs, _traces(), window=(0, 0, 64, 64),
                      geom=geom, edge=0, verbose=False)
    with pytest.raises(ValueError, match='does not fit'):
        correct_multi(multi, locs, _traces(), window=(0, 0, 80, 80),
                      geom=geom, edge=0, verbose=False)


def test_correct_multi_requires_a_window():
    multi, geom = _multi_and_geometry()
    locs = {'g': _locs(30, 30), 'r': _locs(30, 30), 'b': _locs(30, 30)}

    with pytest.raises(TypeError):
        correct_multi(multi, locs, _traces(), geom=geom, verbose=False)


def test_plot_drift_fov_marks_common_area_and_margins():
    import matplotlib.pyplot as plt

    from smlm.multiplex import plot_drift_fov

    multi, _ = _multi_and_geometry()
    drift = _traces()

    fig, axes = plot_drift_fov(multi, drift)

    assert len(axes) == len(drift)
    assert len(axes[0].patches) == 1          # only the common-area outline
    labels = {text.get_text() for text in axes[0].texts}
    assert {'top 0 px', 'bottom 2 px', 'left 0 px', 'right 4 px'} <= labels
    assert 'common area' in axes[0].get_title()
    plt.close(fig)


def test_plot_drift_fov_transpose_keeps_labels_on_original_sides():
    import matplotlib.pyplot as plt

    from smlm.multiplex import plot_drift_fov

    multi, _ = _multi_and_geometry()
    fig, axes = plot_drift_fov(multi, _traces(), transpose=True)

    assert len(axes[0].patches) == 1
    labels = {text.get_text() for text in axes[0].texts}
    assert {'top 0 px', 'bottom 2 px', 'left 0 px', 'right 4 px'} <= labels
    # The stacked image is now wider than tall.
    x0, x1 = axes[0].get_xlim()
    y0, y1 = axes[0].get_ylim()
    assert (x1 - x0) > (y1 - y0)
    plt.close(fig)


def test_plot_multi_drift_shares_one_time_axis():
    import matplotlib.pyplot as plt

    from smlm.multiplex import plot_multi_drift

    fig, axes = plot_multi_drift(_traces(), ex_time=0.1)

    assert len(axes) == 2
    assert axes[0].get_shared_x_axes().joined(axes[0], axes[1])
    assert len(axes[0].lines) == len(_traces())
    assert len(axes[1].lines) == len(_traces())
    assert axes[1].get_xlabel() == 'time (s)'
    np.testing.assert_allclose(axes[0].lines[0].get_xdata(), [0.0, 0.1, 0.2])
    plt.close(fig)


def test_plot_multi_drift_requires_equal_frame_counts():
    from smlm.multiplex import plot_multi_drift

    drift = _traces()
    drift['g'] = DriftTrace(dx=np.zeros(5), dy=np.zeros(5))
    with pytest.raises(ValueError, match='different frame counts'):
        plot_multi_drift(drift)


def test_plot_drift_fov_requires_a_stack():
    from smlm.multiplex import plot_drift_fov

    with pytest.raises(TypeError, match='MultiImageStack'):
        plot_drift_fov((64, 64), _traces())


def test_plot_multichannel_crop_can_show_frame_zero():
    import matplotlib.pyplot as plt

    from smlm.multiplex import plot_multichannel_crop

    multi, _ = _multi_and_geometry(height=8, width=8)
    multi = multi.with_data(
        np.stack([np.full((3, 8, 8), level, dtype=np.float32)
                  for level in (0.0, 1.0, 2.0)]))

    fig = plot_multichannel_crop(multi, y=(2, 6), x=(1, 5), frame=0)

    assert 'frame 0' in fig.axes[0].get_title()
    assert fig.axes[0].get_title().endswith('(4x4)')
    plt.close(fig)


def test_unselected_output_channel_still_requires_registration(monkeypatch):
    multi, geom = _multi_and_geometry()
    incomplete_geom = DualViewGeometry(
        windows=geom.windows,
        calib={'r': geom.calibs['r']},
        calib_windows=geom.calib_windows,
        reference='g',
        check_coverage=False,
    )
    locs = {'g': _locs(30, 30), 'r': _locs(30, 30), 'b': _locs(30, 30)}

    class FakeEstimator:
        def __init__(self, **kwargs):
            pass

        def estimate(self, merged, n_frames):
            return DriftTrace.zeros(n_frames)

    monkeypatch.setattr(drift_module, 'CloudDriftEstimator', FakeEstimator)
    with pytest.raises(ValueError, match="channel 'b' has no registration"):
        estimate_drift_multi(
            multi,
            locs,
            geom=incomplete_geom,
            drift_reference='g',
            estimate_channels=('g',),
            uncertainty_px=0.1,
            verbose=False,
        )
