import matplotlib
matplotlib.use('Agg')

import numpy as np
import pytest

from smlm.datatypes import ImageStack, Localizations, MultiImageStack
from smlm.multiplex import (
    DualViewCalibration,
    MultichannelROI,
    plot_multichannel_crop,
    plot_multichannel_localizations,
)


def test_multichannel_roi_crops_equal_windows_and_preserves_metadata():
    data = np.arange(3 * 6 * 8).reshape(3, 6, 8)
    image = ImageStack(data=data, pixel_size=100, metadata={'source': 'test'})
    roi = MultichannelROI(
        y=(1, 4), x=(2, 5), offsets=((0, 0), (1, 2)),
        channels=('ref', 'mov'))

    multi = roi(image)

    assert multi.shape == (2, 3, 3, 3)
    np.testing.assert_array_equal(multi.data[0], data[:, 1:4, 2:5])
    np.testing.assert_array_equal(multi.data[1], data[:, 2:5, 4:7])
    np.testing.assert_array_equal(multi.offsets, [[1, 2], [2, 4]])
    assert multi.channels == ['ref', 'mov']
    assert multi.ref_name == 'ref'
    assert multi.pixel_size == 100
    assert multi.metadata['source'] == 'test'


def test_multichannel_roi_rejects_out_of_bounds_shift():
    image = ImageStack(data=np.zeros((2, 5, 6)))
    roi = MultichannelROI(y=(2, 5), x=(1, 5), offsets=((0, 0), (1, 0)))

    with pytest.raises(ValueError, match='outside image'):
        roi(image)


def test_multichannel_roi_requires_zero_reference_offset():
    with pytest.raises(ValueError, match=r'offsets\[0\].*\(0, 0\)'):
        MultichannelROI(y=(0, 2), x=(0, 2), offsets=((1, 0),))


def test_show_roi_returns_figure_with_each_window():
    image = ImageStack(data=np.zeros((2, 6, 8)))
    roi = MultichannelROI(
        y=(1, 3), x=(2, 5), offsets=((0, 0), (1, 2)))

    fig = roi.show_roi(image)

    assert len(fig.axes[0].patches) == 2
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_multichannel_stack_crop_updates_each_channel_origin_and_reg():
    calibration = DualViewCalibration(
        poly_a=np.array([0.0, 1.0, 0.0]),
        poly_b=np.array([0.0, 0.0, 1.0]),
        center=np.array([0.0, 0.0]),
        scale=1.0,
        affine=np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]),
        poly_order=1,
        roi={'type': 'MultichannelROI', 'windows': {
            'g': [10, 20, 8, 9],
            'r': [10, 276, 8, 9],
        }},
        reference='g',
        moving='r',
    )
    data = np.arange(2 * 3 * 8 * 9).reshape(2, 3, 8, 9)
    multi = MultiImageStack(
        data=data,
        channels=['g', 'r'],
        reference=0,
        reg={'r': calibration},
        offsets=np.array([[10, 20], [10, 276]]),
        metadata={'source': 'test'},
    )

    cropped = multi.crop(y=(2, 7), x=(3, 8))

    assert cropped.shape == (2, 3, 5, 5)
    np.testing.assert_array_equal(cropped.data, data[:, :, 2:7, 3:8])
    np.testing.assert_array_equal(
        cropped.offsets, [[12, 23], [12, 279]])
    assert cropped.reg['r'] is calibration
    assert cropped.metadata['analysis_crops'] == [
        {'y': [2, 7], 'x': [3, 8]}]
    np.testing.assert_allclose(
        cropped.to_reference([[1.0, 3.0]], 'r'), [[1.0, 3.0]])
    assert cropped.geometry(check_coverage=False).calib_delta('r').tolist() \
        == [3.0, 2.0]


def test_multichannel_stack_crop_rejects_invalid_window():
    multi = MultiImageStack(data=np.zeros((2, 1, 5, 6)))

    with pytest.raises(ValueError, match='outside local size'):
        multi.crop(y=(0, 5), x=(1, 7))

    with pytest.raises(ValueError, match='start, stop'):
        multi.crop(y=(0, 4, 5), x=(0, 4))


def test_plot_multichannel_crop_shows_adjacent_channels_and_each_window():
    multi = MultiImageStack(
        data=np.zeros((2, 1, 8, 6)),
        channels=['g', 'r'],
        offsets=np.array([[0, 10], [0, 266]]),
    )

    fig = plot_multichannel_crop(multi, y=(1, 7), x=(2, 5))

    assert len(fig.axes) == 1
    assert len(fig.axes[0].patches) == 2
    assert [patch.get_x() for patch in fig.axes[0].patches] == [1.5, 7.5]
    assert fig.axes[0].get_title().endswith('(6x3)')
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_multichannel_localizations_overlays_each_channel():
    multi = MultiImageStack(
        data=np.zeros((2, 1, 8, 8)),
        channels=['g', 'r'],
    )
    locs = {
        'g': Localizations(frame=np.array([0]), x=np.array([2.0]),
                           y=np.array([3.0]), photons=np.array([100.0]),
                           bg=np.array([1.0])),
        'r': Localizations(frame=np.array([0]), x=np.array([4.0]),
                           y=np.array([5.0]), photons=np.array([120.0]),
                           bg=np.array([1.0])),
    }

    fig = plot_multichannel_localizations(multi, locs)

    assert len(fig.axes) == 1
    assert len(fig.axes[0].collections) == 2
    np.testing.assert_array_equal(
        fig.axes[0].collections[0].get_paths()[0].vertices,
        [[-0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [0.5, -0.5]])
    np.testing.assert_array_equal(
        fig.axes[0].collections[0].get_offsets(), [[2.0, 3.0]])
    np.testing.assert_array_equal(
        fig.axes[0].collections[1].get_offsets(), [[12.0, 5.0]])
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_multichannel_localizations_orders_by_channel_name():
    multi = MultiImageStack(
        data=np.zeros((2, 1, 8, 8)),
        channels=['g', 'r'],
    )
    green = Localizations(frame=np.array([0]), x=np.array([2.0]),
                          y=np.array([3.0]), photons=np.array([100.0]),
                          bg=np.array([1.0]))
    red = Localizations(frame=np.array([0]), x=np.array([4.0]),
                        y=np.array([5.0]), photons=np.array([120.0]),
                        bg=np.array([1.0]))

    fig = plot_multichannel_localizations(multi, {'r': red, 'g': green})

    np.testing.assert_array_equal(
        fig.axes[0].collections[0].get_offsets(), [[2.0, 3.0]])
    np.testing.assert_array_equal(
        fig.axes[0].collections[1].get_offsets(), [[12.0, 5.0]])
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_plot_multichannel_localizations_requires_channel_mapping():
    multi = MultiImageStack(
        data=np.zeros((2, 1, 8, 8)),
        channels=['g', 'r'],
    )
    green = Localizations(frame=np.array([0]), x=np.array([2.0]),
                          y=np.array([3.0]), photons=np.array([100.0]),
                          bg=np.array([1.0]))

    with pytest.raises(TypeError, match='must be a .* mapping'):
        plot_multichannel_localizations(multi, [green, green])
    with pytest.raises(ValueError, match='missing localization sets'):
        plot_multichannel_localizations(multi, {'g': green})
