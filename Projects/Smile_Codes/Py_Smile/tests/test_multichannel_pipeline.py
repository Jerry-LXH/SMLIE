import importlib

import numpy as np

from smlm.datatypes import ImageStack
from smlm.multiplex import (
    DualViewCalibration,
    MultiChannelCalibrationBundle,
    MultichannelROI,
    PairSet,
    fit_multichannel_pairs,
    pair_multichannel_files,
    multichannel_stack_from_calibration,
)


pipeline = importlib.import_module('smlm.multiplex.pipeline')


def test_pair_multichannel_files_merges_pairs_per_channel_and_file(monkeypatch):
    points = np.array([
        [0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0],
        [0.0, 2.0], [2.0, 0.0], [1.0, 3.0], [3.0, 1.0],
        [0.0, 3.0], [3.0, 0.0],
    ])
    image = ImageStack(data=np.zeros((2, 4, 8)))

    class FakeCamera:
        def __init__(self, **kwargs):
            pass

        def process(self, raw):
            return image

    class FakeDetector:
        def __init__(self, k_sigma, **kwargs):
            self.points = points

        def detect(self, stack):
            return self.points

    class FakeLocalizer:
        def __init__(self, **kwargs):
            pass

        def localize(self, stack, candidates):
            return candidates

    monkeypatch.setattr(pipeline, 'load_raw', lambda path: path)
    monkeypatch.setattr(pipeline, 'EMCCD', FakeCamera)
    monkeypatch.setattr(pipeline, 'LocalMaxDetector', FakeDetector)
    monkeypatch.setattr(pipeline, 'MLELocalizer', FakeLocalizer)

    roi = MultichannelROI(
        y=(0, 4), x=(0, 4),
        offsets=((0, 0), (0, 4)),
        channels=('g', 'r'),
    )
    pairs = pair_multichannel_files(
        ['one.sif', 'two.sif'], roi, verbose=False)

    assert len(pairs) == 2
    assert pairs[0] is None
    assert len(pairs[1]) == 20
    assert pairs[1].files() == ['one.sif', 'two.sif']

    reverse_pairs = pair_multichannel_files(
        ['one.sif', 'two.sif'], roi, reference='r', verbose=False)

    assert len(reverse_pairs[0]) == 20
    assert reverse_pairs[1] is None
    assert reverse_pairs[0].files() == ['one.sif', 'two.sif']


def test_fit_multichannel_pairs_uses_reference_slot_and_roi(monkeypatch):
    calls = []

    def fake_fit_transform(pairs, **kwargs):
        result = {'pairs': pairs, **kwargs}
        calls.append(result)
        return result

    monkeypatch.setattr(pipeline, 'fit_transform', fake_fit_transform)
    roi = MultichannelROI(
        y=(10, 20), x=(30, 40),
        offsets=((0, 0), (2, 5), (4, 9)),
        channels=('g', 'r', 'b'),
    )
    pairs_g = PairSet(
        [[1, 2], [3, 4]], [[1.1, 2.2], [3.1, 4.2]],
        file=['calib.sif', 'calib.sif'])
    pairs_b = PairSet(
        [[5, 6], [7, 8]], [[5.1, 6.2], [7.1, 8.2]],
        file=['calib.sif', 'calib.sif'])

    calibrations = fit_multichannel_pairs(
        [pairs_g, None, pairs_b], roi, verbose=False)

    assert calibrations[1] is None
    assert calibrations[0] is calls[0]
    assert calibrations[2] is calls[1]
    assert [call['reference'] for call in calls] == ['r', 'r']
    assert [call['moving'] for call in calls] == ['g', 'b']
    assert [call['coarse_offset'] for call in calls] == [(-2, -5), (2, 4)]
    assert calls[0]['roi']['windows'] == roi.windows
    assert calls[0]['roi']['files'] == ['calib.sif']


def test_multichannel_stack_from_calibration_attaches_reg_and_roi(tmp_path):
    roi = MultichannelROI(
        y=(0, 4), x=(0, 4),
        offsets=((0, 0), (0, 4)),
        channels=('g', 'r'),
    )
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
        roi={
            'type': 'MultichannelROI',
            'windows': {
                'g': [0, 0, 4, 4],
                'r': [0, 4, 4, 4],
            },
        },
        reference='g',
        moving='r',
    )
    bundle = MultiChannelCalibrationBundle.from_fit(
        [None, calibration], roi, ['calibration.sif'])
    calibration_path = tmp_path / 'calibration.npz'
    bundle.save(calibration_path)
    image = ImageStack(data=np.arange(64).reshape(2, 4, 8), pixel_size=106.6)

    multi = multichannel_stack_from_calibration(image, calibration_path)

    assert multi.channels == ['g', 'r']
    assert multi.ref_name == 'g'
    assert set(multi.reg) == {'r'}
    assert multi.channel_offset('r') == (0, 4)
    assert multi.pixel_size == 106.6
    assert multi.metadata['calibration_bundle'] == str(calibration_path)
    assert multi.metadata['roi']['windows'] == {
        'g': [0, 0, 4, 4],
        'r': [0, 4, 4, 4],
    }
    np.testing.assert_array_equal(multi.data[0], image.data[:, :, :4])
    np.testing.assert_array_equal(multi.data[1], image.data[:, :, 4:])
    np.testing.assert_array_equal(
        multi.to_reference([[1.0, 2.0]], 'r'), [[1.0, 2.0]])
    assert multi.geometry(check_coverage=False).calib_windows == {
        'g': (0, 0, 4, 4),
        'r': (0, 4, 4, 4),
    }
