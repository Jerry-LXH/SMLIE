import numpy as np

from smlm.multiplex import (
    DualViewCalibration,
    MultiChannelCalibrationBundle,
    MultichannelROI,
)


def _calibration(moving, reference, offset):
    return DualViewCalibration(
        poly_a=np.array([1.0, 0.1, 0.2]),
        poly_b=np.array([2.0, 0.3, 0.4]),
        center=np.array([10.0, 20.0]),
        scale=5.0,
        affine=np.array([
            [1.0, 0.0],
            [0.0, 1.0],
            [0.5, -0.5],
        ]),
        poly_order=1,
        roi={'files': ['std.sif']},
        metrics={'poly': {'rms': 0.1}},
        coarse_offset=offset,
        reference=reference,
        moving=moving,
    )


def test_multichannel_bundle_round_trip(tmp_path):
    roi = MultichannelROI(
        y=(10, 30), x=(20, 50),
        offsets=((0, 0), (2, 5), (4, 9)),
        channels=('g', 'r', 'b'),
    )
    calibrations = [
        _calibration('g', 'r', (-2, -5)),
        None,
        _calibration('b', 'r', (2, 4)),
    ]
    source_files = [tmp_path / 'std_one.sif', tmp_path / 'std_two.sif']
    bundle = MultiChannelCalibrationBundle.from_fit(
        calibrations, roi, source_files)
    path = tmp_path / 'multichannel_calibration.npz'

    assert bundle.save(path) == path
    loaded = MultiChannelCalibrationBundle.load(path)

    assert loaded.channels == ('g', 'r', 'b')
    assert loaded.reference == 'r'
    assert loaded.calibrations[1] is None
    assert loaded.roi['coordinate_anchor'] == 'g'
    assert loaded.roi['reference'] == 'r'
    assert loaded.roi['windows'] == {
        'g': [10, 20, 20, 30],
        'r': [12, 25, 20, 30],
        'b': [14, 29, 20, 30],
    }
    assert loaded.source_files == [
        {'path': str(path), 'name': path.name} for path in source_files
    ]
    for original, restored in zip(calibrations, loaded.calibrations):
        if original is None:
            assert restored is None
            continue
        np.testing.assert_array_equal(restored.poly_a, original.poly_a)
        np.testing.assert_array_equal(restored.affine, original.affine)
        assert restored.metrics == original.metrics
        assert restored.roi == original.roi
        assert restored.coarse_offset == original.coarse_offset
        np.testing.assert_array_equal(
            restored.transform([[1.0, 2.0]]),
            original.transform([[1.0, 2.0]]),
        )
