import numpy as np
import pytest

from smlm.multiplex import pair_multichannels, pair_channels


def test_pair_multichannels_aligns_each_pair_to_selected_reference():
    reference = np.array([
        [10.0, 10.0],
        [30.0, 20.0],
        [50.0, 40.0],
        [70.0, 55.0],
    ])
    shift_g = np.array([3.0, -2.0])
    shift_b = np.array([-4.0, 5.0])
    locs = [reference - shift_g, reference, reference - shift_b]

    pairs = pair_multichannels(
        locs,
        channels=('g', 'r', 'b'),
        reference='r',
        tols=(0.5,),
        min_pairs=4,
    )

    assert len(pairs) == 3
    assert len(pairs[0]) == len(pairs[2]) == 4
    assert pairs[1] is None
    np.testing.assert_allclose(pairs[0].shift, shift_g)
    np.testing.assert_allclose(pairs[2].shift, shift_b)
    np.testing.assert_array_equal(pairs[0].ref, reference)
    np.testing.assert_array_equal(pairs[2].ref, reference)


def test_pair_multichannels_defaults_to_first_channel():
    reference = np.array([[10.0, 10.0], [30.0, 20.0], [50.0, 40.0]])
    moving = reference - np.array([2.0, -1.0])

    pairs = pair_multichannels([reference, moving], tols=(0.5,), min_pairs=3)

    assert pairs[0] is None
    assert len(pairs[1]) == 3


def test_pair_multichannels_validates_channel_labels_and_reference():
    locs = [np.zeros((1, 2)), np.ones((1, 2))]

    with pytest.raises(ValueError, match='channels length must match locs length'):
        pair_multichannels(locs, channels=('g',))
    with pytest.raises(ValueError, match='unknown reference channel'):
        pair_multichannels(locs, channels=('g', 'r'), reference='b')


def test_pair_channels_retains_shift_when_below_minimum_pair_count():
    ref = np.array([[10.0, 10.0], [30.0, 20.0]])
    mov = ref - np.array([2.0, -1.0])

    pairs, shift = pair_channels(mov, ref, tols=(0.5,), min_pairs=3)

    assert len(pairs) == 0
    np.testing.assert_allclose(pairs.shift, shift)
