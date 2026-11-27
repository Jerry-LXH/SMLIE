import zipfile

import numpy as np
import pytest

from smlm.datatypes import ImageStack, Localizations
from smlm.io import result_path, run_cached
from smlm.localize import pipeline as localize_pipeline
from smlm.multiplex import pipeline as multiplex_pipeline


def _locs():
    return Localizations(
        frame=np.array([0, 1]),
        x=np.array([2.0, 3.0]),
        y=np.array([4.0, 5.0]),
        photons=np.array([20.0, 30.0]),
        bg=np.array([1.0, 1.0]),
        sigma=np.array([1.0, 1.0]),
    )


def _stack(crop_offset=(4, 5)):
    return ImageStack(
        data=np.zeros((3, 8, 9), dtype=np.float32),
        pixel_size=106.6,
        metadata={
            'data_file': 'sample.sif',
            'channel': 'g',
            'crop_offset': crop_offset,
            'analysis_crops': [{'y': [1, 7], 'x': [2, 8]}],
            'calibration_bundle': 'calibration.npz',
        },
    )


def test_localize_stack_reads_or_overwrites_explicit_cache(
        tmp_path, monkeypatch, capsys):
    pipeline_runs = []
    filter_calls = []

    class FakeDetector:
        def __init__(self, k_sigma, **kwargs):
            self.k_sigma = k_sigma
            self.kwargs = kwargs

    class FakeLocalizer:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

    class FakePipeline:
        def __init__(self, detector, localizer):
            self.detector = detector
            self.localizer = localizer

        def run(self, stack):
            pipeline_runs.append(stack)
            return _locs()

    monkeypatch.setattr(localize_pipeline, 'LocalMaxDetector', FakeDetector)
    monkeypatch.setattr(localize_pipeline, 'MLELocalizer', FakeLocalizer)
    monkeypatch.setattr(localize_pipeline, 'ClassicalPipeline', FakePipeline)

    def fake_filter(locs, **kwargs):
        filter_calls.append(kwargs)
        return locs[np.array([True, False])]

    monkeypatch.setattr(localize_pipeline, 'filter_locs', fake_filter)

    cache_dir = tmp_path / 'loc_cache'
    cache_dir.mkdir()
    cache_base = cache_dir / 'sample_g'
    stack = _stack()
    filt_kw = {'photons': (15, 25)}

    computed = localize_pipeline.localize_stack(
        stack, filt_kw=filt_kw, src_path=cache_base,
        read_cache=False, verbose=False)
    cache_file = result_path(cache_base)
    assert len(pipeline_runs) == 1
    assert len(computed) == 1
    assert len(filter_calls) == 1
    assert cache_file.endswith('.locs.npz')
    with zipfile.ZipFile(cache_file) as archive:
        assert all(item.compress_type == zipfile.ZIP_STORED
                   for item in archive.infolist())

    loaded = localize_pipeline.localize_stack(
        stack, k_sigma=2.0, filt_kw={'photons': (20, 30)},
        src_path=cache_base, read_cache=True, verbose=True)
    assert len(pipeline_runs) == 1
    assert len(computed) == len(loaded) == 1
    assert len(filter_calls) == 1
    assert '[cache] loaded 1 filtered locs' in capsys.readouterr().out
    np.testing.assert_array_equal(loaded.x, computed.x)

    _, meta = run_cached(stack, cache_base, read_cache=True)
    assert meta['meta']['filter_config'] == {'photons': [15, 25]}

    changed_metadata = _stack(crop_offset=(0, 0))
    metadata_locs = localize_pipeline.localize_stack(
        changed_metadata, src_path=cache_base, read_cache=True)
    np.testing.assert_array_equal(metadata_locs.x, computed.x)

    changed_pixels = np.zeros_like(stack.data)
    changed_pixels[0, 0, 0] = 1
    changed_stack = ImageStack(
        data=changed_pixels, pixel_size=stack.pixel_size,
        metadata=stack.metadata)
    with pytest.raises(ValueError, match='input image shape, dtype, or first frame'):
        localize_pipeline.localize_stack(
            changed_stack, src_path=cache_base, read_cache=True)

    changed_dtype = ImageStack(
        data=np.asarray(stack.data, dtype=np.float64),
        pixel_size=stack.pixel_size,
        metadata=stack.metadata,
    )
    with pytest.raises(ValueError, match='input image shape, dtype, or first frame'):
        localize_pipeline.localize_stack(
            changed_dtype, src_path=cache_base, read_cache=True)

    different_shape = ImageStack(
        data=np.zeros((4, 8, 9), dtype=np.float32),
        pixel_size=stack.pixel_size,
        metadata=stack.metadata,
    )
    with pytest.raises(ValueError, match='input image shape, dtype, or first frame'):
        localize_pipeline.localize_stack(
            different_shape, src_path=cache_base, read_cache=True)

    pipeline_runs_before_no_read = len(pipeline_runs)
    assert run_cached(stack, cache_dir / 'not_created',
                      read_cache=False) == (None, None)
    assert len(pipeline_runs) == pipeline_runs_before_no_read

    localize_pipeline.localize_stack(
        stack, filt_kw=filt_kw, src_path=cache_base,
        read_cache=False, verbose=False)
    assert len(pipeline_runs) == 2
    assert len(filter_calls) == 2


def test_run_cached_read_requires_existing_path():
    with pytest.raises(ValueError, match='src_path is required'):
        run_cached(_stack(), None, read_cache=True)

    with pytest.raises(FileNotFoundError, match='cache not found'):
        run_cached(_stack(), 'missing_cache', read_cache=True)

