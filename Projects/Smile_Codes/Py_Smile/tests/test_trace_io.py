"""池化 npz 的读取与"无影像"追加字段。"""
import numpy as np
import pytest

import smlm.trace as tr


def _write_pooled(path, n_frames, T=None, *, with_n_frames=True):
    """写一个和池化输出同格式的 npz (短的补 NaN)。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    T = int(T or max(n_frames))
    I = np.full((len(n_frames), T), np.nan, np.float32)
    B = np.full((len(n_frames), T), np.nan, np.float32)
    for i, n in enumerate(n_frames):
        I[i, :n] = np.linspace(10, 100, n)
        B[i, :n] = 5.0
    kw = dict(intensity=I, bg_flux=B, emitter_id=np.arange(len(n_frames)) + 100,
              file_id=np.zeros(len(n_frames), int),
              files=np.array(['fovA', 'fovB']),
              t_s=np.arange(T) * 0.1, exposure_s=0.1)
    if with_n_frames:
        kw['n_frames'] = np.asarray(n_frames, int)
    np.savez_compressed(path, **kw)
    return path


def test_load_pooled_trims_nan_padding(tmp_path):
    p = _write_pooled(tmp_path / 'G' / 'traces_pooled.npz', [4, 6, 5])
    traces, meta = tr.load_pooled(p)

    assert [len(t['frames']) for t in traces] == [4, 6, 5]
    assert [t['intensity'].size for t in traces] == [4, 6, 5]
    for t in traces:
        assert not np.isnan(t['intensity']).any()
        assert not np.isnan(t['bg_flux']).any()
        assert t['t_s'].size == t['intensity'].size
    assert meta['exposure_s'] == pytest.approx(0.1)
    assert np.array_equal(meta['n_frames'], [4, 6, 5])
    assert meta['n_columns'] == 6
    assert traces[0]['emitter_id'] == 100
    assert traces[0]['channel'] == 'G'          # 默认取父目录名


def test_load_pooled_falls_back_to_last_finite(tmp_path):
    p = _write_pooled(tmp_path / 'G' / 'traces_pooled.npz', [4, 6, 5],
                      with_n_frames=False)
    traces, _ = tr.load_pooled(p, channel='SPARK_6mW')
    assert [len(t['frames']) for t in traces] == [4, 6, 5]
    assert traces[0]['channel'] == 'SPARK_6mW'


def test_changepoint_on_untrimmed_pooled_data_is_silently_wrong(tmp_path):
    """补 NaN 的 trace 直接喂给 segment 会静默退化成一个暗段 —— 这是截断的理由。"""
    n = 300
    y = np.concatenate([np.full(n // 2, 5.0), np.full(n - n // 2, 500.0)])
    padded = np.concatenate([y, np.full(100, np.nan)])
    good, _ = tr.changepoint.segment(y, min_seg=5)
    bad, bounds = tr.changepoint.segment(padded, min_seg=5)
    assert np.unique(good).size == 2            # 截断后能分出亮/暗
    assert np.unique(bad).size == 1            # 带 NaN 只剩一级
    assert bounds.size == 2


def _step_trace(n=400, on_at=100, off_at=200):
    y = np.full(n, 5.0)
    y[on_at:off_at] = 500.0
    return tr.Trace(emitter_id=0, frames=np.arange(n), x=1.0, y=1.0,
                    channel='t', intensity=y)


def test_apply_fillers_without_image():
    traces = [_step_trace(), _step_trace(on_at=50, off_at=120)]
    out = tr.apply_fillers(traces, [tr.ChangepointStates(bleach_frames=100)])

    assert out == traces                       # 同一批对象, 原地修改
    for t in traces:
        assert 'state' in t and 'noise' in t
        assert t['n_on'] >= 1
    assert tr.select_traces(traces, n_on=(1, None)) == [0, 1]


def test_apply_fillers_rejects_nan():
    t = _step_trace()
    t['intensity'] = t['intensity'].copy()
    t['intensity'][300:] = np.nan
    with pytest.raises(ValueError, match='NaN'):
        tr.apply_fillers([t], [tr.ThresholdStates(high=50.0, bleach_frames=50)])


def test_apply_fillers_rejects_filler_that_needs_image():
    t = _step_trace()
    with pytest.raises(ValueError, match='影像块'):
        tr.apply_fillers([t], [tr.ApertureRaw()])
