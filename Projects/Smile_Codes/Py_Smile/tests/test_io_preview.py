"""preview_movie: 内联预览与 save_movie 同源。"""
import base64
import os
import tempfile

import numpy as np
import pytest

imageio = pytest.importorskip('imageio')

from smlm.datatypes import Emitter, ImageStack, Localizations
from smlm.io import preview_movie, save_movie
from smlm.io.video import _split_trace
from smlm.trace import Trace
from smlm.viz import emitter_roi, preview_emitters, save_emitters

T, H, W = 24, 40, 40
RNG = np.random.default_rng(0)


def _setup():
    stack = RNG.normal(50.0, 5.0, (T, H, W)).astype(np.float32)
    yy, xx = np.mgrid[0:H, 0:W]
    on = np.zeros(T)
    on[4:18] = 1.0
    stack += (600 * on[:, None, None]
              * np.exp(-((yy - 20) ** 2 + (xx - 22) ** 2) / (2 * 1.5 ** 2))
              ).astype(np.float32)
    image = ImageStack(stack, pixel_size=100.0)
    em = Emitter(frame=np.arange(T), x=np.full(T, 22.0), y=np.full(T, 20.0),
                 loc_idx=np.arange(T))
    locs = Localizations(frame=np.arange(T), x=np.full(T, 22.0),
                         y=np.full(T, 20.0), photons=np.full(T, 600.0),
                         bg=np.zeros(T), sigma=np.full(T, 1.5),
                         lp=np.full(T, 0.2))
    trace = Trace(emitter_id=0, frames=np.arange(T),
                  intensity=np.where(on > 0, 600.0, 0.0),
                  x=22.0, y=20.0)
    return image, em, locs, trace


def _temp_previews():
    return set(f for f in os.listdir(tempfile.gettempdir())
               if f.startswith('smlm_preview_'))


def test_preview_inlines_mp4_and_leaves_no_temp_file():
    image, em, locs, trace = _setup()
    roi = emitter_roi(image, em, half=8, frames=slice(0, 12))
    before = _temp_previews()

    html = preview_movie(roi, width=360, exposure_time=0.1, speed=5,
                         overlay=locs, trace=trace, verbose=False
                         )._repr_html_()

    assert html.startswith('<video width="360"')
    assert 'data:video/mp4;base64,' in html
    payload = base64.b64decode(html.split('base64,')[1].split('"')[0])
    assert payload[4:8] == b'ftyp'                 # 真的是 mp4
    assert len(payload) > 5000
    assert _temp_previews() == before              # 临时文件已清掉


def test_preview_can_save_at_the_same_time():
    image, em, locs, trace = _setup()
    roi = emitter_roi(image, em, half=8)
    path = os.path.join(tempfile.mkdtemp(), 'em.mp4')

    preview_movie(roi, path, exposure_time=0.1, speed=5, overlay=locs,
                  trace=trace, verbose=False)

    assert os.path.isfile(path) and os.path.getsize(path) > 0
    n_frames = sum(1 for _ in imageio.get_reader(path))   # len(reader) 对 mp4 不可靠
    assert n_frames == roi.data.shape[0]
    os.unlink(path)


def test_trace_panel_changes_the_output():
    """传 trace 会多出右侧面板 —— 说明确实走的是同一套渲染。"""
    image, em, locs, trace = _setup()
    roi = emitter_roi(image, em, half=8, frames=slice(0, 12))
    d = tempfile.mkdtemp()
    a, b = os.path.join(d, 'a.mp4'), os.path.join(d, 'b.mp4')

    save_movie(roi, a, exposure_time=0.1, speed=5, overlay=locs, verbose=False)
    save_movie(roi, b, exposure_time=0.1, speed=5, overlay=locs, trace=trace,
               verbose=False)
    with imageio.get_reader(a) as ra, imageio.get_reader(b) as rb:
        wa = ra.get_meta_data()['size'][0]
        wb = rb.get_meta_data()['size'][0]
    assert wb > wa                                  # 右侧多了曲线图
    for p in (a, b):
        os.unlink(p)


def _setup_many(n=3):
    image, _, locs, _ = _setup()
    emitters = [Emitter(frame=np.arange(T), x=np.full(T, 22.0),
                        y=np.full(T, 20.0), loc_idx=np.arange(T))
                for _ in range(n)]
    traces = [Trace(emitter_id=i, frames=np.arange(T),
                    intensity=np.full(T, 600.0)) for i in range(n)]
    return image, emitters, locs, traces


def test_preview_emitters_builds_one_video_per_index():
    image, emitters, locs, traces = _setup_many(3)

    html = preview_emitters(image, traces, emitters, locs, indices=[0, 2],
                            half=6, frames=slice(0, 8), width=200,
                            exposure_time=0.1, speed=5)._repr_html_()

    assert html.count('<video') == 2
    assert '#0' in html and '#2' in html and '#1' not in html


def test_preview_emitters_tolerates_missing_trace_and_locs():
    image, emitters, _, _ = _setup_many(2)

    html = preview_emitters(image, None, emitters, None, indices=1, half=6,
                            frames=slice(0, 6), width=200,
                            exposure_time=0.1, speed=5)._repr_html_()

    assert html.count('<video') == 1 and '#1' in html


def test_save_emitters_writes_one_file_per_index(tmp_path):
    image, emitters, locs, traces = _setup_many(3)

    paths = save_emitters(image, traces, emitters, locs, indices=[0, 2],
                          out_dir=str(tmp_path),
                          pattern='{tag}_em{idx:03d}.mp4', tag='toy',
                          half=6, frames=slice(0, 8), exposure_time=0.1,
                          speed=5, verbose=False)

    assert [os.path.basename(p) for p in paths] == ['toy_em000.mp4',
                                                    'toy_em002.mp4']
    assert all(os.path.getsize(p) > 0 for p in paths)


def test_split_trace_accepts_dict_and_trace():
    frames = np.arange(5)
    values = np.linspace(0, 1, 5)
    f1, v1 = _split_trace({'frames': frames, 'intensity': values})
    f2, v2 = _split_trace(Trace(emitter_id=0, frames=frames, intensity=values))
    np.testing.assert_array_equal(f1, frames)
    np.testing.assert_array_equal(v1, values)
    np.testing.assert_array_equal(f2, frames)
    np.testing.assert_array_equal(v2, values)
