"""落盘: save_bundle 的目录结构 / 表内容 / trace 往返。"""
import csv
import os

import numpy as np
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

import smlm.fret as fr                                  # noqa: E402
from smlm.datatypes import Trace                        # noqa: E402

N = 200
J, LAG = 100, 2
EX = 0.1


def _trace(eid, intensity, x=10.0, y=20.0):
    v = np.asarray(intensity, float)
    n = v.size
    flux = np.full(n, 7.0)                              # 背景环扣掉的光, 每帧常数
    return Trace(emitter_id=eid, frames=np.arange(n), x=float(x), y=float(y),
                 channel='t', intensity=v, raw=v + flux, bg=np.zeros(n),
                 bg_std=np.zeros(n), bg_flux=flux, sw=5.0, sw2=8.0, area=50.0,
                 radius=2.0, neighbours=0, sigma=1.2)


def _channels(n_mol=3):
    """供体(绿)在台阶后亮起; 受体(红)在台阶后灭掉。"""
    rng = np.random.default_rng(0)
    td, ta = [], []
    for i in range(n_mol):
        a = np.full(N, 300.0)
        a[J:] = 0.0
        d = np.full(N, 5.0)
        d[J + LAG:] = 160.0
        td.append(_trace(i, d + rng.normal(0, 3, N)))
        ta.append(_trace(i, a + rng.normal(0, 5, N)))
    return td, ta


def _records(keeps=(True, False, True)):
    return [fr.PairRecord(emitter_id=i, keep=k,
                          reasons=[] if k else ['red'])
            for i, k in enumerate(keeps)]


def _results(ids=(0, 2)):
    return [fr.FretResult(emitter_id=i, frame=J, E=0.9 - 0.05 * k,
                          E_don=0.88, IA_pre=300.0, ID_pre=5.0)
            for k, i in enumerate(ids)]


def _cal():
    return fr.FretCalibration(gamma=1.1, alpha=0.09, beta=-2.0, detail={})


def test_save_bundle_writes_every_file(tmp_path):
    td, ta = _channels()
    out = tmp_path / 'fov1_fret'
    written = fr.save_bundle(out, records=_records(), results=_results(),
                             traces_d=td, traces_a=ta, cal=_cal(),
                             ex_time=EX, title='fov1', verbose=False)

    assert set(written) == {'dir', 'traces_all', 'traces_good', 'fret_csv',
                            'gallery_png', 'e_png'}
    for key, name in fr.FILES.items():
        p = os.path.join(out, name)
        assert os.path.getsize(p) > 0, name
        assert written[key] == p


def test_fret_csv_keeps_gamma_alpha_beta_on_every_row(tmp_path):
    td, ta = _channels()
    fr.save_bundle(tmp_path, records=_records(), results=_results(),
                   traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX,
                   verbose=False)
    with open(tmp_path / fr.FILES['fret_csv'], encoding='utf-8-sig') as f:
        rows = list(csv.reader(f))
    assert rows[0] == fr.io.COLS
    assert len(rows) == 1 + 2                       # 表头 + 两个分子
    body = rows[1:]
    assert [r[0] for r in body] == ['0', '2']
    assert all(float(r[-3]) == 1.1 for r in body)   # gamma
    assert all(float(r[-2]) == 0.09 for r in body)  # alpha


def test_trace_tables_hold_all_and_only_kept(tmp_path):
    td, ta = _channels()
    fr.save_bundle(tmp_path, records=_records(), results=_results(),
                   traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX,
                   verbose=False)
    allt = fr.load_traces_csv(tmp_path / fr.FILES['traces_all'])
    good = fr.load_traces_csv(tmp_path / fr.FILES['traces_good'])
    assert list(allt['emitter_id']) == [0, 1, 2]
    assert list(good['emitter_id']) == [0, 2]
    assert allt['I_acceptor'].shape == (3, N)
    # 时间轴用 ex_time
    assert allt['t_s'][0, 1] == pytest.approx(EX)


def test_traces_roundtrip_matches_the_source(tmp_path):
    td, ta = _channels()
    fr.save_bundle(tmp_path, records=_records(), results=_results(),
                   traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX,
                   verbose=False)
    got = fr.load_traces_csv(tmp_path / fr.FILES['traces_good'])
    for k, i in enumerate((0, 2)):
        assert np.allclose(got['I_donor'][k], td[i].intensity, atol=0.05)
        assert np.allclose(got['I_acceptor'][k], ta[i].intensity, atol=0.05)
        assert np.allclose(got['bg_donor'][k], td[i].bg_flux)
        assert got['x_px'][k, 0] == pytest.approx(td[i].x)


def test_save_bundle_without_results_still_writes_tables(tmp_path):
    td, ta = _channels()
    written = fr.save_bundle(tmp_path, records=_records(), results=[],
                             traces_d=td, traces_a=ta, cal=_cal(),
                             ex_time=EX, verbose=False)
    assert 'fret_csv' in written and 'traces_all' in written
    assert 'e_png' not in written                   # 没有 E 就不画
    with open(tmp_path / fr.FILES['fret_csv'], encoding='utf-8-sig') as f:
        assert len(list(csv.reader(f))) == 1        # 只剩表头


def test_save_bundle_without_any_kept_pair_skips_the_gallery(tmp_path):
    td, ta = _channels()
    written = fr.save_bundle(tmp_path, records=_records((False, False, False)),
                             results=[], traces_d=td, traces_a=ta,
                             cal=_cal(), ex_time=EX, verbose=False)
    assert 'gallery_png' not in written and 'traces_good' not in written
    assert 'traces_all' in written                  # 未筛选的照旧要留


def test_save_traces_csv_reads_bg_flux_from_the_trace(tmp_path):
    """回归: 曾经去读 trace.meta, 而 Trace 早已没有 meta 字段。"""
    td, ta = _channels(1)
    p = tmp_path / 't.csv'
    fr.save_traces_csv(td, ta, np.arange(N) * EX, p, verbose=False)
    got = fr.load_traces_csv(p)
    assert np.allclose(got['bg_donor'][0], 7.0)     # 不是 nan, 也不是崩
    assert np.allclose(got['bg_acceptor'][0], 7.0)


def test_load_traces_csv_rejects_ragged_tables(tmp_path):
    p = tmp_path / 'ragged.csv'
    p.write_text(','.join(fr.TRACE_COLS) + '\n'
                 + '0,0,0,1,2,10,20,7,7\n'
                 + '0,1,0.1,1,2,10,20,7,7\n'
                 + '1,0,0,1,2,10,20,7,7\n', encoding='utf-8')
    with pytest.raises(ValueError, match='帧数不一致'):
        fr.load_traces_csv(p)
