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
        if key == 'rejected_png':
            continue          # 只有 include_rejected=True 才写
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


def _event(eid=0):
    return fr.StepEvent(emitter_id=eid, frame=J, lag=LAG,
                        IA_pre=300.0, IA_post=0.0, ID_pre=5.0, ID_post=160.0,
                        dIA=300.0, dID=155.0, dID_long=155.0)


def test_fret_csv_carries_acceptor_state_columns(tmp_path):
    """fret.csv 带上判据量, 下游不必回头再筛一遍 (只取 red_label=='steady')。"""
    td, ta = _channels()
    recs = _records()
    for k, (label, duty, spread) in ((0, ('blinking', 0.75, 1.4)),
                                     (2, ('steady', 1.0, 0.9))):
        recs[k].event = _event(k)         # 有事件才会把 qc 写进去
        recs[k].red_label, recs[k].red_duty = label, duty
        recs[k].spread, recs[k].n_d_pre = spread, 1

    fr.save_bundle(tmp_path, records=recs, results=_results(),
                   traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX,
                   verbose=False)
    with open(tmp_path / fr.FILES['fret_csv'], encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    assert [r['red_label'] for r in rows] == ['blinking', 'steady']
    assert float(rows[0]['red_duty']) == pytest.approx(0.75)
    assert float(rows[0]['spread']) == pytest.approx(1.4)
    assert rows[1]['n_d_pre'] == '1'
    assert float(rows[0]['gamma']) == pytest.approx(1.1)     # 列序没被打乱


def test_fret_csv_qc_columns_are_blank_without_diagnostics(tmp_path):
    """没有 PairRecord 诊断时列还在, 只是空 —— 下游按名字读不会时有时无。"""
    td, ta = _channels()
    fr.save_bundle(tmp_path, records=_records(), results=_results(),
                   traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX,
                   verbose=False)
    with open(tmp_path / fr.FILES['fret_csv'], encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    for c in fr.io.QC_COLS:
        assert all(r[c] == '' for r in rows), c


def test_save_bundle_suffix_keeps_manual_round_separate(tmp_path):
    """人工复筛要落在同一目录里, 文件名必须能和机器筛的分开。"""
    assert fr.filenames('_manual')['e_png'] == 'E_distribution_manual.png'
    assert fr.filenames('_manual')['gallery_png'] == 'good_traces_all_manual.png'
    assert fr.filenames() == fr.FILES

    td, ta = _channels()
    written = fr.save_bundle(
        tmp_path, records=_records(), results=_results(),
        traces_d=td, traces_a=ta, cal=_cal(), ex_time=EX, suffix='_manual',
        include_all=False, verbose=False)
    assert os.path.basename(written['fret_csv']) == 'fret_manual.csv'
    assert os.path.basename(written['traces_good']) == 'good_traces_manual.csv'
    assert os.path.basename(written['e_png']) == 'E_distribution_manual.png'
    assert not os.path.exists(tmp_path / 'E_distribution.png')   # 没碰机器筛的


def test_save_bundle_draws_removed_traces_when_asked(tmp_path):
    """人工删掉的 trace 要留一张图, 否则事后无从复核。"""
    td, ta = _channels()
    recs = _records()
    for k in (0, 2):
        recs[k].event = _event(k)
    recs[1].event = _event(1)          # 中间这条被人工删掉
    recs[1].reasons = ['manual']

    plain = fr.save_bundle(tmp_path / 'a', records=recs, results=_results(),
                           traces_d=td, traces_a=ta, ex_time=EX, verbose=False)
    assert 'rejected_png' not in plain
    assert not (tmp_path / 'a' / fr.FILES['rejected_png']).exists()

    written = fr.save_bundle(tmp_path / 'b', records=recs, results=_results(),
                             traces_d=td, traces_a=ta, ex_time=EX,
                             include_rejected=True, verbose=False)
    p = tmp_path / 'b' / fr.FILES['rejected_png']
    assert written['rejected_png'] == str(p)
    assert p.stat().st_size > 0


def test_gallery_picks_traces_by_emitter_id_not_position(tmp_path):
    """人工复筛只载入一部分分子, 下标 ≠ emitter_id, 画图仍要取对的那条。"""
    td, ta = _channels()
    recs = _records()
    keep = [recs[0], recs[2]]
    for r in keep:
        r.event = _event(r.emitter_id)
    # 只给 id 0 和 2 两条 trace: 列表下标 1 对应 id 2
    sub_d = [td[0], td[2]]
    sub_a = [ta[0], ta[2]]

    fig = fr.viz.plot_pair_gallery(keep, sub_a, sub_d, which='kept',
                                   n=2, ex_time=EX)
    titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    assert {t.split()[0] for t in titles} == {'#0', '#2'}


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


def test_save_bundle_can_skip_the_unfiltered_table(tmp_path):
    td, ta = _channels()
    written = fr.save_bundle(tmp_path, records=_records(), results=_results(),
                             traces_d=td, traces_a=ta, cal=_cal(),
                             ex_time=EX, include_all=False, verbose=False)
    assert 'traces_all' not in written
    assert not os.path.exists(tmp_path / fr.FILES['traces_all'])
    assert 'traces_good' in written              # 合格的仍然要留


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
