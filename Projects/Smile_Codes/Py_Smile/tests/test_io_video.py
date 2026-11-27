"""_FigureRenderer 的 numpy 合成路径必须与逐帧 canvas.draw() 逐像素一致。

这里的"参考实现"就是旧版逐帧重绘; 只要它们对得上, 加速就不会改变视频内容。
"""
import numpy as np
import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from smlm.io.video import _FigureRenderer, _quantize  # noqa: E402

H = W = 25
RNG = np.random.default_rng(0)


def _renderer(cm='gray', h=H, w=W, trace=False, **kw):
    if trace:
        kw.setdefault('trace_x', np.arange(24.0))
        kw.setdefault('trace_y', np.arange(24.0) * 10.0)
    return _FigureRenderer(h, w, cm, 0.0, 3000.0, overlay_kw=None,
                           x_label='Frame', y_label='I', **kw)


def _panel_diff(r, frame, pts=None, j=0):
    fast = r(frame, pts, j).astype(int)
    ref = r._render_reference(frame, pts, j).astype(int)
    r0, r1, c0, c1 = r.panel
    return np.abs(fast - ref)[r0:r1, c0:c1]


@pytest.mark.parametrize('cm', ['gray', 'viridis'])
def test_image_panel_matches_canvas_draw_exactly(cm):
    """图像面板 (真正有数据变化的部分) 必须逐像素相同。"""
    r = _renderer(cm)
    for scale in (0.0, 0.03, 0.5, 1.0):
        frame = RNG.uniform(0.0, 1.0, (H, W)) * 3000.0 * scale
        if RNG.random() < 0.5:                 # 夹带 NaN / 越界值
            frame = frame.copy()
            frame[3, 4] = np.nan
            frame[5, 6] = 9e9
        assert _panel_diff(r, frame).max() == 0


def test_image_panel_matches_at_one_to_one_scale():
    """不走 trace 面板时画布就是 25x25 (1:1), 也必须一致。"""
    r = _renderer()
    assert r.base.shape[:2] == (H, W)
    assert _panel_diff(r, RNG.uniform(0.0, 3000.0, (H, W))).max() == 0


def test_overlay_markers_land_on_the_right_pixels():
    """定位点精灵 (alpha 混合 + 对位 + 像素吸附) 必须与 scatter 重合。

    matplotlib 把 marker 吸附到最近像素中心, 所以这里可以要求几乎逐字节一致
    (只剩 alpha 混合的 ±1 取整差)。
    """
    for kw in (None, dict(c='lime', marker='+', s=60, linewidths=1.5)):
        r = _FigureRenderer(H, W, 'gray', 0.0, 3000.0, overlay_kw=kw,
                            x_label='F', y_label='I')
        frame = np.full((H, W), 1500.0)
        # 12.5/16.5 这种"整数+半像素"是 matplotlib 吸附规则的平局点
        pts = np.array([[12.0, 8.0], [3.0, 21.0], [20.0, 12.0]])
        fast = r(frame, pts, 0).astype(int)
        ref = r._render_reference(frame, pts, 0).astype(int)
        # 位置/形状必须一致; 只剩 alpha 混合的取整差 (大 marker 边缘 ±2)
        assert np.abs(fast - ref).max() <= 2


def test_only_the_cursor_column_differs_when_trace_is_shown():
    """带 trace 时唯一的差异是光标样式 (实线 vs 虚线+抗锯齿)。"""
    r = _renderer(trace=True)
    frame = RNG.uniform(0.0, 3000.0, (H, W))
    d = _panel_diff(r, frame)
    assert d.max() == 0                        # 图像面板仍然精确

    fast = r(frame, None, 5).astype(int)
    ref = r._render_reference(frame, None, 5).astype(int)
    full = np.abs(fast - ref).max(axis=2)
    tx, tb, _, _ = r._tr
    cursor = int(round(tx * float(r.trace_x[5]) + tb))
    cols = np.arange(full.shape[1])
    off_cursor = np.abs(cols - cursor) > 2
    assert full[:, off_cursor].max() <= 1
    assert full[:, ~off_cursor].max() > 16     # 差异确实集中在光标上
    assert full.mean() < 1.0


def test_trace_labels_stay_out_of_the_image_panel():
    """y 轴标签/刻度标签不能压到图像面板上 (否则会被合成层擦掉)。"""
    r = _renderer(trace=True)
    r0, r1, c0, c1 = r.panel
    fig_w = r.fig.get_figwidth() * r.dpi
    left = min(a.get_window_extent().x0
               for a in [r.ax_tr.yaxis.label] + list(r.ax_tr.get_yticklabels())
               if a.get_text())
    assert left >= c1
    assert left < fig_w


@pytest.mark.parametrize('h,w', [(1, 1), (1, 9), (9, 1), (4, 4), (7, 7)])
def test_degenerate_sizes_fall_back_to_redraw(h, w):
    r = _FigureRenderer(h, w, 'gray', 0.0, 1.0, overlay_kw=None)
    assert not r.fast
    frame = RNG.uniform(0.0, 1.0, (h, w))
    fast = r(frame, None, 0)
    assert np.array_equal(fast, r._render_reference(frame, None, 0))


def test_quantize_matches_matplotlib_index_lookup():
    """索引必须是 floor(g*256) 而不是 uint8(g*255): 后者会差 ±1。"""
    idx = _quantize(np.array([0.0, 1500.0, 1499.9, 3000.0, -5.0]),
                    0.0, 3000.0)
    assert idx.tolist() == [0, 128, 127, 255, 0]
