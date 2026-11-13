# smlm/fret/widgets.py
"""交互复核界面(ipywidgets 层)。

只有这一层依赖 ipywidgets; 状态在 ``review.ReviewSession``, 渲染在 ``viz``。

相比旧实现改了三件事
--------------------
1. **不再用 ``Output`` 装 matplotlib 图** —— 图会显示两遍, 且 ``clear_output``
   在 VS Code 里不可靠, 切条目时旧图会叠上去。现在 trace 也渲染成 PNG,
   塞进第二个 ``W.Image``。
2. **视频帧保持原生分辨率**(41x84), 前端用 CSS 放大 —— 单帧从 ~100 KB 降到
   2~5 KB, 播放不再卡。
3. **trace 里画出了估计的背景**(点线), 方便判断环形背景是不是被邻居污染了。
"""
from __future__ import annotations

import ipywidgets as W

_PIXELATED_CSS = '''<style>
.fret-pix img { image-rendering: pixelated; image-rendering: crisp-edges; }
</style>'''


def make_ui(session, *, width_px=430, fps_ms=90, scale_hint=5):
    """给一个 ``ReviewSession`` 建界面(不 display, 返回 widget)。"""
    img = W.Image(format='png', layout=W.Layout(border='1px solid #ccc'))
    img.add_class('fret-pix')
    img.layout.width = f'{width_px}px'

    img_trace = W.Image(format='png')
    img_trace.layout.width = '760px'

    sld = W.IntSlider(min=0, max=1, value=0, description='frame',
                      continuous_update=True, layout=W.Layout(width='340px'))
    play = W.Play(min=0, max=1, value=0, interval=fps_ms, step=1, repeat=True)
    lbl = W.HTML()
    lbl_trace = W.HTML()

    b_keep = W.Button(description='✔ 保留', button_style='success')
    b_rej = W.Button(description='✘ 不要', button_style='danger')
    b_back = W.Button(description='← 上一条')
    b_end = W.Button(description='结束并汇总')

    st = {'pngs': [], 'last': -1}

    # ---------- 帧切换 ----------
    def _render(k):
        p = st['pngs']
        if not p or not (0 <= k < len(p)) or k == st['last']:
            return
        st['last'] = k
        img.value = p[k]

    def _on_sld(ch):
        _render(int(ch['new']))

    def _on_play(ch):
        _render(int(ch['new']))

    sld.observe(_on_sld, names='value')
    play.observe(_on_play, names='value')

    # ---------- 换一条 ----------
    def _show():
        ev = session.current()
        done = ev is None
        for b in (b_keep, b_rej, b_back):
            b.disabled = done
        if done:
            st['pngs'], st['last'] = [], -1
            img.value = b''
            img_trace.value = b''
            lbl.value = f'<b>已审完</b> &nbsp; {session.summary()}'
            lbl_trace.value = ''
            return

        i = int(ev.emitter_id)
        pngs, frames = session.movie(ev)
        st['pngs'], st['last'] = pngs, -1
        sld.max = play.max = max(0, len(pngs) - 1)
        sld.value = play.value = 0
        _render(0)

        r_d = session.trace_d[i].meta.get('radius', 5)
        r_a = session.trace_a[i].meta.get('radius', 4)
        lbl.value = (
            f"<b>{session.index + 1} / {len(session)}</b> &nbsp; "
            f"emitter <b>#{i}</b> &nbsp; "
            f"bleach @ {ev.frame * session.ex_time:.1f} s &nbsp; "
            f"ΔI_A={ev.dIA:.0f} &nbsp; ΔI_D={ev.dID:.0f} &nbsp; lag={ev.lag} &nbsp; "
            f"<span style='color:#060'>(左=供体 r={r_d})</span> "
            f"<span style='color:#a00'>(右=受体 r={r_a})</span> &nbsp; "
            f"每帧 {session.ex_time * max(1, (frames[1] - frames[0])):.2f} s"
            if len(frames) > 1 else
            f"<b>{session.index + 1} / {len(session)}</b> &nbsp; emitter #{i}")

        img_trace.value = session.trace_png(ev)
        lbl_trace.value = ('<span style="font-size:11px;color:#555">'
                           '实线 = 孔径积分(已扣背景)　点线 = 扣掉的背景　'
                           '虚线 = MLE 光子数(不受环形背景影响)</span>')

    # ---------- 按钮 ----------
    def _decide(keep):
        session.decide(keep)
        _show()

    def _back(b):
        session.back()
        _show()

    def _end(b):
        print(session.summary())
        print(f'保留 {session.keep_ids}')

    b_keep.on_click(lambda b: _decide(True))
    b_rej.on_click(lambda b: _decide(False))
    b_back.on_click(_back)
    b_end.on_click(_end)

    ui = W.VBox([W.HTML(_PIXELATED_CSS), lbl,
                 W.HBox([img, W.VBox([sld, play])]),
                 img_trace, lbl_trace,
                 W.HBox([b_keep, b_rej, b_back, b_end])])
    _show()
    return ui
