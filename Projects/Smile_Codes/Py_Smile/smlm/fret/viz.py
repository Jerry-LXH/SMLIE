# smlm/fret/viz.py
"""复核界面用的渲染: 直接吐 PNG bytes。

为什么不把 matplotlib 图丢进 ``ipywidgets.Output``
-----------------------------------------------
inline 后端里 figure 的归属取决于创建时刻的 capture 上下文; 在 ``Output`` 之外
创建、又在里面 ``plt.show()``, 很容易被 display 两次。再加上 ``clear_output`` 在
VS Code 里不可靠, 切换条目时旧图会留在原地叠上去。

所以这里一律"渲染成 PNG → 塞进 ``W.Image``": 确定、无重复、快。

视频帧也**不做放大**: 41x84 的灰度噪声图放大 5 倍后 PNG 压不动(单帧 60~120 KB),
200 帧按 80 ms 播放就是 1~2 MB/s 的通信量 —— 这才是卡的原因。保持原生分辨率
(2~5 KB/帧), 放大交给前端 CSS。
"""
from __future__ import annotations

import io
import numpy as np


# ---------------------------------------------------------------------------
#  图像块
# ---------------------------------------------------------------------------
def crop_centered(vol, cy, cx, half):
    """(T,H,W) -> (T,2h+1,2h+1), 越界补 0。"""
    vol = np.asarray(vol)
    F, H, W = vol.shape
    n = 2 * int(half) + 1
    r0, c0 = int(round(cy)) - int(half), int(round(cx)) - int(half)
    out = np.zeros((F, n, n), vol.dtype if vol.dtype.kind == 'f' else np.float32)
    rr0, rr1 = max(0, r0), min(H, r0 + n)
    cc0, cc1 = max(0, c0), min(W, c0 + n)
    if rr1 > rr0 and cc1 > cc0:
        out[:, rr0 - r0:rr1 - r0, cc0 - c0:cc1 - c0] = vol[:, rr0:rr1, cc0:cc1]
    return out


def to_uint8(a, pct=99.7):
    hi = float(np.percentile(a, pct)) if a.size else 1.0
    return np.clip(a / (hi if hi > 0 else 1.0) * 255.0, 0, 255).astype(np.uint8)


def _draw_ring(a, cy, cx, rad, val=255):
    """在 uint8 图上画半径 rad 的方框(越界自动裁)。"""
    H, W = a.shape
    r0, r1 = max(0, cy - rad), min(H - 1, cy + rad)
    c0, c1 = max(0, cx - rad), min(W - 1, cx + rad)
    a[r0, c0:c1 + 1] = val
    a[r1, c0:c1 + 1] = val
    a[r0:r1 + 1, c0] = val
    a[r0:r1 + 1, c1] = val


def _png(arr_gray, scale=1):
    from PIL import Image
    im = Image.fromarray(arr_gray, 'L')
    if scale != 1:
        im = im.resize((im.width * scale, im.height * scale), Image.NEAREST)
    buf = io.BytesIO()
    im.save(buf, format='PNG', optimize=False)
    return buf.getvalue()


def render_movie_pngs(img_d, img_a, centre_d, centre_a, *, half=20, step=10,
                      radius_d=5, radius_a=4, gap=2, pct=99.7):
    """双通道小视频帧(左=供体/绿, 右=受体/红), **原生分辨率**。

    Returns
    -------
    (pngs, frames)
    """
    g = crop_centered(img_d, centre_d[0], centre_d[1], half)
    r = crop_centered(img_a, centre_a[0], centre_a[1], half)
    g8, r8 = to_uint8(g, pct), to_uint8(r, pct)
    n = 2 * int(half) + 1
    sep = np.full((n, int(gap)), 255, np.uint8)

    frames = np.arange(0, g.shape[0], int(step))
    pngs = []
    for f in frames:
        a = np.concatenate([g8[f], sep, r8[f]], axis=1).copy()
        _draw_ring(a, half, half, radius_d)                    # 供体孔径
        _draw_ring(a, half, n + int(gap) + half, radius_a)     # 受体孔径
        pngs.append(_png(a))
    return pngs, frames


# ---------------------------------------------------------------------------
#  trace
# ---------------------------------------------------------------------------
def render_trace_png(trace_a, trace_d, ev=None, *, ex_time=1.0, mle_a=None,
                     mle_d=None, pre_win=100, post_win=400, cal=None,
                     figsize=(7.6, 3.0), dpi=105, title=None):
    """双通道 trace + **估计的背景** + MLE 对照。

    实线 = 孔径积分(已扣背景); 点线 = 从孔径里扣掉的背景(photon/frame);
    虚线 = 逐帧 MLE 光子数(不受环形背景影响, 用来判断背景是否可信)。
    """
    import matplotlib.pyplot as plt

    ia = np.asarray(trace_a.intensity, float)
    idd = np.asarray(trace_d.intensity, float)
    n = min(ia.size, idd.size)
    ia, idd = ia[:n], idd[:n]
    t = np.arange(n) * float(ex_time)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    if ev is not None:
        j, lag = int(ev.frame), int(ev.lag)
        ax.axvspan(t[max(0, j - pre_win)], t[min(n - 1, j)], color='0.90', lw=0)
        ax.axvspan(t[min(n - 1, j + lag)], t[min(n - 1, j + lag + post_win)],
                   color='0.80', lw=0)
        ax.axvline(t[min(j, n - 1)], color='k', ls='--', lw=0.8)

    # 孔径 + 背景 + MLE
    ax.plot(t, idd, color='green', lw=0.8, label='donor aperture')
    bg_d = trace_d.meta.get('bg_flux')
    if bg_d is not None:
        ax.plot(t, np.asarray(bg_d)[:n], color='green', lw=0.9, ls=':',
                alpha=0.9, label='donor bg (subtracted)')
    if mle_d is not None:
        ax.plot(t, np.asarray(mle_d)[:n], color='green', lw=0.5, ls='--',
                alpha=0.75, label='donor MLE')

    ax.plot(t, ia, color='red', lw=0.8, alpha=0.9, label='acceptor aperture')
    bg_a = trace_a.meta.get('bg_flux')
    if bg_a is not None:
        ax.plot(t, np.asarray(bg_a)[:n], color='red', lw=0.9, ls=':',
                alpha=0.9, label='acceptor bg (subtracted)')
    if mle_a is not None:
        ax.plot(t, np.asarray(mle_a)[:n], color='red', lw=0.5, ls='--',
                alpha=0.75, label='acceptor MLE')

    ax.set_xlabel('time (s)', fontsize=7)
    ax.set_ylabel('photons / frame', fontsize=7)
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=5.5, loc='upper right', frameon=False, ncol=2)
    if title is None and ev is not None:
        title = (f'#{ev.emitter_id}   A {ev.IA_pre:.0f}→{ev.IA_post:.0f}   '
                 f'D {ev.ID_pre:.0f}→{ev.ID_post:.0f}   lag={ev.lag}')
    if title:
        ax.set_title(title, fontsize=7)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    plt.close(fig)
    return buf.getvalue()


# ---------------------------------------------------------------------------
#  汇总图(静态, 给 notebook 用)
# ---------------------------------------------------------------------------
def plot_gamma(events, cal, ax=None, raw=True):
    """ΔI_A vs ΔI_D 散点 + γ 三条估计线。"""
    import matplotlib.pyplot as plt
    dD = np.array([e.dID for e in events], float)
    dA = np.array([e.dIA for e in events], float)
    ok = np.isfinite(dD) & np.isfinite(dA)
    dD, dA = dD[ok], dA[ok]
    ax = ax or plt.gca()
    ax.scatter(dD, dA, s=16, alpha=0.75, c='purple')
    if dD.size:
        xs = np.linspace(0, float(np.nanmax(dD)) * 1.05, 50)
        for g, c, ls, lab in ((cal.gamma, 'r', '-', 'γ_med'),
                              (cal.gamma_rma, 'b', '-.', 'γ_rma'),
                              (cal.gamma_ols, 'g', '--', 'γ_ols')):
            if np.isfinite(g):
                ax.plot(xs, g * xs, color=c, ls=ls, lw=1.1,
                        label=f'{lab}={g:.2f}')
        lim = [0, max(float(np.nanmax(dD)), float(np.nanmax(dA))) * 1.05]
        ax.plot(lim, lim, 'k--', lw=0.8, label='γ=1')
    ax.set_xlabel('ΔI_D (photons)')
    ax.set_ylabel('ΔI_A (photons)')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    return ax
