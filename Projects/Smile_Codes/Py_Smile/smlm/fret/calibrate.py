# smlm/fret/calibrate.py
"""γ / α / β 的全局标定, 以及可落盘的 ``FretCalibration``。

三个量的物理含义
----------------
``α`` 串色: 受体通道里的 ``I_A = S_A + α·I_D``。受体漂白、供体去淬灭后
``I_A^post = α·I_D^post``, 所以逐分子的 ``I_A^post / I_D^post`` 中位数就是 α。

``γ`` 通道增益比, 逐分子 ``ΔI_A / ΔI_D`` 的**中位数**。

* ``ΔI_A`` 直接用测到的红降 ``I_A^pre − I_A^post``, **串色不扣** —— 所以得到的是
  ``γ − α``; 要把串色扣掉就 ``add_alpha=True`` (等价于 ``+ α``)。
* 两个量都只在**漂白台阶附近**读数 (``local_step``, 受体侧 ±``win`` 帧), 不用远端
  平台, 也不用整条 trace 的平均: 远端会被别的闪烁/漂白污染, 全长平均会把背景漂移
  算成信号。供体一侧在 ``4·win`` 的长窗里取**亮帧**中位数 (抗闪烁), 否则一个
  恰好灭着的 4 帧窗口会把 ``ΔI_D`` 做小、把 γ 抬高。

三个估计量里推荐 ``γ_med``: 回归斜率(ols/rma)在两边都带噪声时偏差很大 ——
ΔI_D 被低估的点会把斜率压低。

``β`` 受体通道的暗电流/杂散底。用物理定义 ``β = median(I_A^post − α·I_D^post)``,
而不是回归截距。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
import os
import numpy as np

from .datatypes import passed, CRITERIA_STRICT

NPZ_VERSION = 1


# ---------------------------------------------------------------------------
#  稳健直线
# ---------------------------------------------------------------------------
def robust_line(x, y, k_clip=2.5, n_iter=6, min_n=6):
    """迭代 MAD 裁剪的最小二乘直线 ``y = slope*x + icpt``。返回 (slope, icpt, mask)。"""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 2:
        return np.nan, np.nan, np.zeros(x.size, bool)
    m = np.ones(x.size, bool)
    sl, ic = np.polyfit(x, y, 1)
    for _ in range(int(n_iter)):
        r = y - (sl * x + ic)
        s = 1.4826 * np.median(np.abs(r - np.median(r)))
        if not np.isfinite(s) or s <= 0:
            break
        m_new = np.abs(r - np.median(r)) < k_clip * s
        if m_new.sum() < min_n or np.array_equal(m_new, m):
            break
        m = m_new
        sl, ic = np.polyfit(x[m], y[m], 1)
    return float(sl), float(ic), m


# ---------------------------------------------------------------------------
#  模型
# ---------------------------------------------------------------------------
@dataclass
class FretCalibration:
    """仪器级的全局常数。与 ``DualViewCalibration`` 同构: 可 report / save / load。"""
    gamma: float = np.nan
    alpha: float = 0.0
    beta: float = 0.0
    n_alpha: int = 0
    n_gamma: int = 0
    gamma_rma: float = np.nan
    gamma_ols: float = np.nan
    detail: dict = field(default_factory=dict)
    version: int = NPZ_VERSION

    # ---------- 诊断 ----------
    def report(self):
        print(f'  FretCalibration: γ={self.gamma:.3f}  α={self.alpha:.3f}  '
              f'β={self.beta:.0f} photon')
        print(f'    n_gamma={self.n_gamma}  n_alpha={self.n_alpha}   '
              f'γ_rma={self.gamma_rma:.3f}  γ_ols={self.gamma_ols:.3f}')
        if not np.isfinite(self.gamma) or self.gamma <= 0:
            print('    [!] γ ≤ 0 或不可用: 这批"红降绿升"事件可能仍以噪声/错配为主。')

    # ---------- 落盘 ----------
    def save(self, path):
        d = os.path.dirname(os.path.abspath(path))
        if d:
            os.makedirs(d, exist_ok=True)
        np.savez_compressed(
            path,
            gamma=float(self.gamma), alpha=float(self.alpha), beta=float(self.beta),
            gamma_rma=float(self.gamma_rma), gamma_ols=float(self.gamma_ols),
            n_alpha=int(self.n_alpha), n_gamma=int(self.n_gamma),
            version=int(self.version),
            detail=np.array(json.dumps(self.detail, default=str)))
        return path

    @classmethod
    def load(cls, path):
        z = np.load(path, allow_pickle=False)
        kw = {k: z[k].item() for k in ('gamma', 'alpha', 'beta', 'gamma_rma',
                                       'gamma_ols', 'n_alpha', 'n_gamma')}
        kw['version'] = int(z['version'])
        if 'detail' in z:
            try:
                kw['detail'] = json.loads(str(z['detail']))
            except Exception:
                kw['detail'] = {}
        return cls(**kw)

    def __repr__(self):
        return (f'FretCalibration(γ={self.gamma:.3f}, α={self.alpha:.3f}, '
                f'β={self.beta:.0f}, n={self.n_gamma})')


# ---------------------------------------------------------------------------
#  台阶附近的局部读数
# ---------------------------------------------------------------------------
def local_step(ya, yd, j, lag, win=4, sd=None):
    """在 **bleach 台阶附近**直接量 ΔI_A / ΔI_D。

    只看紧贴台阶的那几帧, 不用远端平台, 也不用整条 trace 的平均 —— 远端平台
    会被前后别的闪烁/漂白污染, 整条 trace 的平均更是把背景漂移算进信号。

    供体一侧用**抗闪烁**口径 (:func:`~smlm.fret.events.post_level`): 在 ``4·win``
    帧的长窗里挑出 ``I_D^pre + 3σ_D`` 以上的亮帧取中位数。供体去淬灭后经常闪,
    固定取 4 帧的中位数会随机偏低, 直接把 γ 的分母做小。

    Parameters
    ----------
    ya, yd : (F,) ndarray
        受体(红) / 供体(绿) 强度 trace。
    j : int
        台阶帧(受体开始掉的**前一帧**的右边界)。
    lag : int
        供体去淬灭的滞后帧。
    win : int
        半窗长(帧)。受体取 ``[j-win, j)`` vs ``[j, j+win)``;
        供体取 ``[j-win, j)`` 作基线, 其后 ``4·win`` 帧作亮态窗。
    sd : float, optional
        供体通道噪声 σ; 不给就用 ``robust_noise(yd)`` 现估。

    Returns
    -------
    dict or None
        ``{IA_pre, IA_post, ID_pre, ID_post, dIA, dID, n_on}``;
        窗口越界或全 nan 返回 None。
    """
    from .events import post_level

    ya = np.asarray(ya, float)
    yd = np.asarray(yd, float)
    n = min(ya.size, yd.size)
    j = int(j)
    lag = int(max(0, lag))
    win = int(max(1, win))
    jd = j + lag
    span = int(4 * win)
    if j - win < 0 or j + win > n or jd + min(win, span) > n:
        return None

    def _m(y, a, b):
        s = y[a:b]
        s = s[np.isfinite(s)]
        return float(np.median(s)) if s.size else np.nan

    IA_pre = _m(ya, j - win, j)
    IA_post = _m(ya, j, j + win)
    ID_pre = _m(yd, j - win, j)

    hi = min(n, jd + span)
    s_d = max(float(robust_noise(yd)) if sd is None else float(sd), 1e-9)
    ID_post, n_on = post_level(yd, jd, hi, ID_pre + 3.0 * s_d,
                               max(2, win // 2))
    if not np.isfinite(ID_post):        # 没有可信亮态 -> 退回纯中位数
        ID_post, n_on = _m(yd, jd, hi), 0
    if not np.all(np.isfinite([IA_pre, IA_post, ID_pre, ID_post])):
        return None
    return dict(IA_pre=IA_pre, IA_post=IA_post, ID_pre=ID_pre, ID_post=ID_post,
                dIA=IA_pre - IA_post, dID=ID_post - ID_pre, n_on=int(n_on))


def gamma_points(events, traces_a=None, traces_d=None, *, win=4,
                 min_snr_d=2.0):
    """给 γ 用的 ``(ids, ΔI_A, ΔI_D)``。

    给了 ``traces_a``/``traces_d`` 就在**台阶附近**用 ``win`` 帧重测, 否则退回
    ``StepEvent`` 里已经量好的 ``dIA/dID``。``ΔI_A`` 一律是**原始**红降, 不扣串色
    (串色没扣就是没扣, 要不要在 γ 上加回 α 由调用方决定)。
    """
    ids, dAs, dDs = [], [], []
    use_tr = traces_a is not None and traces_d is not None
    for e in events:
        sd = float(e.sig_d) if np.isfinite(e.sig_d) else 0.0
        if use_tr:
            loc = local_step(traces_a[int(e.emitter_id)].intensity,
                             traces_d[int(e.emitter_id)].intensity,
                             e.frame, e.lag, win,
                             sd=sd if sd > 0 else None)
            if loc is None:
                continue
            da, dd = loc['dIA'], loc['dID']
        else:
            da, dd = float(e.dIA), float(e.dID)
        if not (np.isfinite(da) and np.isfinite(dd)) or dd <= 0:
            continue
        if dd <= float(min_snr_d) * sd:      # 分母是噪声的点不要
            continue
        ids.append(int(e.emitter_id))
        dAs.append(da)
        dDs.append(dd)
    return (np.asarray(ids, int), np.asarray(dAs, float),
            np.asarray(dDs, float))


# ---------------------------------------------------------------------------
#  估计
# ---------------------------------------------------------------------------
def estimate_alpha(events, min_snr_on_d=3.0):
    """α = median(I_A^post / I_D^post): 受体已灭、供体仍亮那一截的串色比。"""
    sel = [e for e in events
           if e.flags.get('post_d') and e.IA_post > 0 and e.ID_post > 0]
    if not sel:
        return float('nan'), 0, np.zeros(0)
    v = np.array([e.IA_post / e.ID_post for e in sel])
    return float(np.median(v)), len(sel), v


def estimate_beta(events, alpha):
    """β = median(I_A^post − α·I_D^post): 受体通道里非串色的残留。"""
    sel = [e for e in events if np.isfinite(e.IA_post)]
    if not sel:
        return 0.0, np.zeros(0)
    b = np.array([e.IA_post - alpha * e.ID_post for e in sel])
    return float(np.median(b)), b


def estimate_gamma(events, *, traces_a=None, traces_d=None, win=4,
                   min_snr_d=2.0, alpha=0.0, add_alpha=False):
    """返回 ``(gamma_med, gamma_rma, gamma_ols, n, mask, ratios)``。

    ``gamma_med`` 是逐分子 ``ΔI_A/ΔI_D`` 的**中位数**; ``ΔI_A`` 直接取测到的红降
    (串色没有扣)。``add_alpha=True`` 时把 α 加回去(等价于扣掉串色的那部分红降)。
    """
    ids, dA, dD = gamma_points(events, traces_a, traces_d, win=win,
                               min_snr_d=min_snr_d)
    if dD.size < 6:
        return (np.nan,) * 3 + (int(dD.size), np.zeros(dD.size, bool), np.zeros(0))

    off = float(alpha) if (add_alpha and np.isfinite(alpha)) else 0.0
    sl, _, mask = robust_line(dD, dA)
    cov = np.cov(dD, dA)[0, 1]
    rma = float(np.sign(cov) * dA.std() / max(dD.std(), 1e-12))
    ratios = dA / dD
    med = float(np.median(ratios))
    return (med + off, rma + off, float(sl) + off, int(dD.size), mask, ratios)


def calibrate(events, *, min_snr_d=2.0, traces_a=None, traces_d=None,
              win=4, add_alpha=False, verbose=True):
    """一次算齐 α / γ / β。``events`` 是 ``events.scan_steps`` 的输出。

    ``traces_a`` / ``traces_d`` 给了之后 γ 只用**台阶附近 ``win`` 帧**的信号重算;
    不给就用 ``StepEvent`` 里 ``sw`` 帧的读数。
    """
    good = [e for e in events if passed(e, CRITERIA_STRICT)]
    pool = good or list(events)

    alpha, n_alpha, av = estimate_alpha(pool)
    al = 0.0 if not np.isfinite(alpha) else alpha
    g_med, g_rma, g_ols, n_g, mask, ratios = estimate_gamma(
        pool, traces_a=traces_a, traces_d=traces_d, win=win,
        min_snr_d=min_snr_d, alpha=al, add_alpha=add_alpha)

    g_ids, dA, dD = gamma_points(pool, traces_a, traces_d, win=win,
                                 min_snr_d=min_snr_d)
    beta, bvals = estimate_beta(pool, al)

    cal = FretCalibration(gamma=g_med, alpha=alpha, beta=beta,
                          n_alpha=int(n_alpha), n_gamma=int(n_g),
                          gamma_rma=g_rma, gamma_ols=g_ols,
                          detail={'n_events': len(events),
                                  'n_good': len(good),
                                  'gamma_win': int(win),
                                  'gamma_add_alpha': bool(add_alpha),
                                  'alpha_p16_p84': [float(np.percentile(av, 16)),
                                                    float(np.percentile(av, 84))]
                                  if av.size else [],
                                  'ratio_p16_p84': [float(np.percentile(ratios, 16)),
                                                    float(np.percentile(ratios, 84))]
                                  if ratios.size else []})

    if verbose:
        print(f'\nα = I_A/I_D (串色) = {alpha:.3f}   n = {n_alpha}' if np.isfinite(alpha)
              else '\n[!] 没有可用于 α 的分子。')
        if np.isfinite(g_med):
            print(f'γ 估计 (台阶附近 ±{win} 帧, n = {n_g}):')
            print(f'  γ_med = {g_med:.3f}    <== 推荐   '
                  f'(逐分子 ΔI_A/ΔI_D 的中位数)')
            print(f'  γ_rma = {g_rma:.3f}     γ_ols = {g_ols:.3f}')
            if np.isfinite(ratios).any():
                p16, p50, p84 = np.nanpercentile(ratios, [16, 50, 84])
                print(f'  比值分布 p16~p84 = {p16:.2f}~{p84:.2f}   '
                      f'(中位 {p50:.2f}, 共 {ratios.size} 个)')
            if not add_alpha and np.isfinite(alpha):
                print(f'  [注] ΔI_A 未扣串色 -> 若要把 α 扣掉: γ = {g_med + al:.3f}')
            if dD.size >= 6:
                k = np.argsort(-np.abs(ratios))
                print('  最极端的 5 个:  #id    ΔI_A    ΔI_D   ratio')
                for t in k[:5]:
                    print(f'                 {g_ids[t]:>5} {dA[t]:8.0f} '
                          f'{dD[t]:8.0f} {ratios[t]:7.2f}')
            print(f'  β (漂白后残留 A 通道) = {beta:.0f} photon')
        else:
            print(f'[!] 合格点只有 {n_g} 个(需 ≥6), 无法稳健估计 γ。')
    return cal
