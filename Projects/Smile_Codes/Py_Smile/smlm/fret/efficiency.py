# smlm/fret/efficiency.py
"""FRET 效率。

    S_A = I_A^pre − α·I_D^pre − β
    E   = S_A / (S_A + γ·I_D^pre)

窗口取法
--------
``pre``  : 台阶左侧最多 ``pre_win`` 帧; 不足就从第 0 帧开始, 并计入 ``n_short``。
           **只平均受体亮着的帧** (``gate_pre_d``): 受体熄灭的间隙里供体去淬灭,
           那几帧的亮度是"没转移"的状态, 混进来会把 I_D^pre 抬高、把 E 压低。
           受体平稳亮着时掩码全 True, 与不做门控等价。
``post`` : 台阶右侧最多 ``post_win`` 帧, 且只平均**绿光 > post_thr_d** 的帧 ——
           受体漂白后供体去淬灭会再闪烁, 暗帧会把 I_D^post 压低。

默认**不做任何截断**: I_D^pre 为负就让它为负, E 算出 >1 或 <0 都照原样留着。
截断/下界属于下游分析的选择, 不该埋在计算里。
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .datatypes import StepEvent, FretResult
from .calibrate import FretCalibration


@dataclass(frozen=True)
class EfficiencySpec:
    pre_win: int = 100          # 台阶左侧最多取多少帧
    post_win: int = 400         # 台阶右侧最多取多少帧
    post_thr_d: float = 100.0   # post 窗口里只平均"绿光 > 该值"的帧
    min_on: int = 10            # post 窗口内至少要有这么多帧在阈值以上
    min_pre: int = 0            # pre 窗最少要有多少帧 (0 = 只要求 ≥2)
    gate_pre_d: bool = True     # pre 窗只平均受体亮着的帧
    clip_n: float = 3.0
    clip_iter: int = 2

    def report(self):
        pre = f'pre≤{self.pre_win} 帧'
        if self.min_pre:
            pre += f' (至少 {self.min_pre})'
        if self.gate_pre_d:
            pre += ' 只取受体亮帧'
        return (f'EfficiencySpec({pre}, post≤{self.post_win} 帧, '
                f'post 只取绿>{self.post_thr_d:.0f}, 至少 {self.min_on} 帧)')


def window_mean(x, n_clip=3.0, n_iter=2):
    """sigma-clip 均值: 既是"平均亮度", 又不会被单个尖峰带跑。"""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    for _ in range(int(n_iter)):
        if x.size < 3:
            break
        m, s = float(x.mean()), float(x.std())
        if s <= 0:
            break
        y = x[np.abs(x - m) <= float(n_clip) * s]
        if y.size < max(2, x.size // 2):
            break
        x = y
    return float(x.mean()) if x.size else np.nan


def fret_efficiency(events, traces_a, traces_d, cal: FretCalibration,
                    spec: EfficiencySpec = None, *, t_bleach=None, verbose=True):
    """对一批事件算 E。``traces_a`` = 受体(红), ``traces_d`` = 供体(绿)。"""
    spec = spec or EfficiencySpec()
    gamma = float(cal.gamma) if np.isfinite(cal.gamma) else float('nan')
    alpha = float(cal.alpha) if np.isfinite(cal.alpha) else 0.0
    beta = float(cal.beta) if np.isfinite(cal.beta) else 0.0

    out, stats = [], {'n_short': 0, 'n_noon': 0, 'n_drop': 0, 'n_pre': 0,
                      'n_gate': 0, 'n_in': len(events)}
    n_tr = min(len(traces_a), len(traces_d))
    bad = [int(ev.emitter_id) for ev in events
           if not 0 <= int(ev.emitter_id) < n_tr]
    if bad:
        # 负下标会静默读到别的 trace, 必须响亮报错
        raise ValueError(
            f'{len(bad)} 个事件的 emitter_id 越界 (合法范围 [0, {n_tr - 1}]): '
            f'{bad[:8]}{" ..." if len(bad) > 8 else ""}\n'
            f'  常见原因: 事件不是从 scan_steps / select_pairs 出来的, '
            f'emitter_id 还留着 find_step 的 -1。')
    for ev in events:
        i, j = int(ev.emitter_id), int(ev.frame)
        ya = np.asarray(traces_a[i].intensity, float)
        yd = np.asarray(traces_d[i].intensity, float)
        n = min(ya.size, yd.size)
        l0, l1 = max(0, j - spec.pre_win), j
        p0, p1 = j, min(n, j + spec.post_win)
        if l1 - l0 < 2 or p1 - p0 < 2:
            continue
        if l1 - l0 < spec.min_pre:
            stats['n_pre'] += 1
            continue
        if l1 - l0 < spec.pre_win:
            stats['n_short'] += 1

        # pre 窗: 只统计受体亮着的帧 (受体熄灭时供体去淬灭, 那几帧不属于"转移态")
        gate = getattr(ev, 'red_on', None) if spec.gate_pre_d else None
        if gate is not None:
            gate = np.asarray(gate, bool)[l0:l1]
        sel = np.isfinite(ya[l0:l1]) if gate is None else gate
        n_gate = int(sel.sum()) if gate is not None else int(l1 - l0)
        if n_gate < 2:
            stats['n_gate'] += 1
            continue

        IA_pre = window_mean(ya[l0:l1][sel], spec.clip_n, spec.clip_iter)
        ID_pre = window_mean(yd[l0:l1][sel], spec.clip_n, spec.clip_iter)

        seg_d, seg_a = yd[p0:p1], ya[p0:p1]
        on = np.isfinite(seg_d) & (seg_d > spec.post_thr_d)
        if on.sum() < spec.min_on:
            stats['n_noon'] += 1
            continue
        ID_post = window_mean(seg_d[on], spec.clip_n, spec.clip_iter)
        IA_post = window_mean(seg_a[on], spec.clip_n, spec.clip_iter)

        SA = IA_pre - alpha * ID_pre - beta
        den = SA + gamma * ID_pre
        if (not np.isfinite(SA)) or (not np.isfinite(den)) or abs(den) < 1e-9 \
                or not np.isfinite(ID_post) or ID_post <= 0:
            stats['n_drop'] += 1
            continue

        sig_a = ev.sig_a if np.isfinite(ev.sig_a) else np.nan
        out.append(FretResult(
            emitter_id=i, frame=j, lag=int(ev.lag),
            n_pre=int(l1 - l0), n_on=int(on.sum()),
            IA_pre=IA_pre, ID_pre=ID_pre, IA_post=IA_post, ID_post=ID_post,
            S_A=float(SA), E=float(SA / den),
            E_don=float(1.0 - ID_pre / ID_post),
            snr_A=float(SA / max(sig_a, 1e-9)) if np.isfinite(sig_a) else np.nan,
            t_bleach=float(t_bleach(j)) if callable(t_bleach) else np.nan))

    if verbose:
        _report(out, stats, spec, gamma, alpha, beta)
    return out, stats


def _report(res, stats, spec, gamma, alpha, beta):
    for k, msg in (('n_pre', f'pre 窗不足 min_pre={spec.min_pre} 帧'),
                   ('n_short', f'pre 不足 {spec.pre_win} 帧'),
                   ('n_gate', 'pre 窗里受体亮着的帧不足 2 (无法定 I_D^pre)'),
                   ('n_noon', f'post 窗内绿 > {spec.post_thr_d:.0f} 的帧不足 '
                              f'{spec.min_on}'),
                   ('n_drop', '分母 ~0 / I_D^post ≤ 0 / 非有限')):
        if stats[k]:
            print(f'  {msg}: {stats[k]} 条')
    if not res:
        print(f'[!] {stats.get("n_in", 0)} 个事件一条都没算出 E。'
              f'上面列的就是卡住的地方; 若一条原因都没列, '
              f'说明窗口长度不够 (pre/post < 2 帧)。')
        return
    E = np.array([r.E for r in res])
    ED = np.array([r.E_don for r in res])
    IDp = np.array([r.ID_pre for r in res])
    print(f'窗口: pre ≤ {spec.pre_win} 帧 (中位 {np.median([r.n_pre for r in res]):.0f})'
          f'   post ≤ {spec.post_win} 帧, 只取绿 > {spec.post_thr_d:.0f}'
          f' (中位 {np.median([r.n_on for r in res]):.0f})')
    print(f'\nE 分布 (n = {len(res)}, 未截断):')
    print(f'  中位 {np.median(E):.3f}   p16~p84 = '
          f'{np.percentile(E, 16):.3f}~{np.percentile(E, 84):.3f}')
    print(f'  范围 {E.min():.3f} ~ {E.max():.3f}')
    print(f'  E>1: {int((E > 1).sum())}   E<0: {int((E < 0).sum())}   '
          f'I_D^pre<0: {int((IDp < 0).sum())}')
    print(f'  (γ={gamma:.2f}, α={alpha:.3f}, β={beta:.0f})')
