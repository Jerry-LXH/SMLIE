"""two_emitter.level2 —— 亮度也未知: 把每个分子的亮度 Φ 当待估参数。

与 level1 的关系
----------------
**模型完全相同**（4 态联合链 + 帧内占空边际化 + 泊松发射），唯一区别是：

===========  ==========================  ==========================
             level1                       level2
===========  ==========================  ==========================
待估参数     ``r = (x1,y1,x2,y2)``        ``θ = (x1,y1,x2,y2, Φ1, Φ2[, β])``
Φ            已知（取自 metadata）          **待估**（``tie_phi=True`` 时两个分子共享一个 Φ）
β            固定（ROI 外圈估计）          ``fit_bg=True`` 时一并估计
E 步         前向-后向                      **完全一样**（状态/占空后验不依赖 Φ）
M 步         4 维 Nelder-Mead              6 维（或 5/7 维）Nelder-Mead
===========  ==========================  ==========================

为什么 E 步不用改
-----------------
完全数据对数似然里

    log p(D,Z,n | θ) = log π(Z_0) + Σ_t [log q(n) ]  + Σ_t log p(D_t | n, θ)
                       └── 与 θ 无关（转移矩阵已知）──┘   └── 只通过 μ 依赖 θ ──┘

而 μ_{t,p} = β + (n1/M)Φ1 P1(r1)_p + (n2/M)Φ2 P2(r2)_p。
状态后验与占空后验 π̄ 只依赖 q 和 D，**不含 Φ**；所以 E 步与 level1 逐字相同，
只是 M 步要多估 2 个参数。EM 的单调性因此原样保留。

已知的耦合 / 陷阱
-----------------
1. **标签交换**：两个分子不可区分 ⟺ (Φ1, r1) ↔ (Φ2, r2)。误差统计必须取两种配对里
   较小的那个（``fit_pair_level2`` 已处理）。
2. **Φ 与 ROI 截断**：P_i 的和 = ROI 圈入比例（<1），而 Φ 定义为"整帧全亮"的光子数。
   两者是自洽的 —— 圈到 ROI 外的光子本来就不会被探测，不影响 Φ 的可辨识性（只要
   ROI 足够大）。
3. **Φ 与 β 的弱耦合**：若 ROI 内几乎没有"两个都暗"的帧，β 与 Φ 会部分退化。
   ``fit_bg=True`` 时建议同时提供已知的 β 初值。
   **实测（d=40 数据集，template PSF，4 对平均）**：

   | β | Φ_fit (真值 1200) | Φ 相对偏差 | d_est (真值 40) | d 偏差 |
   |---|---|---|---|---|
   | 2.000（能量守恒值） | 1238.1 | **+3.18%** | 40.19 | +0.19 nm |
   | 2.200 | 1228.1 | +2.34% | 39.95 | **−0.05 nm** |
   | 2.462（边缘估计） | 1216.1 | +1.34% | 39.66 | −0.34 nm |
   | 2.700 | 1206.1 | +0.51% | 39.39 | −0.61 nm |

   → Φ 与 β 沿一条**ridge** 反向滑动（β 每 +0.1，Φ 约 −0.4%），
   **d 对 β 的近优值 (≈2.2) 略高于能量守恒值 (2.0)**。
   结论：**一般应固定 β**（用能量守恒值），不要指望 ``fit_bg=True`` 同时救 Φ。

4. **★ Φ 的可辨识上限 = 实际占空**：模型的占空被 ``k_on/(k_on+k_off)`` **硬编码**
   （本例 0.5），但有限长闪烁序列的**实际**占空会漂移。本数据集实测
   $\langle S\rangle = 1516.58$ ⟹ 等效占空 **0.522**（+4.4%；$T$=1500、相关时间
   $\sim1/(k_{on}+k_{off})$=0.1 s ⟹ 期望涨落 ±0.029，即 0.8σ，完全正常）。
   模型无法吸收这个偏差，于是 Φ 把它整个吃掉 —— 这就是 β=2.0 时 Φ 偏高 +3.2% 的
   **全部**原因（+3.2% ≈ +4.4% 的占空过剩，其余被位置自由度吸收）。
   ⟹ **想得到无偏 Φ 必须同时估计速率（占空）**，这是 Level 2b 的内容；
   在占空固定的前提下，Φ 的误差下界就是 $O(\pm 3\%)$。

5. **Φ 与间距的耦合（本文最关心）**：PSF 形状/宽度失配会同时污染 Φ 和 d ——
   所以 level2 默认 ``psf_model='template'``（用仿真器自己的 PSF 核）。

★ **重要：亮度比值 ↔ 间距 存在强退化，必须用多起点**

实测（合成数据，$d=40$ nm $\ll \sigma_{\rm PSF}=99$ nm，$T$=1500 帧，$M$=4）：

- **Φ 的"和"被牢牢约束**（$\Phi_1+\Phi_2$ 守恒），但**"怎么分"几乎不被约束**：
  图像主要只锁定 $\Phi_1\Phi_2/(\Phi_1+\Phi_2)^2 \cdot d^2$ 这类组合，
  于是 $(\Phi_1,\Phi_2,d)$ 沿一条**极平坦的山谷**滑动。
- 后果是**强烈的局部最优**：从对称初值出发能精确落在真值
  （$\Phi\to(1201,799)$、$d\to39.98$ nm，$\Delta\log L$ 仅 +1.3 ≈ 1σ），
  但从不对称初值出发会掉进伪解（$\Phi\to(832,1567)$、$d\to46.6$ nm，
  $\Delta\log L=-5975$，即 2ΔlogL≈12000，统计上被**强烈排除**）。
- 因此 ``fit`` 默认 **``multi_start=3``**：分别从
  $(\hat\Phi,\hat\Phi)$、$(0.7\hat\Phi,1.3\hat\Phi)$、$(1.3\hat\Phi,0.7\hat\Phi)$
  出发（$\hat\Phi=$ ``estimate_phi_sum()``/2，由总光子数粗估），**取末轮 logL 最高者**。
  返回字典里的 ``starts`` 给出每个起点的终点 ``(Φ1, Φ2, d, logL)``——
  **主要看 $d$ 是否一致**（一致即稳健）；``logL_margin`` 是"最优比次优高多少"，
  接近 0 且 $d$ 不一致才说明山谷内多点等价、结果不可信
  （$d$ 一致时 margin=0 是良性的：多起点收敛到同一点）。
- ``tie_phi=True``（两分子共享一个 Φ）没有这个退化，**推荐在"亮度可视为相同"时优先使用**。
"""
from __future__ import annotations

import numpy as np

from .level1 import JointHMMTwoEmitter, build_pair_hmm

# 归一化尺度: 让 Nelder-Mead 的 6 个维度量级相当
_SCALE_R = 20.0        # nm
_SCALE_PHI = 200.0     # 光子/帧
_SCALE_BG = 0.5        # 光子/px/帧


class JointHMMTwoEmitter2(JointHMMTwoEmitter):
    """亮度未知的 4 态联合 HMM（继承 level1，仅扩参 + 换 M 步）。

    参数向量
    --------
    ``θ = [x1, y1, x2, y2] (+ [Φ1] 或 [Φ1, Φ2]) (+ [β])``

    构造参数（在 level1 的基础上新增）
    ----------------------------------
    phi2     : 分子 2 的 Φ 初值；默认 = ``phi_photons``
    tie_phi  : True 时两个分子共享一个 Φ（θ 少一维）
    fit_bg   : True 时把 β 也放进 θ 一起估
    """

    def __init__(self, *args, phi2=None, tie_phi=False, fit_bg=False, **kw):
        super().__init__(*args, **kw)
        self.phi2 = float(self.phi if phi2 is None else phi2)
        self.tie_phi = bool(tie_phi)
        self.fit_bg = bool(fit_bg)
        self.theta_ = None

    # ------------------------------------------------------------ θ 打包
    @property
    def n_phi(self) -> int:
        return 1 if self.tie_phi else 2

    @property
    def theta_size(self) -> int:
        return 4 + self.n_phi + (1 if self.fit_bg else 0)

    def scales(self) -> np.ndarray:
        """θ 各分量的物理尺度（用于归一化到 O(1) 再交给 Nelder-Mead）。"""
        s = [_SCALE_R] * 4
        s += [_SCALE_PHI] * self.n_phi
        if self.fit_bg:
            s += [_SCALE_BG]
        return np.asarray(s, float)

    def pack_theta(self, r, phi1=None, phi2=None, bg=None) -> np.ndarray:
        r = np.asarray(r, float).ravel()[:4]
        p1 = self.phi if phi1 is None else float(phi1)
        p2 = (p1 if self.tie_phi else
              (self.phi2 if phi2 is None else float(phi2)))
        th = [r[0], r[1], r[2], r[3], p1]
        if not self.tie_phi:
            th.append(p2)
        if self.fit_bg:
            th.append(self.bg if bg is None else float(bg))
        return np.asarray(th, float)

    def unpack_theta(self, theta):
        """θ → ``(r(4,), Φ1, Φ2, β)``。"""
        th = np.asarray(theta, float).ravel()
        r = th[:4]
        p1 = float(th[4])
        if self.tie_phi:
            p2, i = p1, 5
        else:
            p2, i = float(th[5]), 6
        bg = float(th[i]) if self.fit_bg else self.bg
        return r, p1, p2, bg

    # ------------------------------------------------------------ 同步状态
    def set_theta(self, theta) -> np.ndarray:
        """把 θ 写回实例状态（``phi/phi2/bg``），返回位置 ``r``。

        E 步（``loglik_table``）读的是实例属性，所以每次前向-后向之前都要同步。
        """
        r, p1, p2, bg = self.unpack_theta(theta)
        self.phi, self.phi2, self.bg = p1, p2, bg
        return r

    def log_likelihood_theta(self, theta) -> float:
        """``log p(D_{1:T} | θ)``（前向算法，先把 θ 同步进实例）。"""
        return self.log_likelihood(self.set_theta(theta))

    # ------------------------------------------------------------ 发射模型
    def _mu_matrix(self, theta) -> np.ndarray:
        """(K, P) 期望像素强度；**参数是 θ（可含 Φ,β）而不是单纯的 r**。"""
        r, p1, p2, bg = self.unpack_theta(theta)
        P1 = self.psf_pixels(r[0:2])
        P2 = self.psf_pixels(r[2:4])
        M = self.M
        n1s, n2s = self._sub_grid
        return (bg
                + np.outer(n1s * (p1 / M), P1.ravel())
                + np.outer(n2s * (p2 / M), P2.ravel()))

    # ------------------------------------------------------------ M 步
    def mstep_theta(self, theta, stats, maxiter: int = 150, xtol: float = 1e-3):
        """``argmax_θ Q(θ; π̄)`` —— 6(或5/7)维 Nelder-Mead（带尺度归一 + 精修）。

        ``Φ, β`` 必须为正，越界返回一个巨大的惩罚值（保持目标函数"处处有定义"，
        避免 NM 单纯形撞到非物理区域后 NaN）。
        """
        from scipy.optimize import minimize
        s = self.scales()
        y0 = np.asarray(theta, float) / s

        def f(y):
            th = y * s
            r, p1, p2, bg = self.unpack_theta(th)
            if not np.isfinite(th).all() or p1 <= 1.0 or p2 <= 1.0 or bg <= 1e-6:
                return 1e30
            return -self.q_value_fast(th, stats)

        def run(y_start, it, xtol_):
            step = np.maximum(np.abs(y_start) * 0.05, 0.02)
            simp = np.vstack([y_start, y_start + np.diag(step)])
            return minimize(f, y_start, method="Nelder-Mead",
                            options=dict(maxiter=int(it), xatol=xtol_,
                                         fatol=1e-3, initial_simplex=simp))

        res = run(y0, maxiter, xtol)
        y = np.asarray(res.x, float)
        res2 = run(y, max(30, int(maxiter) // 3), xtol * 0.2)      # 精修
        if res2.fun <= res.fun:
            res = res2
        return np.asarray(res.x, float) * s

    def mstep_q(self, *a, **kw):        # 防止误用 level1 的 4 维 M 步
        raise TypeError("level2 用 mstep_theta(theta, stats)；"
                        "mstep_q 只适用于 4 维的 level1")

    # ------------------------------------------------------------ 初值 & EM
    def initial_theta(self, r_init=None, phi1=None, phi2=None, bg=None):
        r = np.asarray(r_init if r_init is not None else self.initial_positions(),
                       dtype=float)
        return self.pack_theta(r, phi1, phi2, bg)

    def estimate_phi_sum(self, r_init=None) -> float:
        """从总光子数粗估 ``Φ1 + Φ2``（两级粗估，用于多起点的中心）。

        ``Σ_{t,p}(D − β) = T·Σ_i Φ_i·E[占空]·(P_i 的 ROI 圈入比例)``，
        取平稳占空 ``p_on``、并用当前位置估出的 PSF 圈入比例代替。
        """
        r = self.centroid_nm() if r_init is None else np.asarray(r_init, float).ravel()[:4]
        if r.size < 4:
            r = np.concatenate([r, r[:2]])
        cap = 0.5 * (self.psf_pixels(r[0:2]).sum() + self.psf_pixels(r[2:4]).sum())
        p_on = float(self.prior1[0])
        tot = float(self.data.sum() - self.bg * self.T * self.h * self.w)
        if tot <= 0 or cap <= 0 or p_on <= 0:
            return max(2.0 * self.phi, 1.0)
        return max(tot / (self.T * p_on * cap), 1.0)

    def initial_thetas(self, r_init=None, n_start: int = 3, phi_init=None,
                       spreads=((1.0, 1.0), (0.7, 1.3), (1.3, 0.7),
                                (0.45, 1.55), (1.55, 0.45))):
        """生成多起点 θ0 列表（围绕 ``Φ̂`` 展开，覆盖亮度比值方向）。

        见模块 docstring：$(\Phi_1,\Phi_2,d)$ 有平坦山谷，必须多起点 + 取最高 logL。
        ``tie_phi=True`` 时只有一个 Φ，退化为单起点。
        """
        r0 = np.asarray(r_init if r_init is not None else self.initial_positions(),
                        dtype=float)
        if self.tie_phi:
            p0 = float(phi_init if phi_init is not None else
                       self.estimate_phi_sum(r0) / 2.0)
            return [self.pack_theta(r0, p0)]
        base = float(phi_init if phi_init is not None else
                     self.estimate_phi_sum(r0) / 2.0)
        return [self.pack_theta(r0, base * a, base * b)
                for a, b in list(spreads)[:max(int(n_start), 1)]]

    def _em(self, th0, n_iter, mstep_maxiter, verbose=False):
        """单次 EM（内部用）：返回 ``(θ, hist, logL 序列)``。"""
        th = np.asarray(th0, float).ravel().copy()
        hist, logLs = [th.copy()], []
        for it in range(int(n_iter)):
            r = self.set_theta(th)                       # 同步 Φ/β 供 E 步使用
            fb = self.forward_backward(r)                # E 步（与 level1 相同）
            logLs.append(fb["logL"])
            stats = self.mstep_stats(fb["pi_bar"])       # 充分统计量只算一次
            th = self.mstep_theta(th, stats, maxiter=mstep_maxiter)   # M 步
            hist.append(th.copy())
            if verbose:
                rr, p1, p2, bg = self.unpack_theta(th)
                print(f"  it{it:02d} logL={fb['logL']:.1f} "
                      f"d={self.separation_nm(rr):6.2f} nm  "
                      f"Φ1={p1:7.1f} Φ2={p2:7.1f} β={bg:.3f}")
        return th, hist, logLs

    def fit(self, theta0=None, n_iter: int = 20, mstep_maxiter: int = 150,
            verbose: bool = False, r_init=None, multi_start: int = 3,
            phi_init=None):
        """EM: E 步（前向-后向）+ M 步（θ 联合最大化 Q），**多起点**取最优。

        返回 ``dict(theta, r, phi1, phi2, bg, logL, hist, d_est_nm,
        logL_margin, logL_all, theta_all)``。

        - ``logL`` 逐轮应单调不降（标准 EM 的保证，可直接当正确性检查）；
        - ``logL_margin`` = 最优末轮 logL − 次优末轮 logL，``starts`` 列出每个起点的
          终点 ``(Φ1, Φ2, d, logL)``；
        - **判读方式**：先看 ``starts`` 里各起点的 ``d`` 是否一致 —— **一致即稳健**
          （即使 Φ 不同，只要 $d$ 相同就是同一个解，余下只是标签/亮度比退化）；
          再看 ``logL_margin``：接近 0 且 $d$ 不一致 ⇒ 山谷内多点统计等价、结果不可信。
          注意 $d$ 一致时 ``logL_margin`` 也可能为 0，那是**良性**的（多起点收敛到同一点）。
        """
        if theta0 is not None:
            starts = [np.asarray(theta0, float).ravel()]
        else:
            starts = self.initial_thetas(r_init=r_init, n_start=multi_start,
                                         phi_init=phi_init)
        best, all_last, all_th = None, [], []
        for th0 in starts:
            th, hist, logLs = self._em(th0, n_iter, mstep_maxiter, verbose)
            last = logLs[-1] if logLs else -np.inf
            all_last.append(last)
            all_th.append(th)
            if best is None or last > best[1]:
                best = (th, last, hist, logLs)
        th, _, hist, logLs = best
        order = np.sort(np.asarray(all_last, float))[::-1]
        margin = float(order[0] - order[1]) if order.size > 1 else np.inf

        self.set_theta(th)
        self.theta_, self.logL_history_ = th, logLs
        r, p1, p2, bg = self.unpack_theta(th)
        self.r_ = r
        starts_rep = []
        for th_i in all_th:
            ri, q1, q2, _ = self.unpack_theta(th_i)
            starts_rep.append(dict(phi1=q1, phi2=q2,
                                   d_nm=self.separation_nm(ri),
                                   r=np.asarray(ri, float).copy()))
        return {"theta": th, "r": r, "phi1": p1, "phi2": p2, "bg": bg,
                "logL": logLs, "hist": hist,
                "d_est_nm": self.separation_nm(r),
                "logL_margin": margin, "logL_all": all_last, "theta_all": all_th,
                "n_start": len(starts), "starts": starts_rep}

    def __repr__(self):
        return (f"JointHMMTwoEmitter2(T={self.T}, ROI={self.h}x{self.w}, M={self.M}, "
                f"dt={self.dt}, Φ=({self.phi:.1f},{self.phi2:.1f}), "
                f"β={self.bg:.3f}, tie_phi={self.tie_phi}, fit_bg={self.fit_bg})")


# =========================================================================== #
# 便捷封装
# =========================================================================== #
def build_pair_hmm2(ds, pair_idx: int, **kw) -> JointHMMTwoEmitter2:
    """构造 level2 模型（默认用 ``psf_model='template'`` 消除 PSF 失配）。"""
    kw.setdefault("psf_model", "template")
    return build_pair_hmm(ds, pair_idx, cls=JointHMMTwoEmitter2, **kw)


def _align_pair(res_r, gt_xy):
    """返回与真值配对一致的 (下标, 交换?) —— 两分子不可区分, 取位置误差更小的配对。"""
    e_ab = np.hypot(res_r[0] - gt_xy[0, 0], res_r[1] - gt_xy[0, 1]) + \
        np.hypot(res_r[2] - gt_xy[1, 0], res_r[3] - gt_xy[1, 1])
    e_ba = np.hypot(res_r[2] - gt_xy[0, 0], res_r[3] - gt_xy[0, 1]) + \
        np.hypot(res_r[0] - gt_xy[1, 0], res_r[1] - gt_xy[1, 1])
    return (e_ab <= e_ba), 0.5 * min(e_ab, e_ba)


def fit_pair_level2(ds, pair_idx: int, half_px: int = 6, bin_frames: int = 1,
                    m_sub: int = 4, n_iter: int = 20, tie_phi: bool = False,
                    fit_bg: bool = False, psf_model: str = "template",
                    phi_init=None, bg=None, mstep_maxiter: int = 150,
                    multi_start: int = 3, verbose: bool = False, **kw) -> dict:
    """对一个 pair 跑 level2（位置 + 亮度 联合估计）。

    ``multi_start`` 默认 3：见模块 docstring —— (Φ1,Φ2,d) 有平坦山谷，
    必须多起点 + 取最高 logL。``tie_phi=True`` 时自动退化为单起点。
    """
    hmm = build_pair_hmm2(ds, pair_idx, half_px=half_px, bin_frames=bin_frames,
                          m_sub=m_sub, tie_phi=tie_phi, fit_bg=fit_bg,
                          psf_model=psf_model, bg=bg, **kw)
    out = hmm.fit(theta0=None, n_iter=n_iter, mstep_maxiter=mstep_maxiter,
                  verbose=verbose, multi_start=multi_start, phi_init=phi_init)

    idx = ds.pair_emitters(pair_idx)
    gt = ds.gt_xyz[idx, :2]
    same, e_pos = _align_pair(out["r"], gt)
    p_fit = (out["phi1"], out["phi2"]) if same else (out["phi2"], out["phi1"])
    p_true = tuple(float(ds.gt_brightness[k]) * ds.dt * int(bin_frames)
                   for k in idx)
    # 稳健性诊断：各起点终点 d 的散布（小 ⇒ 稳健）；标签对齐后各起点的 Φ
    d_st = np.asarray([s["d_nm"] for s in out["starts"]], float)
    ph_st = np.asarray([(s["phi1"], s["phi2"]) if same else
                        (s["phi2"], s["phi1"]) for s in out["starts"]], float)
    out.update(dict(pair_idx=pair_idx, hmm=hmm, gt=gt,
                    d_true_nm=ds.pair_gt_separation_nm(pair_idx),
                    emitter_err_mean_nm=e_pos,
                    phi_fit=p_fit, phi_true=p_true,
                    d_start_std_nm=float(np.std(d_st)),
                    phi_start_std=tuple(float(v) for v in np.std(ph_st, axis=0)),
                    phi_err_rel=tuple(float(f / t - 1.0) for f, t in
                                      zip(p_fit, p_true))))
    return out


def summarize_level2(datasets, half_px: int = 6, bin_frames: int = 1, m_sub: int = 4,
                     n_iter: int = 20, n_pairs: int | None = None,
                     tie_phi: bool = False, fit_bg: bool = False,
                     psf_model: str = "template", mstep_maxiter: int = 150,
                     multi_start: int = 3, **kw) -> list[dict]:
    """逐 pair 跑 level2 并汇总（每个数据集一行）。"""
    rows = []
    for ds in datasets:
        npair = ds.n_pairs if n_pairs is None else min(int(n_pairs), ds.n_pairs)
        d_est, d_true = [], []
        phi_f, phi_t, p_rel = [], [], []
        epos, dsp, mgn = [], [], []
        for p in range(npair):
            res = fit_pair_level2(ds, p, half_px=half_px, bin_frames=bin_frames,
                                  m_sub=m_sub, n_iter=n_iter, tie_phi=tie_phi,
                                  fit_bg=fit_bg, psf_model=psf_model,
                                  mstep_maxiter=mstep_maxiter,
                                  multi_start=multi_start, **kw)
            d_est.append(res["d_est_nm"])
            d_true.append(res["d_true_nm"])
            epos.append(res["emitter_err_mean_nm"])
            dsp.append(res["d_start_std_nm"])
            mgn.append(res["logL_margin"])
            phi_f.extend(res["phi_fit"])
            phi_t.extend(res["phi_true"])
            p_rel.extend(res["phi_err_rel"])
        d_est = np.asarray(d_est, float)
        d_true = np.asarray(d_true, float)
        rows.append({
            "separation_nm": ds.separation_nm,
            "n_pairs": npair,
            "bin_frames": int(bin_frames), "m_sub": int(m_sub),
            "d_est_mean_nm": float(np.mean(d_est)),
            "d_bias_nm": float(np.mean(d_est - d_true)),
            "d_rmse_nm": float(np.sqrt(np.mean((d_est - d_true) ** 2))),
            "d_start_std_nm": float(np.mean(dsp)),
            "logL_margin_min": float(np.min(mgn)),
            "emitter_err_mean_nm": float(np.mean(epos)),
            "phi_fit_mean": float(np.mean(phi_f)),
            "phi_true_mean": float(np.mean(phi_t)),
            "phi_bias_rel": float(np.mean(np.asarray(phi_f) / np.asarray(phi_t) - 1.0)),
            "phi_scatter_rel": float(np.std(np.asarray(phi_f) / np.asarray(phi_t) - 1.0)),
        })
    return rows
