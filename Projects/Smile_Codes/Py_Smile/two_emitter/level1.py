"""two_emitter.level1 —— 状态未知、转移矩阵已知: 两个双态发射体的联合推断。

模型(推导见 `test_two_emitter.ipynb` 第 7 节)
------------------------------------------------
单分子:  ``z^(i)_t ∈ {0=on, 1=off}``, 已知速率 ``k_on/k_off`` → 帧转移 ``P = exp(Q·dt)``。
联合:    ``Z_t = (z^(1)_t, z^(2)_t)`` 只有 4 个取值; 两分子独立 ⇒ ``Ã = P ⊗ P``。

★ 帧内跳变: 一帧并不"瞬时", 分子可能在帧中跳变。真正决定发射的是**帧内 on 占空**
   ``n_i/M`` (把一帧切成 M 个亚步, Δ=dt/M, 数出其中 on 的亚步数 n_i),
   而不是帧边界上的状态。给定端点 ``(a,b)`` 的占空条件分布

       ``q_{a_i b_i}(n_i) = P(n_i | 帧初状态 a_i, 帧末状态 b_i)``

   可以用一个 2 状态 × (M+1) 计数的小 DP 精确算出。于是

       ``p(D_t | a, b) = Σ_{n1,n2} q_{a1,b1}(n1)·q_{a2,b2}(n2)·p(D_t | n1, n2)``

   注意 ``p(D_t | ·, ·)`` 依赖**相邻两个**状态 → 这是"成对势"链,
   前向-后向仍是 T×4 的递归, 只需把 ``Ã_{ab}`` 换成 ``φ_t(a,b) = Ã_{ab}·p(D_t|a,b)``。
   取 ``M=1`` 就退化成常见的"以帧状态近似"的模型(且帧内占空被强制为 0 或 1),
   这也是 level0 里 ``states`` 判据出偏置的根源。

发射模型
--------
``μ_{t,p}(n1,n2) = β + (n1/M)·Φ₁·P₁(r₁)_p + (n₂/M)·Φ₂·P₂(r₂)_p``

``P_i`` 是 ROI 内**未归一化**的像素积分高斯(其和 = ROI 内圈入比例), ``Φ_i`` = 整帧全亮
时发射的光子数。噪声用 EMCCD 的"泊松 + 过噪因子 F"的高斯近似: ``var = F²·μ``。

推断
----
- ``log_likelihood(r)``   : 前向算法
- ``forward_backward(r)`` : 状态后验 γ/ξ + 每帧每个分子的后验 on 占空 ŵ
- ``fit(...)``            : EM (E=前向-后向, M=位置用残差图的 MLE 闭式更新)
- ``fit_gradient(...)``   : 直接对 −logL 做数值梯度优化(交叉验证)
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.special import erf, logsumexp

from .level0 import mle_localize_image

_SQRT2 = float(np.sqrt(2.0))


def _softmax_mean_rows(lw: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """对 ``lw`` (T, K) 的每一行做 softmax 后求 ``vals`` (K,) 的加权平均。

    稳健处理全 ``-inf`` 的行(返回 0), 避免 ``-inf - (-inf) = nan``。
    """
    finite = np.isfinite(lw)
    any_f = finite.any(axis=1, keepdims=True)
    mx = np.where(any_f, np.max(np.where(finite, lw, -np.inf), axis=1, keepdims=True), 0.0)
    w = np.where(finite, np.exp(lw - mx), 0.0)
    s = w.sum(axis=1, keepdims=True)
    num = (w * vals[None, :]).sum(axis=1, keepdims=True)
    return (num / np.maximum(s, 1e-300)).ravel()


# =========================================================================== #
# 单分子 2 态
# =========================================================================== #
def transition_2state(k_on: float, k_off: float, dt: float) -> np.ndarray:
    """2 态帧转移矩阵 ``P(dt)`` (状态顺序 ``[on, off]``)。"""
    Q = np.array([[-float(k_off), float(k_off)],
                  [float(k_on), -float(k_on)]], dtype=float)
    return expm(Q * float(dt))


def stationary_2state(k_on: float, k_off: float) -> np.ndarray:
    """2 态平稳分布 ``[P(on), P(off)]``。"""
    p_on = float(k_on) / (float(k_on) + float(k_off))
    return np.array([p_on, 1.0 - p_on])


def oncount_given_endpoints(P_sub: np.ndarray, M: int,
                            on_index: int = 0) -> np.ndarray:
    """``q[a, b, n]`` = P(帧内 on 亚步数 = n, 帧末状态 = b | 帧初状态 = a)。

    约定: 状态只在亚步边界更新; 亚步 ``k`` (k=1..M) 期间的发射由该亚步**开始时**的
    状态决定。因此 ``Σ_{b,n} q[a,b,n] = 1`` 且 ``Σ_n q[a,b,n] = P_sub[a,b]``。
    """
    P = np.asarray(P_sub, dtype=float)
    ns = P.shape[0]
    z = np.zeros(ns)
    z[on_index] = 1.0
    q = np.zeros((ns, ns, M + 1))
    for a in range(ns):
        al = np.zeros((ns, M + 1))
        al[a, 0] = 1.0
        for _ in range(M):
            new = np.zeros((ns, M + 1))
            for s in range(ns):
                inc = int(z[s])
                new[:, inc:] += np.outer(P[s], al[s, : M + 1 - inc])
            al = new
        q[a] = al
    return q


# =========================================================================== #
# 数据准备
# =========================================================================== #
def pair_roi_photons(ds, pair_idx: int, half_px: int = 6, bin_frames: int = 1,
                     photons: bool = True):
    """取某个 pair 的 ROI 图像 (T', h, w)。

    ``bin_frames > 1`` 时把连续 ``bin_frames`` 帧**相加**(等效于降低帧率、曝光×m),
    用来研究"帧内跳变"的影响 —— 同一份数据, 只改分析分辨率。
    返回 ``(roi_data, (r0, c0))``, ``roi_data`` 为光子数(float)。
    """
    data = ds.photons() if photons else ds.images.astype(np.float32)
    roi = ds.roi(pair_idx, half_px=half_px)
    sub = np.asarray(data[:, roi[0], roi[1]], dtype=float)      # (T, h, w)
    if bin_frames > 1:
        m = int(bin_frames)
        T = sub.shape[0] // m
        sub = sub[: T * m].reshape(T, m, *sub.shape[1:]).sum(axis=1)
    return sub, (roi[0].start, roi[1].start)


# =========================================================================== #
# 4 态联合 HMM
# =========================================================================== #
class JointHMMTwoEmitter:
    """两个双态发射体的 4 态联合 HMM(状态未知、转移矩阵已知)。

    参数
    ----
    data : (T, h, w)
        ROI 内的光子数图像(每帧曝光 ``dt``)。
    pixel_size_nm, sigma_nm : float
        像素尺寸与 PSF 高斯等效 σ。
    phi_photons : float
        单个发射体"整帧全亮"发射的光子数 (``= brightness × dt``)。
    dt : float
        帧曝光时间 (s) —— 注意这是**分析**分辨率下的, 若做了帧合并应传入合并后的。
    k_on, k_off : float
        已知的 2 态速率 (1/s)。
    m_sub : int
        每帧切分的亚步数 M。``M=1`` = 帧分辨率近似; M 越大帧内占空建模越精确。
    bg : float, 可选
        每帧每像素背景(光子)。None 时用 ROI 中位数初始化。
    F_EM : float
        EMCCD 过噪因子(默认 √2)。
    roi_origin : (r0, c0)
        ROI 在全视场中的起点, 用于把位置换算回全视场 nm 坐标。
    """

    #: 4 个联合状态 ``(z1, z2)``。约定与 :func:`transition_2state` 一致:
    #: **z = 0 表示 on, z = 1 表示 off** (``on_index=0``)。
    #: 于是联合索引 ``s = 2*z1 + z2`` 对应
    #: ``0=(0,0) 都亮, 1=(0,1) 只分子1亮, 2=(1,0) 只分子2亮, 3=(1,1) 都暗``。
    #: 亮度由 ``q[a>>1, b>>1, n1] · q[a&1, b&1, n2]`` 给出, 即 ``z_i=0`` 一侧取 ``q`` 的
    #: 第 0 行(其 ``Σ_n q[0,b,n] = P_sub[0,b]``, 而 ``P_sub[0,0]=0.909`` = 保持 on)。
    #: (实测校验: 真值"都亮"的帧上 ``γ[:,0]≈0.80``, "都暗"的帧上 ``γ[:,3]≈0.97``。)
    STATES = ((0, 0), (0, 1), (1, 0), (1, 1))

    def __init__(self, data, pixel_size_nm: float, sigma_nm: float,
                 phi_photons: float, dt: float, k_on: float, k_off: float,
                 m_sub: int = 1, bg: float | None = None,
                 F_EM: float = _SQRT2, roi_origin=(0, 0),
                 noise: str = "poisson", psf_model: str = "airy",
                 na: float | None = None, wavelength_nm: float | None = None,
                 psf_kernel=None):
        self.data = np.asarray(data, dtype=float)
        if self.data.ndim != 3:
            raise ValueError("data 应为 (T, h, w)")
        if noise not in ("poisson", "gauss"):
            raise ValueError("noise 只能是 'poisson' 或 'gauss'")
        if psf_model not in ("airy", "gauss", "template"):
            raise ValueError("psf_model 只能是 'airy' / 'gauss' / 'template'")
        if psf_model == "template" and psf_kernel is None:
            raise ValueError("psf_model='template' 需要 psf_kernel=(psf_obj, K)"
                             "（用 build_pair_hmm 会自动接好）")
        self.noise = noise
        self.psf_model = psf_model
        self.psf_kernel = psf_kernel            # ('template' 模式) (psf_obj, K)
        self.na = None if na is None else float(na)
        self.wavelength_nm = None if wavelength_nm is None else float(wavelength_nm)
        self.T, self.h, self.w = self.data.shape
        self.pixel_size_nm = float(pixel_size_nm)
        self.sigma_px = float(sigma_nm) / self.pixel_size_nm
        self.phi = float(phi_photons)           # 分子1的"整帧全亮"光子数
        self.phi2 = self.phi                    # 分子2; level1 恒相等, level2 可独立估计
        self.dt = float(dt)
        self.k_on, self.k_off = float(k_on), float(k_off)
        self.M = max(int(m_sub), 1)
        self.F_EM = float(F_EM)
        self.roi_origin = (int(roi_origin[0]), int(roi_origin[1]))
        # ---- 转移: 单分子 + 联合 ----
        self.P_frame = transition_2state(self.k_on, self.k_off, self.dt)
        self.P_sub = transition_2state(self.k_on, self.k_off, self.dt / self.M)
        self.prior1 = stationary_2state(self.k_on, self.k_off)
        self.Atilde = np.kron(self.P_frame, self.P_frame)          # (4,4)
        self.prior_joint = np.kron(self.prior1, self.prior1)       # (4,)
        self.log_Atilde = np.log(self.Atilde)

        # ---- 帧内占空分布 q[a,b,n] ----
        q = oncount_given_endpoints(self.P_sub, self.M)
        self.q = q
        with np.errstate(divide="ignore"):
            self.log_q = np.where(q > 0, np.log(np.where(q > 0, q, 1.0)), -np.inf)

        # ---- (n1,n2) 扁平网格 (K = (M+1)²) ----
        _g1, _g2 = np.meshgrid(np.arange(self.M + 1), np.arange(self.M + 1),
                               indexing="ij")
        self._sub_grid = (_g1.ravel().astype(float), _g2.ravel().astype(float))
        self._rtab = None                       # Airy 径向查表缓存
        self._airy_pad = 6                      # Airy 像素积分时的外扩像素数

        # ---- 背景/缓存 ----
        self.bg = float(bg) if bg is not None else self.estimate_bg()
        self._D2 = self.data ** 2
        self._sumD = self.data.sum(axis=(1, 2))                    # (T,)
        self._npix = self.h * self.w
        self.logL_history_: list[float] = []

    # ------------------------------------------------------------------ 背景
    def bg_mask(self):
        """ROI 外圈(离中心 >= (min(h,w)-1)/2 - 1 像素)的像素掩码 —— 基本无 PSF。"""
        yy, xx = np.mgrid[0:self.h, 0:self.w]
        cy, cx = (self.h - 1) / 2.0, (self.w - 1) / 2.0
        r = np.maximum(np.abs(yy - cy), np.abs(xx - cx))
        return r >= (min(self.h, self.w) - 1) / 2.0 - 1.0

    def estimate_bg(self, mask=None) -> float:
        """用 ROI 外圈像素的均值估计每帧每像素背景 (光子)。"""
        m = self.bg_mask() if mask is None else np.asarray(mask, bool)
        if not np.any(m):
            return max(float(np.median(self.data)), 1e-6)
        return max(float(self.data[:, m].mean()), 1e-6)

    # ------------------------------------------------------------------ PSF
    def _cdf(self, u):
        return 0.5 * (1.0 + erf(np.asarray(u, dtype=float) / _SQRT2))

    def _gauss_pixels(self, r_nm) -> np.ndarray:
        """ROI 内**未归一化**的像素积分高斯 (其和 = ROI 内的圈入比例)。"""
        s = self.sigma_px
        ux = float(r_nm[0]) / self.pixel_size_nm - self.roi_origin[1]   # 列方向
        uy = float(r_nm[1]) / self.pixel_size_nm - self.roi_origin[0]   # 行方向
        xc = np.arange(self.w) + 0.5
        yc = np.arange(self.h) + 0.5
        Ex = self._cdf((xc + 0.5 - ux) / s) - self._cdf((xc - 0.5 - ux) / s)
        Ey = self._cdf((yc + 0.5 - uy) / s) - self._cdf((yc - 0.5 - uy) / s)
        return Ey[:, None] * Ex[None, :]

    def _radial_table(self, n: int = 4000):
        """Airy 径向强度 ``I(r) = [2J₁(x)/x]²`` 的一维查表 (构造一次, 复用)。"""
        if self._rtab is None:
            if self.na is None or self.wavelength_nm is None:
                raise ValueError("psf_model='airy' 需要同时给出 na 与 wavelength_nm")
            from scipy.special import j1
            pad = self._airy_pad
            r_max = np.hypot(self.h + 2 * pad, self.w + 2 * pad) * self.pixel_size_nm
            rr = np.linspace(0.0, r_max, int(n))
            x = 2.0 * np.pi * self.na * rr / self.wavelength_nm
            safe = np.where(x > 1e-9, x, 1.0)
            self._rtab = (rr, np.where(x > 1e-9, (2.0 * j1(safe) / safe) ** 2, 1.0))
        return self._rtab

    def _airy_pixels(self, r_nm, sub: int = 4, pad: int | None = None) -> np.ndarray:
        """ROI 内像素积分的**圆孔径 Airy** 斑 (未归一化, 和 = ROI 内圈入比例)。

        用 ``sub×sub`` 子样在放大后的窗口上积分并按总能量归一 (再裁回 ROI),
        使发射模型与仿真里真实的 Airy PSF 一致。径向强度走一维查表 + 线性插值,
        避免每次调用都算 Bessel 函数。
        """
        pad = self._airy_pad if pad is None else int(pad)
        px = self.pixel_size_nm
        H, W = self.h + 2 * pad, self.w + 2 * pad
        xs = (np.arange(W * sub) + 0.5) / sub - pad        # 局部像素坐标(可负)
        ys = (np.arange(H * sub) + 0.5) / sub - pad
        ux = float(r_nm[0]) / px - self.roi_origin[1]
        uy = float(r_nm[1]) / px - self.roi_origin[0]
        dx = (xs - ux) * px
        dy = (ys - uy) * px
        R = np.hypot(dy[:, None], dx[None, :])
        rr, I = self._radial_table()
        prof = np.interp(R, rr, I)
        prof = prof.reshape(H, sub, W, sub).mean(axis=(1, 3))      # 像素积分
        tot = prof.sum()
        if tot > 0:
            prof = prof / tot
        return prof[pad:pad + self.h, pad:pad + self.w]

    def psf_pixels(self, r_nm) -> np.ndarray:
        """归一化到"总能量=1"的像素化 PSF, 截到 ROI; 其和 = ROI 内圈入比例。

        ``psf_model`` 三选一:

        - ``'airy'``     解析圆孔径 Airy —— 与仿真器的 FFT 核有 0.53% RMS 形状差
                         （表现为模型 PSF 偏宽 2–3%），是残余间距偏差地板的主因；
        - ``'gauss'``    像素积分高斯 —— 仅作对照；
        - ``'template'`` **仿真器自己的过采样核**，逐点一致、无失配（见 ``two_emitter.psf``）。
        """
        if self.psf_model == "template":
            from .psf import roi_template
            psf_obj, K = self.psf_kernel
            return roi_template(psf_obj, K, r_nm, self.roi_origin,
                                self.h, self.w, self.pixel_size_nm)
        if self.psf_model == "airy":
            return self._airy_pixels(r_nm)
        return self._gauss_pixels(r_nm)

    # ------------------------------------------------- 帧发射似然（对数）
    def loglik_table(self, r) -> np.ndarray:
        """``logE[t, n1, n2] = log p(D_t | n1, n2, r)``。

        ``noise='poisson'`` (默认)
            泊松对数似然 ``Σ_p [D_p·log μ_p − μ_p]``。EMCCD 的过噪只影响方差,
            不影响均值, 因此泊松得分函数仍然无偏 —— 用它做 MLE 是**一致**的,
            而且不会引入伪吸引 (见下)。
        ``noise='gauss'``
            高斯近似 ``Σ_p [-(D−μ)²/(F²μ) − log(F²μ)]``(与 ``smlm.localize``
            的 MLELocalizer 一致)。注意其中的 ``−½Σlog μ`` 项会**偏好峰值更高**的
            构型, 当两个 PSF 重叠时等于给"间距变小"一个额外的奖励 —— 在 σ 不可调
            (固定)的双体联合拟合里会带来 ~10% 的间距低估, 故只作对照。
        """
        r = np.asarray(r, dtype=float).ravel()
        P1 = self.psf_pixels(r[0:2])
        P2 = self.psf_pixels(r[2:4])
        M = self.M
        K = (M + 1) * (M + 1)
        n1s, n2s = np.meshgrid(np.arange(M + 1), np.arange(M + 1),
                               indexing="ij")
        n1s = n1s.ravel() / M
        n2s = n2s.ravel() / M
        p1 = P1.ravel()
        p2 = P2.ravel()

        # --- 只算 K 个 (n1,n2) 的 log μ / 1/μ / Σμ; 然后一次 gemm 得到 (T,K) ---
        logmu = np.empty((K, p1.size))
        invmu = np.empty((K, p1.size))
        smu = np.empty(K)
        for k in range(K):
            mu = self.bg + n1s[k] * self.phi * p1 + n2s[k] * self.phi2 * p2
            logmu[k] = np.log(mu)
            invmu[k] = 1.0 / mu
            smu[k] = mu.sum()

        D = self.data.reshape(self.T, -1)
        out = (D @ logmu.T) - smu[None, :]                 # Σ D log μ − Σ μ

        if self.noise == "gauss":
            F2 = self.F_EM ** 2
            S2 = (D ** 2) @ invmu.T                        # Σ D²/μ
            out = (-0.5 / F2 * (S2 - 2.0 * self._sumD[:, None] + smu[None, :])
                   - 0.5 * logmu.sum(axis=1)[None, :]
                   - 0.5 * self._npix * np.log(F2))
        return out.reshape(self.T, M + 1, M + 1)

    # ------------------------------------------------- 成对势 φ_t(a,b)
    def _phi_log(self, r, want_ratios: bool = False):
        """返回 ``(logphi, ratios, logE)``。

        ★ 注意: 单分子占空分布 ``q[a,b,n]`` 已经包含了帧转移概率
        (``Σ_n q[a,b,n] = P(dt)_{ab}``), 两个分子独立 ⇒ 乘积即
        ``Ã_{ab}·P(n1,n2 | a,b)``。因此这里 **不能再乘一次 ``Ã_{ab}``**, 否则
        转移概率被重复计入(会把后验压平、EM 发散)。

        ``logphi[t,a,b] = log Σ_{n1,n2} q[a1,b1,n1]·q[a2,b2,n2]·p(D_t | n1,n2)``。

        实现上利用**可分离性**降低复杂度(先各自沿 n1 / n2 归约一次, 再组合),
        比逐 (a,b) 对整张 (T,M+1,M+1) 做 logsumexp 快约一个数量级:
          ``B[a1,b1](n2) = logsumexp_{n1}[logq₁ + logE]``
          ``C[a2,b2](n1) = logsumexp_{n2}[logq₂ + logE]``
        ``ratios`` = (m1, m2), 每个 (T,4,4): 给定 (a,b) 的后验平均占空。
        """
        logE = self.loglik_table(r)
        T, M = self.T, self.M
        vals = np.arange(M + 1) / M                       # (M+1,)
        logphi = np.empty((T, 4, 4))
        m1 = np.zeros((T, 4, 4))
        m2 = np.zeros((T, 4, 4))

        # 先各自沿一个轴归约 (n1: axis=1, n2: axis=2)
        B = {(a1, b1): logsumexp(logE + self.log_q[a1, b1][None, :, None], axis=1)
             for a1 in (0, 1) for b1 in (0, 1)}
        C = {(a2, b2): logsumexp(logE + self.log_q[a2, b2][None, None, :], axis=2)
             for a2 in (0, 1) for b2 in (0, 1)}

        lq = self.log_q
        for a in range(4):
            a1, a2 = a >> 1, a & 1
            for b in range(4):
                b1, b2 = b >> 1, b & 1
                BB = B[(a1, b1)] + lq[a2, b2][None, :]     # (T, M+1) 对 n2
                logphi[:, a, b] = logsumexp(BB, axis=1)
                if want_ratios:
                    m2[:, a, b] = _softmax_mean_rows(BB, vals)
                    m1[:, a, b] = _softmax_mean_rows(
                        C[(a2, b2)] + lq[a1, b1][None, :], vals)
        return logphi, (m1, m2), logE

    # ------------------------------------------------- 前向 / 前向-后向
    def _alpha_beta(self, r, want_ratios: bool = False):
        """前向-后向的公共实现。

        索引约定(frames ``D_1..D_T`` ↔ 0-based ``logphi[0..T-1]``; 状态 ``Z_0..Z_T``):

        - ``logalpha[k]`` = ``log p(D_1..D_k, Z_k)`` —— 长度 **T+1**, ``logalpha[0]=log π``;
        - ``logbeta[k]``  = ``log p(D_{k+1}..D_T | Z_k)`` —— 长度 **T+1**, ``logbeta[T]=0``;
        - ``logL = logsumexp(logalpha[T])``;
        - ``γ_k`` (长度 T+1) = P(Z_k | D);  ``ξ_t`` (长度 T) = P(Z_t, Z_{t+1} | D) ——
          帧 t 的发射由 ``(Z_t, Z_{t+1})`` 决定, 正好与 ``logphi[t]`` 对齐。
        """
        logphi, ratios, logE = self._phi_log(r, want_ratios=want_ratios)
        logpi = np.log(self.prior_joint)
        T = self.T

        # ---- 前向 / 后向递归 (4 状态, 手写 logsumexp 比 scipy 快 ~10 倍) ----
        logalpha = np.full((T + 1, 4), -np.inf)
        logalpha[0] = logpi
        exp, log = np.exp, np.log
        for k in range(1, T + 1):
            z = logalpha[k - 1][:, None] + logphi[k - 1]
            m = z.max(axis=0)
            logalpha[k] = m + log(exp(z - m).sum(axis=0))
        logL = float(logsumexp(logalpha[T]))

        logbeta = np.full((T + 1, 4), -np.inf)
        logbeta[T] = 0.0
        for k in range(T - 1, -1, -1):
            z = logphi[k] + logbeta[k + 1][None, :]
            m = z.max(axis=1)
            logbeta[k] = m + log(exp(z - m[:, None]).sum(axis=1))

        loggamma = logalpha + logbeta - logL                    # (T+1, 4)
        gamma = np.exp(loggamma)
        logxi = (logalpha[:-1][:, :, None] + logphi
                 + logbeta[1:][:, None, :] - logL)              # (T, 4, 4)
        xi = np.exp(logxi)
        return {"logL": logL, "gamma": gamma, "xi": xi,
                "logalpha": logalpha, "logbeta": logbeta,
                "loggamma": loggamma, "ratios": ratios,
                "logE": logE, "logphi": logphi}

    def _pi_bar(self, logalpha, logbeta, logL, logE) -> np.ndarray:
        """帧内占空的后验 ``π̄_t(n1,n2) = P(n_t = n | D)`` —— M 步的权重。

        直接对**联合** ``(Z_t, Z_{t+1}, n1, n2)`` 做 softmax:

            ``log P(a,b,n|D) ∝ logα_{t-1}(a) + log q₁(n₁) + log q₂(n₂)
                                + log p(D_t|n) + log β_t(b) − log L``

        注意 ``q[a,b,n] = P(n, Z_{t+1}=b | Z_t=a)`` **已经含了转移概率** ``Ã_{ab}``,
        所以这里不再乘 ``Ã``(否则转移概率被重复计入 —— 与 ``_phi_log`` 同一个坑)。

        再对 ``(a,b)`` 用 ``logaddexp`` 求和得到 ``log π̄(n)``; 归一化后
        ``Σ_n π̄(n) = 1`` 自动成立。

        为什么不用 ``ξ_t(a,b)/p(D_t|a,b)`` 的代数化简: ``q`` 在 M 较大时动态范围
        极大, 化简式里 ``q₁q₂`` 的连乘很容易下溢; 这里的写法只做加法与
        ``logaddexp``, 数值上稳得多。
        """
        T, M = self.T, self.M
        lq = self.log_q
        acc = np.full((T, M + 1, M + 1), -np.inf)
        for a in range(4):
            la = logalpha[:-1, a]
            for b in range(4):
                base = la + logbeta[1:, b] - logL        # 不含 Ã: 它在 q 里
                lw = (base[:, None, None]
                      + lq[a >> 1, b >> 1][None, :, None]     # 分子 1 的占空 q₁(n₁)
                      + lq[a & 1, b & 1][None, None, :]       # 分子 2 的占空 q₂(n₂)
                      + logE)
                np.logaddexp(acc, lw, out=acc)
        mx = acc.max(axis=(1, 2), keepdims=True)
        W = np.exp(acc - mx)
        s = W.sum(axis=(1, 2), keepdims=True)
        return W / np.maximum(s, 1e-300)

    def forward(self, r):
        """前向算法 (log 域)。返回 ``(logL, logalpha (T+1,4))``。"""
        out = self._alpha_beta(r, want_ratios=False)
        return out["logL"], out["logalpha"]

    def log_likelihood(self, r) -> float:
        """``log p(D_{1:T} | r, A)`` (前向算法)。"""
        return self._alpha_beta(r, want_ratios=False)["logL"]

    def forward_backward(self, r):
        """前向-后向 (E 步)。返回 dict:

        - ``logL``            : 完整序列似然;
        - ``gamma[k]`` (T+1,4): ``P(Z_k = · | D)``;
        - ``xi[t]`` (T,4,4)   : ``P(Z_t = ·, Z_{t+1} = · | D)`` (与帧 t 对齐);
        - ``wfrac[t,i]`` (T,2): 第 t 帧第 i 个分子的后验 on 占空 ``E[n_i/M | D]``;
        - ``pi_bar[t]`` (T,M+1,M+1): 帧内占空的后验 ``π̄_t(n1,n2)``, M 步要用。
        """
        out = self._alpha_beta(r, want_ratios=True)
        m1, m2 = out["ratios"]
        xi = out["xi"]
        wfrac = np.stack([(xi * m1).sum(axis=(1, 2)),
                          (xi * m2).sum(axis=(1, 2))], axis=1)
        pi_bar = self._pi_bar(out["logalpha"], out["logbeta"], out["logL"],
                              out["logE"])
        return {"logL": out["logL"], "gamma": out["gamma"], "xi": xi,
                "wfrac": wfrac, "pi_bar": pi_bar,
                "logalpha": out["logalpha"], "logbeta": out["logbeta"]}

    # ------------------------------------------------- M 步 (最大化 Q)
    def mstep_stats(self, pi_bar):
        """M 步的**充分统计量**: ``Ď_{n,p} = Σ_t π̄_{t,n}·D_{t,p}``。

        有了它就无需在每次 M 步试位时重新扫全部帧 ——
        ``Q(r) = Σ_{n,p}[Ď_{n,p}log μ_{n,p}] − Σ_n m_n Σ_p μ_{n,p}``,
        每次试位只花在 (K, P) 的小矩阵上, 比逐帧重算快 ~10 倍。
        """
        flat = pi_bar.reshape(self.T, -1)                 # (T, K)
        D = self.data.reshape(self.T, -1)                 # (T, P)
        Dbar = flat.T @ D                                 # (K, P)
        m_n = flat.sum(axis=0)                            # (K,)
        if self.noise == "gauss":
            D2bar = flat.T @ (D ** 2)
            return Dbar, m_n, D2bar, self._sumD
        return Dbar, m_n

    def _mu_matrix(self, r):
        """(K, P) 的期望像素强度矩阵, K = (M+1)²。"""
        P1 = self.psf_pixels(r[0:2])
        P2 = self.psf_pixels(r[2:4])
        M = self.M
        n1s, n2s = self._sub_grid
        MU = (self.bg
              + np.outer(n1s * (self.phi / M), P1.ravel())
              + np.outer(n2s * (self.phi2 / M), P2.ravel()))
        return MU

    def q_value_fast(self, r, stats) -> float:
        """用充分统计量计算 Q(r) (比 ``q_value`` 快约一个数量级)。"""
        MU = self._mu_matrix(r)
        if self.noise == "poisson":
            Dbar, m_n = stats
            return float((Dbar * np.log(MU)).sum()
                         - (m_n * MU.sum(axis=1)).sum())
        Dbar, m_n, D2bar, sumD = stats
        F2 = self.F_EM ** 2
        # Q = −(1/2F²)Σ Ď/μ + (1/F²)ΣD − (1/2F²)Σ m_n Σμ − (P·T/2)log F²
        #     − ½Σ m_n Σ log μ      (其中 (1/F²)ΣD 与 r 无关, 只影响 Q 的绝对值)
        return float((-0.5 / F2 * (D2bar / MU).sum()
                      + sumD.sum() / F2
                      - 0.5 / F2 * (m_n * MU.sum(axis=1)).sum()
                      - 0.5 * (m_n * np.log(F2 * MU).sum(axis=1)).sum()))

    def q_value(self, r, pi_bar) -> float:
        """期望完整数据对数似然 ``Q(r) = Σ_t Σ_n π̄_t(n)·log p(D_t | n, r)``。

        (等价于 :meth:`q_value_fast` 配上 ``mstep_stats(pi_bar)``, 但慢一些;
        留作自查。)
        """
        return self.q_value_fast(r, self.mstep_stats(pi_bar))

    def mstep_q(self, r, pi_bar=None, stats=None, maxiter: int = 80,
                xtol: float = 0.05, method: str = "Nelder-Mead",
                polish: bool = True):
        """M 步: ``argmax_r Q(r; π̄)`` —— 保证 EM 单调不降 ``logL``。"""
        from scipy.optimize import minimize
        if stats is None:
            stats = self.mstep_stats(pi_bar)
        f = lambda rr: -self.q_value_fast(rr, stats)

        def run(x0, it, simp=None):
            opts = (dict(maxiter=int(it), xatol=xtol, fatol=1e-3)
                    if method == "Nelder-Mead" else dict(maxiter=int(it)))
            if simp is not None and method == "Nelder-Mead":
                opts["initial_simplex"] = simp
            return minimize(f, np.asarray(x0, float), method=method, options=opts)

        res = run(r, maxiter)
        if polish and method == "Nelder-Mead":
            x = np.asarray(res.x, float)
            scale = np.maximum(np.abs(x) * 1e-4, 0.5)
            simp = np.vstack([x, x + np.diag(scale)])
            res2 = run(x, max(20, maxiter // 2), simp=simp)
            if res2.fun <= res.fun:
                res = res2
        return np.asarray(res.x, dtype=float)

    def _mstep_positions(self, r, wfrac, half_px=None):
        """(备用) 残差图 + 高斯 MLE 的启发式更新 —— 快但有模型失配偏置。"""
        P1 = self.psf_pixels(r[0:2])
        P2 = self.psf_pixels(r[2:4])
        mu1 = wfrac[:, 0][:, None, None] * self.phi * P1[None]     # (T,h,w)
        mu2 = wfrac[:, 1][:, None, None] * self.phi2 * P2[None]

        S1 = (self.data - self.bg - mu2).sum(axis=0)
        S2 = (self.data - self.bg - mu1).sum(axis=0)

        out = np.array(r, dtype=float)
        for k, S in enumerate((S1, S2)):
            est = mle_localize_image(np.clip(S, 0.0, None), self.pixel_size_nm,
                                     F_EM=self.F_EM)
            if est is None:
                continue
            out[2 * k] = est["x_nm"] + self.roi_origin[1] * self.pixel_size_nm
            out[2 * k + 1] = est["y_nm"] + self.roi_origin[0] * self.pixel_size_nm
        return out

    def _mstep_params(self, r, wfrac, fit_phi=True, fit_bg=True):
        """闭式更新背景 β 与亮度 Φ(加权最小二乘)。"""
        P1 = self.psf_pixels(r[0:2])
        P2 = self.psf_pixels(r[2:4])
        mu1 = wfrac[:, 0][:, None, None] * P1[None]               # (T,h,w)
        mu2 = wfrac[:, 1][:, None, None] * P2[None]
        if fit_phi:
            num1 = (wfrac[:, 0][:, None, None] * P1[None]
                    * (self.data - self.bg - self.phi * mu2)).sum()
            den1 = ((wfrac[:, 0] ** 2)[:, None, None] * P1[None] ** 2).sum()
            num2 = (wfrac[:, 1][:, None, None] * P2[None]
                    * (self.data - self.bg - self.phi * mu1)).sum()
            den2 = ((wfrac[:, 1] ** 2)[:, None, None] * P2[None] ** 2).sum()
            est = [float(num1 / den1) if den1 > 0 else None,
                   float(num2 / den2) if den2 > 0 else None]
            est = [e for e in est if e is not None and np.isfinite(e)]
            if est:
                self.phi = max(float(np.mean(est)), 1e-6)
        if fit_bg:
            resid = self.data - self.phi * (mu1 + mu2)
            self.bg = max(float(resid.mean()), 1e-6)

    # ------------------------------------------------- 初值
    def summed_signal(self) -> np.ndarray:
        """时间积分信号图 ``Σ_t D_t − T·β`` (clip 到非负) —— 已扣背景。"""
        return np.clip(self.data.sum(axis=0) - self.bg * self.T, 0.0, None)

    def centroid_nm(self, img=None) -> np.ndarray:
        """时间积分图的强度质心 (全视场 nm)。

        ★ 必须**先扣掉背景**: 背景在 ROI 内是均匀的, 但 ROI 窗口相对发射体并不
        严格对称, 未扣背景时会有一个把质心拉向窗口中心的偏置(可到十几 nm)。
        """
        img = self.summed_signal() if img is None else np.asarray(img, float)
        img = np.clip(img, 0.0, None)
        tot = img.sum()
        if tot <= 0:
            return np.array([np.nan, np.nan])
        yc = np.arange(self.h) + 0.5
        xc = np.arange(self.w) + 0.5
        cy = float((img.sum(axis=1) * yc).sum() / tot)
        cx = float((img.sum(axis=0) * xc).sum() / tot)
        return np.array([(cx + self.roi_origin[1]) * self.pixel_size_nm,
                         (cy + self.roi_origin[0]) * self.pixel_size_nm])

    def principal_axis(self):
        """时间积分图的协方差主轴 (单位向量) 与标准差 (nm)。

        两个发射体在积分图上表现为一个沿连线拉长的斑 -> 主轴方向就是 pair 的取向,
        因此可以用它把 (间距 × 取向) 的初值搜索从 2D 网格压到很小的范围。
        """
        img = self.summed_signal()
        tot = img.sum()
        yc = np.arange(self.h) + 0.5
        xc = np.arange(self.w) + 0.5
        cy = float((img.sum(axis=1) * yc).sum() / tot)
        cx = float((img.sum(axis=0) * xc).sum() / tot)
        dy = (yc - cy)[:, None] * self.pixel_size_nm
        dx = (xc - cx)[None, :] * self.pixel_size_nm
        vyy = float((img * dy ** 2).sum() / tot)
        vxx = float((img * dx ** 2).sum() / tot)
        vxy = float((img * dy * dx).sum() / tot)
        evals, evecs = np.linalg.eigh(np.array([[vyy, vxy], [vxy, vxx]]))
        v = evecs[:, -1]                     # (dy, dx) 方向
        return np.array([v[1], v[0]]), float(np.sqrt(max(evals[-1], 0.0)))

    def initial_positions(self, d_list=(0.0, 20.0, 40.0, 60.0, 90.0, 130.0, 180.0, 240.0),
                          ang_off_deg=(-45.0, -22.5, 0.0, 22.5, 45.0),
                          r0=None) -> np.ndarray:
        """初值: 由积分图主轴 + 粗格点搜索 (间距 × 主轴附近取向) 选似然最大者。"""
        r0 = self.centroid_nm() if r0 is None else np.asarray(r0, float)
        axis, _ = self.principal_axis()
        th0 = float(np.arctan2(axis[1], axis[0]))
        best, best_L = None, -np.inf
        for d in d_list:
            for off in ang_off_deg:
                th = th0 + np.deg2rad(off)
                h = 0.5 * float(d)
                r = np.array([r0[0] - h * np.cos(th), r0[1] - h * np.sin(th),
                              r0[0] + h * np.cos(th), r0[1] + h * np.sin(th)])
                L = self.log_likelihood(r)
                if L > best_L:
                    best_L, best = L, r
        return best

    # ------------------------------------------------- EM
    def fit(self, r_init=None, n_iter: int = 20, mstep: str = "q",
            fit_phi: bool = False, fit_bg: bool = False,
            mstep_maxiter: int = 80, verbose: bool = False):
        """EM: E 步 = 前向-后向求状态后验; M 步更新连续位置。

        ``mstep='q'`` (默认)
            最大化期望完整数据对数似然 ``Q(r; π̄)`` —— 标准 EM, **保证 logL 单调不降**;
        ``mstep='residual'``
            残差图 + 高斯 MLE 的启发式更新 —— 快, 但有 PSF 模型失配偏置(对照用)。

        返回 ``dict(r=估计位置(nm), logL=逐轮似然, hist=逐轮位置)``。
        """
        r = np.asarray(r_init if r_init is not None else self.initial_positions(),
                       dtype=float)
        hist, logLs = [r.copy()], []
        for it in range(int(n_iter)):
            fb = self.forward_backward(r)
            logLs.append(fb["logL"])
            if fit_phi or fit_bg:
                self._mstep_params(r, fb["wfrac"], fit_phi=fit_phi, fit_bg=fit_bg)
            if mstep == "q":
                stats = self.mstep_stats(fb["pi_bar"])       # 充分统计量只算一次
                r = self.mstep_q(r, stats=stats, maxiter=mstep_maxiter)
            else:
                r = self._mstep_positions(r, fb["wfrac"])
            hist.append(r.copy())
            if verbose:
                print(f"  it{it:02d}  logL={fb['logL']:.1f}  "
                      f"d={np.hypot(r[0] - r[2], r[1] - r[3]):.2f} nm")
        self.logL_history_ = logLs
        self.r_ = r
        return {"r": r, "logL": logLs, "hist": hist}

    def fit_gradient(self, r_init=None, method: str = "Nelder-Mead",
                     options=None):
        """直接对 −logL 做数值优化(不依赖 EM) —— 用作交叉验证。"""
        from scipy.optimize import minimize
        r0 = np.asarray(r_init if r_init is not None else self.initial_positions(),
                        dtype=float)

        def neg(r):
            # 极轻的正则(nm^-2)避免数值上跑到离谱位置; 对 logL(~1e5) 完全可忽略
            return -self.log_likelihood(r) + 1e-8 * float(np.sum(r ** 2))

        res = minimize(neg, r0, method=method, options=options or {"maxiter": 800})
        self.r_ = np.asarray(res.x, dtype=float)
        return {"r": self.r_, "logL": float(-res.fun), "success": bool(res.success)}

    # ------------------------------------------------- 汇报
    def separation_nm(self, r=None) -> float:
        r = np.asarray(r if r is not None else self.r_, dtype=float)
        return float(np.hypot(r[0] - r[2], r[1] - r[3]))

    def __repr__(self):
        return (f"JointHMMTwoEmitter(T={self.T}, ROI={self.h}x{self.w}, M={self.M}, "
                f"dt={self.dt}, Φ={self.phi:.1f}, β={self.bg:.3f})")


# =========================================================================== #
# 便捷封装
# =========================================================================== #
def build_pair_hmm(ds, pair_idx: int, half_px: int = 6, bin_frames: int = 1,
                   m_sub: int = 1, k_on: float | None = None,
                   k_off: float | None = None, psf_model: str = "airy",
                   noise: str = "poisson", cls=None,
                   **kw) -> JointHMMTwoEmitter:
    """从数据集与 manifest 配置构造某个 pair 的 4 态联合 HMM。

    速率默认从 ``ds.meta['photophysics']['rate_matrix']`` 取 (2 态情形:
    ``k_off = K[0,1]``, ``k_on = K[1,0]``); ``psf_model='airy'`` 时用 meta 里的
    NA / λ 构造与仿真一致的 Airy PSF。
    """
    data, origin = pair_roi_photons(ds, pair_idx, half_px=half_px,
                                    bin_frames=bin_frames)
    meta = ds.meta
    if k_on is None or k_off is None:
        Km = np.asarray(meta["photophysics"]["rate_matrix"], dtype=float)
        k_off = float(Km[0, 1]) if k_off is None else k_off
        k_on = float(Km[1, 0]) if k_on is None else k_on
    dt = ds.dt * int(bin_frames)
    phi = float(meta.get("photons_per_on_frame", 1.0)) * int(bin_frames)
    psf_meta = meta.get("psf", {}) or {}
    sigma = float(psf_meta.get("psf_sigma_nm", meta.get("psf_sigma_nm", 100.0)))
    if psf_model == "template" and "psf_kernel" not in kw:
        # 用 dataset metadata 复现仿真器的 PSF 核（全局带缓存, FFT 只算一次）
        from .psf import get_kernel
        kw["psf_kernel"] = get_kernel(psf_meta)
    if cls is None:
        cls = JointHMMTwoEmitter
    return cls(data, pixel_size_nm=ds.pixel_size, sigma_nm=sigma,
               phi_photons=phi, dt=dt, k_on=k_on, k_off=k_off,
               m_sub=m_sub, roi_origin=origin,
               psf_model=psf_model, noise=noise,
               na=psf_meta.get("na"), wavelength_nm=psf_meta.get("wavelength"),
               **kw)


def fit_pair_level1(ds, pair_idx: int, half_px: int = 6, bin_frames: int = 1,
                    m_sub: int = 4, n_iter: int = 20, init=None,
                    mstep: str = "q", mstep_maxiter: int = 80, **kw) -> dict:
    """对一个 pair 跑 level1 联合推断。返回位置/间距/似然等。

    ``m_sub=4`` 是默认值: 实测 M≥2 已能修掉大部分"帧内跳变"带来的偏置,
    M=4 之后基本饱和, 而耗时随 M 增长(见 notebook 里的 M 扫描)。
    """
    hmm = build_pair_hmm(ds, pair_idx, half_px=half_px, bin_frames=bin_frames,
                         m_sub=m_sub, **kw)
    out = hmm.fit(r_init=init, n_iter=n_iter, mstep=mstep,
                  mstep_maxiter=mstep_maxiter)
    r = out["r"]
    idx = ds.pair_emitters(pair_idx)
    gt = ds.gt_xyz[idx, :2]
    # 两个发射体不可区分 -> 用两种标签配对里误差更小的那个
    e_ab = np.hypot(r[0] - gt[0, 0], r[1] - gt[0, 1]) + \
        np.hypot(r[2] - gt[1, 0], r[3] - gt[1, 1])
    e_ba = np.hypot(r[2] - gt[0, 0], r[3] - gt[0, 1]) + \
        np.hypot(r[0] - gt[1, 0], r[1] - gt[1, 1])
    return {
        "pair_idx": pair_idx, "hmm": hmm, "r": r,
        "d_est_nm": hmm.separation_nm(r),
        "d_true_nm": ds.pair_gt_separation_nm(pair_idx),
        "gt": gt, "logL": out["logL"], "r_hist": out["hist"],
        "emitter_err_mean_nm": float(0.5 * min(e_ab, e_ba)),
        "bg": hmm.bg, "phi": hmm.phi,
    }


def summarize_level1(datasets, half_px: int = 6, bin_frames: int = 1, m_sub: int = 4,
                     n_iter: int = 20, n_pairs: int | None = None,
                     mstep: str = "q", mstep_maxiter: int = 80, **kw) -> list[dict]:
    """对一组数据集逐 pair 跑 level1, 汇总间距误差 (每个数据集一行)。"""
    rows = []
    for ds in datasets:
        npair = ds.n_pairs if n_pairs is None else min(int(n_pairs), ds.n_pairs)
        d_est, d_true, err = [], [], []
        for p in range(npair):
            res = fit_pair_level1(ds, p, half_px=half_px, bin_frames=bin_frames,
                                  m_sub=m_sub, n_iter=n_iter, mstep=mstep,
                                  mstep_maxiter=mstep_maxiter, **kw)
            d_est.append(res["d_est_nm"])
            d_true.append(res["d_true_nm"])
            err.append(res["emitter_err_mean_nm"])
        d_est = np.asarray(d_est)
        d_true = np.asarray(d_true)
        err = np.asarray(err)
        rows.append({
            "separation_nm": ds.separation_nm,
            "n_pairs": npair,
            "bin_frames": int(bin_frames), "m_sub": int(m_sub),
            "dt_eff": ds.dt * int(bin_frames),
            "d_est_mean_nm": float(np.mean(d_est)),
            "d_est_std_nm": float(np.std(d_est)),
            "d_bias_nm": float(np.mean(d_est - d_true)),
            "d_rmse_nm": float(np.sqrt(np.mean((d_est - d_true) ** 2))),
            "emitter_err_mean_nm": float(np.mean(err)),
        })
    return rows
