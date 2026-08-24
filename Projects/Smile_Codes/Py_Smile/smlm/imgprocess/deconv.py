# smlm/imgprocess/deconv.py
"""反卷积类 Transform: RL3D / 稀疏反卷积 (split-Bregman) / Wiener(OTF 重加权)。

- `RL3D`           : 加速 Richardson-Lucy (对应 RIED 的 RL3D / RLdeconv)。
- `SparseDeconv`   : split-Bregman TV+L1 稀疏反卷积 (对应 RIED 的 Spasedeconv)。
- `OtfReweight`    : Wiener 型 OTF 重加权 (对应 SOFI 的 otf_reweight, Stein 2015 式 9)。
同时提供纯函数 (`rl_deconv`, `sparse_deconv`, `otf_reweight`, `scalar_psf` ...)
供直接调用或与其它 Transform 组合。
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.fft import fft2, ifft2, fftn, ifftn, fftfreq, fftshift, ifftshift
from scipy.special import j0, j1, struve
from .base import Transform

# ======================================================================
# PSF / OTF 工具
# ======================================================================

def psf2otf(psf, shape):
    """Matlab `psf2otf` 等价: 把 PSF 中心对齐到 (0,0)、补零到 shape 后 FFT。"""
    psf = np.asarray(psf, dtype=np.float32)
    otf = np.zeros(shape, dtype=np.float32)
    otf[tuple(slice(0, p) for p in psf.shape)] = psf
    for axis, axis_size in enumerate(psf.shape):
        otf = np.roll(otf, -int(axis_size / 2), axis=axis)
    return fftn(otf)


def scalar_psf(shape, pixelsize, na, wavelength, z=0.0, nsub=2000):
    """RIED `Generate_PSF` 等价的标量衍射积分 PSF (Debye-Wolf, 可离焦)。

    与 matlab `kernel()`/`Generate_PSF.m` 同一公式:
        u   = 8π z sin²(asin NA/2) / λ
        IP(r) = ∫₀¹ 2 exp(i u p²/2) J0(2π r NA/λ p) dp
        PSF  = |IP|² (再归一化到总和 1)
    其中 p 为归一化光瞳半径, r 为像面径向距离 (与 pixelsize 同单位)。

    - z=0 (焦面): 用 Struve 闭式解 ∫₀^x J0 = xJ0 + (πx/2)(J1·H0 - J0·H1),
      向量化、快 (与数值积分一致到 ~1e-13);
    - z≠0 (离焦): 退化为对 p∈[0,1] 的分块梯形数值积分 (nsub 采样数)。

    ⚠️ 已知差异: matlab `Generate_PSF.m` 里有 `idx = s1 <= 1` (s1 为 nm 单位
    半径), 实测会让 PSF 只剩中心 1 像素 (≈ delta, RL 退化为无操作)。本实现
    **不复刻该掩膜**, 在全网格上积分, 得到真实的衍射 PSF。
    """
    ny, nx = shape
    y = (np.arange(ny) - ny // 2) * pixelsize
    x = (np.arange(nx) - nx // 2) * pixelsize
    R = np.hypot(*np.meshgrid(y, x, indexing="ij"))
    v = 2 * np.pi * na * R / wavelength
    if z == 0:
        xv = np.asarray(v, dtype=float)
        ip = np.empty_like(xv)
        nz = xv > 0
        ip[~nz] = 2.0
        xz = xv[nz]
        intg = xz * j0(xz) + (np.pi * xz / 2) * (j1(xz) * struve(0, xz) - j0(xz) * struve(1, xz))
        ip[nz] = 2 * intg / xz
    else:
        sin2 = (1 - (1 - na ** 2)) / 2
        u = 8 * np.pi * z * sin2 / wavelength
        p = np.linspace(0, 1, nsub)
        phase = np.exp(1j * u * p ** 2 / 2)
        ip = np.zeros(R.shape, dtype=complex)
        for c0 in range(0, nx, 16):              # 分块控制峰值内存
            c1 = min(c0 + 16, nx)
            integ = 2 * phase[None, None, :] * j0(v[:, c0:c1, None] * p[None, None, :])
            ip[:, c0:c1] = np.trapezoid(integ, x=p, axis=-1)
    psf = np.abs(ip) ** 2
    s = psf.sum()
    return (psf / s).astype(np.float32) if s > 0 else psf.astype(np.float32)


def otf_components(shape, order, pixelsize, na, wavelength):
    """OTF 重加权所需的频域量 (供 otf_reweight 与可视化共用)。

    返回 (otf, ideal, k, kc):
      - otf  : 真实的 m 阶 SOFI OTF = |FT{PSF^m}|, 归一化到 [0,1];
      - ideal: 理想共焦 OTF (均匀圆盘自卷积), 支撑到 m*kc;
      - k    : 空间频率模长 (1/长度, 与 pixelsize 同单位);
      - kc   : 宽场截止频率 = 2*na/wavelength。
    """
    ny, nx = shape
    ky = fftfreq(ny, d=pixelsize)
    kx = fftfreq(nx, d=pixelsize)
    k = np.hypot(*np.meshgrid(ky, kx, indexing="ij"))
    kc = 2 * na / wavelength                         # 宽场截止频率

    psf = scalar_psf((ny, nx), pixelsize, na, wavelength)   # 唯一 PSF 生成器
    otf = np.abs(fft2(np.fft.ifftshift(psf ** order)))
    otf /= otf.max()

    u = np.clip(k / (order * kc), 0, 1)
    ideal = (2 / np.pi) * (np.arccos(u) - u * np.sqrt(np.maximum(1 - u ** 2, 0)))
    ideal[k >= order * kc] = 0.0
    return otf, ideal, k, kc


def otf_reweight(img, order, pixelsize, na, wavelength, eps=0.1):
    """m 阶 SOFI 的 OTF 重加权 —— 正则化逆滤波 (Wiener 型反卷积, Stein 2015 式 9)。

        I_ideal(k) = U_ISM(2k/m) / (eps + OTF_wf^{(*m)}(k)) * I_sofi(k)

    它是"真信息恢复": 只在 |k| < m*kc 内放大已被 U^m 衰减的高频, 不生成
    支撑域外的新频率; 但受 SNR 限制, eps 越小越锐利也越脏。
    """
    img = np.asarray(img, dtype=np.float32)
    otf, ideal, _, _ = otf_components(img.shape, order, pixelsize, na, wavelength)
    out = np.real(ifft2(fft2(img) * (ideal / (eps + otf))))
    return np.ascontiguousarray(out, dtype=np.float32)


# ======================================================================
# Richardson-Lucy (加速版, 对应 RIED 的 RL3D / RLdeconv)
# ======================================================================

def _rl_core(yk, data, otf, xk, vk, it):
    """单步加速 RL 更新 (yk/xk/vk 对应 Matlab deblur_coreRL)。"""
    rel = fftn(data / np.maximum(np.real(ifftn(otf * fftn(yk))), 1e-5))
    xk_update = xk
    xk = np.maximum(yk * np.real(ifftn(np.conj(otf) * rel)), 1e-5)
    vk_update = vk
    vk = np.maximum(xk - yk, 1e-5)
    if it == 1:
        alpha = 0.0
    else:
        denom = np.sum(vk_update * vk_update) + np.finfo(float).eps
        alpha = np.sum(vk_update * vk) / denom
        alpha = max(min(alpha, 1.0), 0.0)
    yk = np.maximum(np.real(xk + alpha * (xk - xk_update)), 1e-5)
    return yk, xk, vk


def rl_deconv(data, psf, iterations):
    """加速 Richardson-Lucy 反卷积 (2D)。

    data: (H, W); psf: 2D 核 (内部会归一化); iterations: 迭代次数。
    """
    kernel = np.asarray(psf, dtype=np.float32)
    kernel = kernel / kernel.sum()
    otf = psf2otf(kernel, np.asarray(data).shape)
    yk = np.asarray(data, dtype=np.float32).copy()
    xk = np.zeros_like(yk)
    vk = np.zeros_like(yk)
    data = np.asarray(data, dtype=np.float32)
    for it in range(1, int(iterations) + 1):
        yk, xk, vk = _rl_core(yk, data, otf, xk, vk, it)
    return np.real(yk)


class RL3D(Transform):
    """Richardson-Lucy 反卷积 Transform (逐帧应用)。

    psf: 2D PSF 核 (RIED 里是 `kernel(pixel, wavelength, NA)` 生成的 PSF);
    iterations: 迭代次数 (RIED 的 iter1 / iter2)。
    """
    def __init__(self, psf, iterations=4):
        self.psf = np.asarray(psf, dtype=np.float32)
        self.iterations = int(iterations)

    def __call__(self, stack):
        data = np.asarray(stack.data, dtype=np.float32)
        out = np.empty_like(data)
        for i in range(data.shape[0]):
            out[i] = np.abs(rl_deconv(data[i], self.psf, self.iterations))
        return replace(stack, data=out)

    def __repr__(self):
        return f"RL3D(psf={self.psf.shape}, iterations={self.iterations})"


# ======================================================================
# 有限差分算子 (对应 RIED Utils/forward_diff.m, back_diff.m)
# ======================================================================

def forward_diff(data, step=1, dim=0):
    """前向差分 (data[k+1]-data[k], 远边界视为 0)。dim ∈ {0,1,2}。

    对应 Matlab `forward_diff(data,step,dim)` 的输出 (含其末尾取反后的结果)。
    """
    d = np.diff(data, axis=dim, append=0)   # append=0 -> 最后一项 = 0 - data[-1]
    return (d / step).astype(data.dtype)


def back_diff(data, step=1, dim=0):
    """后向差分 (data[k]-data[k-1], 近边界视为 0)。dim ∈ {0,1,2}。

    对应 Matlab `back_diff(data,step,dim)`。
    """
    d = np.diff(data, axis=dim, prepend=0)  # prepend=0 -> 首项 = data[0] - 0
    return (d / step).astype(data.dtype)


def _shrink(v, t):
    """软阈值 (shrinkage): sign(v) * max(|v|-t, 0)。对应 split-Bregman 的 d 更新。"""
    return np.sign(v) * np.maximum(np.abs(v) - t, 0.0)


# ======================================================================
# 稀疏反卷积 (split-Bregman TV+L1, 对应 RIED 的 Spasedeconv)
# ======================================================================

def sparse_deconv(input_img, mu, sigma_t, l10, iter_bregman,
                  psf=None, iter=5, scale=1.0):
    """split-Bregman 稀疏反卷积 (对应 Matlab Spasedeconv)。

    参数:
      input_img   : (H, W) 待反卷积图像 (RIED 里是熵加权相关后的 EC_FI);
      mu          : fidelity 权重 (RIED 的 fidelity, 默认 50);
      sigma_t     : z 方向权重 (RIED 传 0, 2D 时无影响);
      l10         : L1 稀疏权重 (RIED 的 sparsity);
      iter_bregman: Bregman 外迭代次数 (RIED 传 100);
      psf         : 末段 RL 反卷积核; iter: 末段 RL 迭代次数 (RIED 的 iter2);
      scale       : PSF 幂次缩放 (RIED 的 scale, 默认 1)。

    返回归一化到 [0,1] 的 (H, W) float32 图像。
    """
    lamda = 1.0
    siranu = float(mu)
    zbei = float(sigma_t)
    input_img = np.asarray(input_img, dtype=np.float32)
    sx, sy = input_img.shape[:2]
    y = np.zeros((sx, sy, 1), dtype=np.float32)
    y[:, :, 0] = input_img
    sizex = y.shape
    y_flag = sizex[2]
    l1 = float(l10)

    # --- 频域分母: Frefft = Σ |FT{差分核}|², divide = mu/lamda + Frefft ---
    kx  = np.zeros((3, 1, 1)); kx[:, 0, 0] = [1, -2, 1]      # 沿 axis0 (x)
    ky  = np.zeros((1, 3, 1)); ky[0, :, 0] = [1, -2, 1]      # 沿 axis1 (y)
    kz  = np.zeros((1, 1, 3)); kz[0, 0, :] = [1, -2, 1]      # 沿 axis2 (z)
    kxy = np.zeros((2, 2, 1)); kxy[:, :, 0] = [[1, -1], [-1, 1]]
    kxz = np.zeros((1, 2, 2)); kxz[0, :, 0] = [1, -1]; kxz[0, :, 1] = [-1, 1]
    kyz = np.zeros((2, 1, 2)); kyz[:, 0, 0] = [1, -1]; kyz[:, 0, 1] = [-1, 1]

    def f(k): return np.fft.fftn(k, s=sizex)

    Frefft = f(ky) * np.conj(f(ky))            # dim1 (y) 二阶差分
    Frefft = Frefft + f(kx) * np.conj(f(kx))   # dim0 (x) 二阶差分
    Frefft = Frefft + (zbei ** 2) * f(kz) * np.conj(f(kz))
    Frefft = Frefft + 2 * f(kxy) * np.conj(f(kxy))
    Frefft = Frefft + 2 * zbei * f(kxz) * np.conj(f(kxz))
    Frefft = Frefft + 2 * zbei * f(kyz) * np.conj(f(kyz))
    divide = (siranu / lamda) + Frefft

    b = [np.zeros(sizex, np.float32) for _ in range(8)]
    x = np.zeros(sizex, np.float32)
    frac = (siranu / lamda) * y
    frac[frac < 0] = 0

    thr = 1.0 / lamda
    for ii in range(1, int(iter_bregman) + 1):
        frac_f = np.fft.fftn(frac)
        if ii > 1:
            x = np.real(np.fft.ifftn(frac_f / divide))
        else:
            x = np.real(np.fft.ifftn(frac_f / (siranu / lamda)))
        frac = (siranu / lamda) * y

        # 1) x 方向二阶差分
        u = back_diff(forward_diff(x, 1, 0), 1, 0)
        d = _shrink(u + b[0], thr); b[0] = b[0] + (u - d)
        frac = frac + back_diff(forward_diff(d - b[0], 1, 0), 1, 0)
        # 2) y 方向二阶差分
        u = back_diff(forward_diff(x, 1, 1), 1, 1)
        d = _shrink(u + b[1], thr); b[1] = b[1] + (u - d)
        frac = frac + back_diff(forward_diff(d - b[1], 1, 1), 1, 1)
        # 3) z 方向二阶差分
        u = back_diff(forward_diff(x, 1, 2), 1, 2)
        d = _shrink(u + b[2], thr); b[2] = b[2] + (u - d)
        frac = frac + (zbei ** 2) * back_diff(forward_diff(d - b[2], 1, 2), 1, 2)
        # 4) xy 交叉
        u = forward_diff(forward_diff(x, 1, 0), 1, 1)
        d = _shrink(u + b[3], thr); b[3] = b[3] + (u - d)
        frac = frac + 2 * back_diff(back_diff(d - b[3], 1, 1), 1, 0)
        # 5) xz 交叉
        u = forward_diff(forward_diff(x, 1, 0), 1, 2)
        d = _shrink(u + b[4], thr); b[4] = b[4] + (u - d)
        frac = frac + 2 * zbei * back_diff(back_diff(d - b[4], 1, 2), 1, 0)
        # 6) yz 交叉
        u = forward_diff(forward_diff(x, 1, 1), 1, 2)
        d = _shrink(u + b[5], thr); b[5] = b[5] + (u - d)
        frac = frac + 2 * zbei * back_diff(back_diff(d - b[5], 1, 2), 1, 1)
        # 7) L1 稀疏 (对 x 本身)
        u = x
        d = _shrink(u + b[6], thr); b[6] = b[6] + (u - d)
        frac = frac + l1 * (d - b[6])

        x[x < 0] = 0

    x[x < 0] = 0
    x = x[:, :, :y_flag]
    cr = x.astype(np.float32)
    cr = cr / cr.max()

    if psf is not None and int(iter) > 0:
        out2d = np.abs(rl_deconv(cr[:, :, 0],
                                 np.asarray(psf, dtype=np.float32) ** scale,
                                 int(iter)))
    else:
        out2d = np.abs(cr[:, :, 0])
    out2d = out2d / out2d.max()
    return np.ascontiguousarray(out2d, dtype=np.float32)


class SparseDeconv(Transform):
    """稀疏反卷积 Transform (对应 RIED 的 Spasedeconv)。

    psf: 后反卷积 PSF; fidelity=mu; sparsity=l10;
    iter_bregman: Bregman 迭代次数 (RIED 固定 100);
    iterations : 末段 RL 迭代次数 (RIED 的 iter2);
    scale      : PSF 幂次缩放 (RIED 的 scale)。
    """
    def __init__(self, psf, fidelity=50, sparsity=0.5, iter_bregman=100,
                 iterations=5, scale=1.0, sigma_t=0.0):
        self.psf = np.asarray(psf, dtype=np.float32)
        self.mu = float(fidelity)
        self.l10 = float(sparsity)
        self.iter_bregman = int(iter_bregman)
        self.iterations = int(iterations)
        self.scale = scale
        self.sigma_t = float(sigma_t)

    def __call__(self, stack):
        data = np.asarray(stack.data, dtype=np.float32)
        out = np.empty_like(data)
        for i in range(data.shape[0]):
            out[i] = sparse_deconv(data[i], self.mu, self.sigma_t, self.l10,
                                   self.iter_bregman, self.psf,
                                   self.iterations, self.scale)
        return replace(stack, data=out)

    def __repr__(self):
        return (f"SparseDeconv(fidelity={self.mu}, sparsity={self.l10}, "
                f"iter_bregman={self.iter_bregman}, rl_iter={self.iterations})")


# ======================================================================
# Wiener 型 OTF 重加权 Transform
# ======================================================================

class OtfReweight(Transform):
    """OTF 重加权 (Wiener 型反卷积)。需要 ImageStack.pixel_size (nm/px)。

    对 m 阶 SOFI 输出: 放大 |k|<m*kc 内被 U^m 衰减的高频 (真信息恢复)。
    """
    def __init__(self, order, na, wavelength, eps=0.1):
        self.order = int(order)
        self.na = float(na)
        self.wavelength = float(wavelength)
        self.eps = float(eps)

    def __call__(self, stack):
        ps = stack.pixel_size
        if ps is None:
            raise ValueError("OtfReweight 需要 ImageStack.pixel_size (nm/px)")
        data = np.asarray(stack.data, dtype=np.float32)
        out = np.empty_like(data)
        for i in range(data.shape[0]):
            out[i] = otf_reweight(data[i], self.order, ps, self.na,
                                  self.wavelength, self.eps)
        return replace(stack, data=out)

    def __repr__(self):
        return f"OtfReweight(order={self.order}, na={self.na}, wl={self.wavelength})"
