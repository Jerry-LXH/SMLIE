# smlm/imgprocess/core.py
"""SOFI / RIED 的核心纯数学 (与具体 Transform 无关, 便于复用与单测)。

两大部分:
  A. 时序累积量 (矩-累积量关系)         -> SOFI 传统 n 阶累积量
  B. 熵加权相关 (熵图 × 交叉累积量)     -> RIED 的 0-lag 二阶累积量 + 可选熵加权

输入输出约定: (T, H, W) float32 数组。
"""
from __future__ import annotations
from math import factorial
import numpy as np

# ======================================================================
# A. 时序累积量 (SOFI 核心)
# ======================================================================

def set_partitions(n):
    """{0..n-1} 的所有集合划分, 每块为非空元组, 元素保持升序。

    递归构造: 把最大元素 n-1 依次并入已有块, 或自成一块。
    """
    if n == 0:
        yield ()
        return
    for part in set_partitions(n - 1):
        for i in range(len(part)):
            yield part[:i] + (part[i] + (n - 1,),) + part[i + 1:]
        yield part + ((n - 1,),)


def cumulant_from_deltas(deltas):
    """由 m 路"已减去各自时间均值"的信号计算 m 阶联合累积量。

    Leonov-Shiryaev / 矩-累积量关系:
        kappa(X_1..X_m) = sum_{pi} (-1)^{|pi|-1} (|pi|-1)! * prod_{B in pi} E[prod_{i in B} X_i]
    因为已经中心化, 含单元素块的划分整体为 0, 直接跳过。

    m=2: <d0 d1>
    m=3: <d0 d1 d2>
    m=4: <d0d1d2d3> - <d0d1><d2d3> - <d0d2><d1d3> - <d0d3><d1d2>
    """
    deltas = list(deltas)
    m = len(deltas)
    total = None
    for part in set_partitions(m):
        if any(len(b) == 1 for b in part):
            continue
        k = len(part)
        coef = ((-1) ** (k - 1)) * factorial(k - 1)
        term = None
        for block in part:
            prod = deltas[block[0]]
            for i in block[1:]:
                prod = prod * deltas[i]
            mom = prod.mean(axis=0)
            term = mom if term is None else term * mom
        total = coef * term if total is None else total + coef * term
    return total


def deltas_from_lags(sub, lags):
    """按滞后 tau 切片并逐路中心化, 得到等长的 (T-maxlag, H, W) 片段。

    lags 必须互不相同 (尤其不能全 0): 若同一时刻自乘, 散粒噪声会在每个
    像素上产生一个正的 shot-noise 基底。
    """
    maxlag = max(lags)
    L = sub.shape[0] - maxlag
    out = []
    for tau in lags:
        d = sub[tau:tau + L]
        out.append(d - d.mean(axis=0, keepdims=True))
    return out


# ======================================================================
# B. 熵加权相关 (RIED 核心)
# ======================================================================

def calprob(trace, bin, maxv, minv):
    """时间轨迹直方图 -> 概率 (对应 RIED Utils/calprob.m)。

    trace: 一维时间序列; bin: 直方图分箱数; maxv/minv: 归一化上下限。
    """
    trace = np.asarray(trace)
    idx = np.round((trace - minv) / (maxv - minv) * bin).astype(np.int64)
    idx[idx < 0] = 0
    idx[idx >= bin] = bin - 1
    out = np.zeros(bin, dtype=np.float64)
    np.add.at(out, idx, 1)
    return out / len(trace)


def x_entropy(d1, d2, bin, maxv, minv):
    """两条时间轨迹的对称交叉熵项 (对应 RIED 的 xEtr.m)。

        xEtr = -(1/2) * Σ_b [ P2[b]*log(P1[b]) + P1[b]*log(P2[b]) ]

    当 d1==d2 时退化为香农熵 -Σ P log P。
    """
    P1 = calprob(d1, bin, maxv, minv)
    P2 = calprob(d2, bin, maxv, minv)
    lp1 = np.where(P1 > 0, np.log(P1), 0.0)
    lp2 = np.where(P2 > 0, np.log(P2), 0.0)
    return -0.5 * float(np.sum(P2 * lp1 + P1 * lp2))


def _hist_field(data, bin, maxv, minv, chunk=200):
    """逐像素时间直方图 -> (H, W, bin) 的 P 与 logP (向量化, 加速 entromap)。

    data: (T, H, W) —— 时间在轴 0。按帧分块累加, 控制峰值内存。
    """
    T, H, W = data.shape
    counts = np.zeros((H * W, bin), dtype=np.float64)
    pix = np.arange(H * W, dtype=np.int32).reshape(H, W)[None, ...]        # (1,H,W)
    for s in range(0, T, chunk):
        sub = data[s:s + chunk]
        idx = np.round((sub - minv) / (maxv - minv) * bin).astype(np.int32)  # (c,H,W)
        np.clip(idx, 0, bin - 1, out=idx)
        flat = pix * bin + idx                                             # (c,H,W), int32
        counts += np.bincount(flat.ravel(), minlength=H * W * bin).reshape(H * W, bin)
    P = (counts / T).reshape(H, W, bin)
    logP = np.zeros_like(P)
    np.log(P, out=logP, where=P > 0)          # P==0 处保持 0, 避免 log(0) 警告
    return P, logP


def _xetr(Pa, Pb, logPa, logPb):
    """两个 (...,bin) 直方图场之间的 xEtr 场: -0.5 Σ (Pb*logPa + Pa*logPb)。"""
    return -0.5 * np.sum(Pb * logPa + Pa * logPb, axis=-1)


def entromap(data, bin, maxv, minv):
    """2×2 邻域交叉熵图 (对应 RIED 的 entromap.m): (T,H,W) -> (2H,2W)。

    每个 2×2 块(i,j)在输出 2×2 位置生成:
      (0,0): xEtr(p(i,j),   p(i,j))       —— 自熵
      (0,1): xEtr(p(i,j),   p(i,j+1))     —— 水平对
      (1,0): xEtr(p(i,j),   p(i+1,j))     —— 垂直对
      (1,1): [xEtr(p(i,j),p(i+1,j+1)) + xEtr(p(i+1,j),p(i,j+1))]/2 —— 对角线对
    """
    data = np.asarray(data)
    H, W = data.shape[1], data.shape[2]
    P, logP = _hist_field(data, bin, maxv, minv)
    out = np.zeros((2 * H, 2 * W), dtype=np.float32)
    # 只对 i=0..H-2, j=0..W-2 的块计算 (Matlab 循环到 x-1, y-1) -> 末行/列保持 0
    out[0:2 * H - 2:2, 0:2 * W - 2:2] = _xetr(P[:-1, :-1], P[:-1, :-1],
                                              logP[:-1, :-1], logP[:-1, :-1])
    out[0:2 * H - 2:2, 1:2 * W - 2:2] = _xetr(P[:-1, :-1], P[:-1, 1:],
                                              logP[:-1, :-1], logP[:-1, 1:])
    out[1:2 * H - 2:2, 0:2 * W - 2:2] = _xetr(P[:-1, :-1], P[1:, :-1],
                                              logP[:-1, :-1], logP[1:, :-1])
    d1 = _xetr(P[:-1, :-1], P[1:, 1:], logP[:-1, :-1], logP[1:, 1:])
    d2 = _xetr(P[:-1, 1:], P[1:, :-1], logP[:-1, 1:], logP[1:, :-1])
    out[1:2 * H - 2:2, 1:2 * W - 2:2] = 0.5 * (d1 + d2)
    return out


def xc_cumulant(data, order=1, offset=0.5):
    """RIED 版交叉累积量 (对应 RIED 的 XCcumulant.m): (T,H,W) -> (2H,2W)。

        cdata = |data - offset * mean_t(data)| ^ order
      (0,0): <cdata[t] * cdata[t+1]>_t                  —— 时间对
      (0,1): <cdata[i,j] * cdata[i,j+1]>                —— 水平对
      (1,0): <cdata[i,j] * cdata[i+1,j]>                —— 垂直对
      (1,1): (<cdata[i,j]cdata[i+1,j+1]> + <cdata[i+1,j]cdata[i,j+1]>)/2

    注意: 这是"偏移均值差 × 相邻像素乘积"的 0-lag 二阶统计量,
    **没有** XC-SOFI 那套平场校正 (RIED 原版如此)。
    """
    data = np.asarray(data, dtype=np.float32)          # (T,H,W)
    H, W = data.shape[1], data.shape[2]
    cdata = np.abs(data - offset * data.mean(axis=0, keepdims=True)) ** order
    out = np.zeros((2 * H, 2 * W), dtype=np.float32)
    out[0::2, 0::2] = np.mean(cdata[:-1, :, :] * cdata[1:, :, :], axis=0)       # 时间对
    out[0::2, 1:2 * W - 2:2] = np.mean(cdata[:, :, :-1] * cdata[:, :, 1:],
                                       axis=0)                                   # 水平对 (所有偶数行)
    out[1:2 * H - 2:2, 0::2] = np.mean(cdata[:, :-1, :] * cdata[:, 1:, :],
                                       axis=0)                                   # 垂直对
    out[1:2 * H - 2:2, 1:2 * W - 2:2] = 0.5 * (
        np.mean(cdata[:, :-1, :-1] * cdata[:, 1:, 1:], axis=0) +
        np.mean(cdata[:, :-1, 1:] * cdata[:, 1:, :-1], axis=0))                  # 对角线对
    return out


def entrocor(data, bin=50, subfactor=0.5, maxv=None, minv=480.0):
    """RIED 熵加权相关 (对应 entrocor.m): EC = 熵图 × 交叉累积量, 归一化。

    data: (T,H,W) 预处理/预反卷积后的栈。返回 (2H,2W)。
    """
    data = np.asarray(data, dtype=np.float32)
    if maxv is None:
        maxv = float(data.max())
    entrom = entromap(data, bin, maxv, minv)
    ecum = xc_cumulant(data, 1, subfactor)
    out = entrom * ecum
    return (out / out.max()).astype(np.float32)
