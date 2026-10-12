# smlm/localize/loc/fisher.py
"""Fisher 信息: 共享的矩阵机制。

设计取舍: 各模型的**参数导数 Jacobian** 是方案相关的, 留在各自的方案模块里
(gaussian2d / astigmatic / doublehelix); 这里只共享"由 Jacobian 组装 Fisher
矩阵"的噪声机制, 避免重复实现噪声模型 W = 1/(F_EM^2 * mu) + 1/(2 mu^2)。
"""
from __future__ import annotations
import numpy as np


def fisher_from_jacobian(J, mu, F_EM):
    """由模型对各参数的一阶导数 Jacobian 计算 Fisher 信息矩阵。

    参数
    ----
    J    : (P,) + roi.shape, 第 i 行是 d model / d theta_i。
    mu   : 模型预测 (与 J 的像素数一致); 内部加 1e-12 防除零。
    F_EM : EMCCD 过噪因子 (非 EMCCD 为 1)。

    返回
    ----
    (P, P) Fisher 信息矩阵。协方差近似为它的逆。
    """
    J = np.asarray(J, dtype=float)
    mu = np.asarray(mu, dtype=float).ravel() + 1e-12
    W = 1.0 / (F_EM ** 2 * mu) + 1.0 / (2 * mu ** 2)
    J2 = J.reshape(J.shape[0], -1)
    return J2 @ (W[:, None] * J2.T)
