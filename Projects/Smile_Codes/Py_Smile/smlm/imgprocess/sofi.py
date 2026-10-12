# smlm/imgprocess/sofi.py
"""SOFI / RIED 类超分辨核心, 重构为 `Transform` 子类。

所有估计器都是 ImageStack 进、单帧 ImageStack 出 (SOFI 输出是二维超分辨图,
包成 1 帧 ImageStack 以便在 Pipeline 里继续接 FourierUpsample / SparseDeconv
等后处理)。算法细节在 `core.py`。

实现:
  - `CumulantSofit` : 传统 n 阶累积量 SOFI (分块中心化 + 块间平均)。
  - `RiedSofit`     : RIED 0-lag 二阶累积量 (XCcumulant), 可选熵加权 (entromap)。

注: 上采样不再内嵌在估计器里 —— 用 `FourierUpsample` 作为 Pipeline 的显式
阶段, 便于"逐步增加 pipeline 长度"的实验。
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from .base import Transform
from .core import (cumulant_from_deltas, entromap, xc_cumulant)
from .interp import fourier_upsample


class Sofit(Transform):
    """SOFI / RIED 估计器基类 (Transform 子类, 鸭子类型可直接喂 ndarray)。

    子类只需实现 `_estimate(data:(T,H,W)) -> (H,W)`。ImageStack 进时输出包成
    1 帧 ImageStack 并把估计器参数写进 metadata['sofi'] 供追溯 (自动按
    info()['upsample'] 缩小 pixel_size)。
    """
    method = 'base'

    def _estimate(self, data: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def info(self) -> dict:
        return {'method': self.method}

    def _apply(self, data):
        return np.asarray(self._estimate(np.asarray(data, dtype=np.float32)),
                          dtype=np.float32)          # (H, W)

    def _wrap(self, stack, img):
        meta = dict(stack.metadata)
        meta['sofi'] = {**self.info()}
        ps = stack.pixel_size
        # 若估计器本身做了上采样, 在这里同步缩小 pixel_size (默认不上采样)
        up = int(meta['sofi'].get('upsample', 1))
        new_ps = (ps / up) if (ps is not None and up > 1) else ps
        return replace(stack, data=img[None], pixel_size=new_ps, metadata=meta)


class CumulantSofit(Sofit):
    """传统 n 阶累积量 SOFI / 也兼容 SACD 式时间累积量 (逐像素时间自累积量)。

    order  : 累积量阶数 m (2-...);
    lags   : 长度为 m 的滞后元组; None -> (0,...,m-1);
    subseq : 分块长度; 块内中心化, 块间累积量平均 (抗非平稳/漂白);
    upsample/taper: 若 >1, 先在每个分块上做傅里叶上采样再算累积量;
    subfactor: None(默认) = 经典 SOFI —— 各路按各自时间均值完全中心化后再求
             m 阶联合累积量 (真累积量);
             0..1 数值 = SACD 式 —— 各路只减 subfactor·(块时间均值), 不再逐路
             中心化 (可配 take_abs), 直接对 (可能带 DC 的) 信号按同一 m 阶
             “联合累积量”公式做滞后乘积 —— 复现 SACDm 的 cumulant() (order
             2/3/4 与其硬编码公式一致)。
    take_abs: False(默认)。True 时对 (data-subfactor·μ) 取绝对值 —— 即 SACD
             的 `|I - subfactor·mean(I)|` 预处理 (需 subfactor 非 None)。
    """
    method = 'cumulant'

    def __init__(self, order=2, lags=None, subseq=500,
                 upsample=1, taper=False, subfactor=None, take_abs=False):
        self.order = int(order)
        self.lags = lags
        self.subseq = int(subseq)
        self.upsample = int(upsample)
        self.taper = bool(taper)
        self.subfactor = None if subfactor is None else float(subfactor)
        self.take_abs = bool(take_abs)
        if self.take_abs and self.subfactor is None:
            raise ValueError("take_abs=True 需同时给 subfactor(SACD式: 减 sub·μ 后取绝对值);"
                             "经典中心化路径默认不加 abs。")

    def info(self):
        return {**super().info(), 'order': self.order, 'lags': self.lags,
                'subseq': self.subseq, 'upsample': self.upsample,
                'subfactor': self.subfactor, 'take_abs': self.take_abs}

    def _estimate(self, data):
        order = self.order
        lags = self.lags if self.lags is not None else tuple(range(order))
        assert len(lags) == order, "lags 的个数必须等于阶数"
        if data.shape[0] <= max(lags):
            raise ValueError(
                f"帧数 {data.shape[0]} 不足: 需要 > max(lags) = {max(lags)}")

        T = data.shape[0]
        subseq = min(self.subseq, T)
        acc, n = None, 0
        for s in range(0, T - subseq + 1, subseq):
            sub = data[s:s + subseq]
            if self.upsample > 1:
                sub = fourier_upsample(sub, self.upsample, taper=self.taper)
            maxlag = max(lags)
            L = sub.shape[0] - maxlag
            X = [sub[tau:tau + L] for tau in lags]
            if self.subfactor is None:
                # 经典 SOFI: 各路独立中心化 (真 m 阶累积量)
                D = [x - x.mean(axis=0, keepdims=True) for x in X]
            else:
                # SACD 式: 各路只减 subfactor·块均值 (不再逐路中心化)
                mu = sub.mean(axis=0, keepdims=True)
                D = [x - self.subfactor * mu for x in X]
                if self.take_abs:
                    D = [np.abs(d) for d in D]
            c = cumulant_from_deltas(D)
            acc = c if acc is None else acc + c
            n += 1
        return acc / n

    def __repr__(self):
        return (f"CumulantSofit(order={self.order}, lags={self.lags}, "
                f"subseq={self.subseq}, upsample={self.upsample}, "
                f"subfactor={self.subfactor}, take_abs={self.take_abs})")


class RiedSofit(Sofit):
    """RIED 0-lag 二阶累积量 (+ 可选熵加权)。

    subfactor : 交叉累积量的偏移均值减法因子 (RIED 的 subfactor, 0-1);
    bin       : 熵图直方图分箱数 (RIED 的 bin, 默认 50);
    entropy   : True=熵图×累积量 (完整 RIED 核心); False=仅 0-lag 二阶累积量;
    minv/maxv : 熵图直方图归一化范围。
    use_abs   : True(默认)=相关前对 (data-subfactor·μ) 取绝对值 (RIED 原版行为);
                False=保留符号 —— 配合 subfactor=1、entropy=False 时即退化为
                去均值信号的滞后1协方差, 与 CumulantSofit(order=2, lags=(0,1)) 基本一致。

    **minv 语义**：
      - `minv=None`(默认): 自动取 `data.min()` —— 直方图覆盖全动态范围。
        这是**正确**的做法: RIED 原版写死 480(按 ADU 尺度), 在光子尺度数据上
        会让大多数像素落入第 0 个直方图箱 → 熵图退化 → 熵加权过度稀疏、丢结构;
      - `minv=480.0`: 显式传 480 以**复现 matlab/RIED 原版行为**(ADU 数据适用)。

    输出是 (2H,2W) 网格 (entromap/XCcumulant 各自带来 2× 虚拟像素)。
    """
    method = 'ried'

    def __init__(self, subfactor=0.5, bin=50, entropy=True, minv=None, use_abs=True):
        self.subfactor = float(subfactor)
        self.bin = int(bin)
        self.entropy = bool(entropy)
        self.minv = None if minv is None else float(minv)
        self.use_abs = bool(use_abs)

    def info(self):
        return {**super().info(), 'subfactor': self.subfactor, 'bin': self.bin,
                'entropy': self.entropy, 'minv': self.minv, 'use_abs': self.use_abs,
                'upsample': 2}   # 输出为 2× 虚拟像素网格 -> 有效像素减半

    def _estimate(self, data):
        maxv = float(data.max())
        # minv=None → 自动取 data.min()(直方图覆盖全范围, 避免 480 造成退化);
        # 显式 minv 且 < maxv 时用之(如 480 复现 matlab); 否则回退 data.min()。
        if self.minv is None or self.minv >= maxv:
            minv = float(data.min())
        else:
            minv = self.minv
        ecum = xc_cumulant(data, 1, self.subfactor, self.use_abs)
        if self.entropy:
            entrom = entromap(data, self.bin, maxv, minv)
            out = entrom * ecum
        else:
            out = ecum
        m = out.max()
        return out / m if m > 0 else out

    def __repr__(self):
        return (f"RiedSofit(subfactor={self.subfactor}, bin={self.bin}, "
                f"entropy={self.entropy}, use_abs={self.use_abs})")
