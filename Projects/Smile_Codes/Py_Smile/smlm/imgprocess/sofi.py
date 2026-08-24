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
from .core import (cumulant_from_deltas, deltas_from_lags,
                   entromap, xc_cumulant)
from .interp import fourier_upsample


class Sofit(Transform):
    """SOFI 估计器基类 (Transform 子类)。

    子类只需实现 `_estimate(data:(T,H,W)) -> (H,W)`。输出包成 1 帧 ImageStack,
    并把估计器参数写进 metadata['sofi'] 供追溯。
    """
    method = 'base'

    def _estimate(self, data: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def info(self) -> dict:
        return {'method': self.method}

    def __call__(self, stack):
        data = np.asarray(stack.data, dtype=np.float32)
        img = np.asarray(self._estimate(data), dtype=np.float32)  # (H, W)
        meta = dict(stack.metadata)
        meta['sofi'] = {**self.info()}
        ps = stack.pixel_size
        # 若估计器本身做了上采样, 在这里同步缩小 pixel_size (默认不上采样)
        up = int(meta['sofi'].get('upsample', 1))
        new_ps = (ps / up) if (ps is not None and up > 1) else ps
        return replace(stack, data=img[None], pixel_size=new_ps, metadata=meta)


class CumulantSofit(Sofit):
    """传统 n 阶累积量 SOFI (逐像素时间自累积量)。

    order  : 累积量阶数 m (2-4);
    lags   : 长度为 m 的滞后元组; None -> (0,...,m-1);
    subseq : 分块长度; 块内中心化, 块间累积量平均 (抗非平稳/漂白);
    upsample/taper: 若 >1, 先在每个分块上做傅里叶上采样再算累积量
            (兼容旧 AC 估计器行为; 也推荐改用 Pipeline 里显式的 FourierUpsample)。
    """
    method = 'cumulant'

    def __init__(self, order=2, lags=None, subseq=500,
                 upsample=1, taper=False):
        self.order = int(order)
        self.lags = lags
        self.subseq = int(subseq)
        self.upsample = int(upsample)
        self.taper = bool(taper)

    def info(self):
        return {**super().info(), 'order': self.order, 'lags': self.lags,
                'subseq': self.subseq, 'upsample': self.upsample}

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
            c = cumulant_from_deltas(deltas_from_lags(sub, lags))
            acc = c if acc is None else acc + c
            n += 1
        return acc / n

    def __repr__(self):
        return (f"CumulantSofit(order={self.order}, lags={self.lags}, "
                f"subseq={self.subseq}, upsample={self.upsample})")


class RiedSofit(Sofit):
    """RIED 0-lag 二阶累积量 (+ 可选熵加权)。

    subfactor : 交叉累积量的偏移均值减法因子 (RIED 的 subfactor, 0-1);
    bin       : 熵图直方图分箱数 (RIED 的 bin, 默认 50);
    entropy   : True=熵图×累积量 (完整 RIED 核心); False=仅 0-lag 二阶累积量;
    minv/maxv : 熵图直方图归一化范围 (默认 minv=480, maxv=数据最大值, 对齐 Matlab)。

    输出是 (2H,2W) 网格 (entromap/XCcumulant 各自带来 2× 虚拟像素)。
    """
    method = 'ried'

    def __init__(self, subfactor=0.5, bin=50, entropy=True, minv=480.0):
        self.subfactor = float(subfactor)
        self.bin = int(bin)
        self.entropy = bool(entropy)
        self.minv = float(minv)

    def info(self):
        return {**super().info(), 'subfactor': self.subfactor, 'bin': self.bin,
                'entropy': self.entropy, 'minv': self.minv,
                'upsample': 2}   # 输出为 2× 虚拟像素网格 -> 有效像素减半

    def _estimate(self, data):
        maxv = float(data.max())
        # RIED 的 minv 是经验下限(默认 480); 若数据整体低于它, 退化为数据最小值,
        # 否则直方图退化为单 bin 且 out.max()==0 会产生 NaN。
        minv = self.minv if self.minv < maxv else float(data.min())
        ecum = xc_cumulant(data, 1, self.subfactor)
        if self.entropy:
            entrom = entromap(data, self.bin, maxv, minv)
            out = entrom * ecum
        else:
            out = ecum
        m = out.max()
        return out / m if m > 0 else out

    def __repr__(self):
        return (f"RiedSofit(subfactor={self.subfactor}, bin={self.bin}, "
                f"entropy={self.entropy})")
