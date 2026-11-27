# =====================================================================
#  过滤器: 基类 + 可组合管线
# =====================================================================
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Iterable
import numpy as np

from ..datatypes import Emitter


class EmitterFilter(ABC):
    @abstractmethod
    def keep(self, emitters) -> np.ndarray: ...          # bool mask

    def __call__(self, emitters):
        m = self.keep(emitters)
        kept = [e for e, k in zip(emitters, m) if k]
        print(f'{type(self).__name__}: 保留 {len(kept)} / {len(emitters)}')
        return kept


class FilterPipeline(EmitterFilter):
    def __init__(self, *filters):
        self.filters = filters
    def keep(self, emitters):                            # 满足接口
        return np.ones(len(emitters), bool)
    def __call__(self, emitters):
        for filt in self.filters:
            emitters = filt(emitters)
        return emitters


class MinDurationFilter(EmitterFilter):
    """短寿命过滤。mode='total' 计总点数, 'consecutive' 计最长连续段。"""
    def __init__(self, min_frames, mode='consecutive'):
        self.min_frames = min_frames; self.mode = mode
    def _longest_run(self, fr):
        if fr.size == 0: return 0
        breaks = np.where(np.diff(fr) > 1)[0]
        seg = np.diff(np.concatenate([[-1], breaks, [fr.size - 1]]))
        return int(seg.max())
    def keep(self, emitters):
        if self.mode == 'consecutive':
            v = np.array([self._longest_run(e.frame) for e in emitters])
        else:
            v = np.array([e.n_locs for e in emitters])
        return v >= self.min_frames


class FirstFrameFilter(EmitterFilter):
    """保留"开始就已亮"的 emitter。

    判定依据是 emitter 的首次出现帧 on_frame 是否落在指定帧号窗口内。

    Parameters
    ----------
    frame_range : (start, end) | None
        帧号窗口(含两端)。保留 on_frame ∈ [start, end] 的 emitter。
        start / end 任一端可为 None 表示该端不设限, 例如:
          (0, 10)     -> 前 10 帧内就亮的分子
          (None, 50)  -> 50 帧及以前出现的分子
          (5, None)   -> 第 5 帧及以后才出现的分子
        None (默认) -> 沿用旧行为: 以数据中最小的 on_frame 为界。
    """
    def __init__(self, frame_range=None):
        self.frame_range = frame_range

    def keep(self, emitters):
        if not emitters:
            return np.array([], bool)
        on = np.array([e.on_frame for e in emitters], dtype=int)

        if self.frame_range is None:                     # 向后兼容: 只保留最早出现的
            return on == on.min()

        try:
            start, end = self.frame_range
        except (TypeError, ValueError):
            raise ValueError(
                "frame_range 需要能解包成两项的 (start, end), "
                "例如 (0, 10); 任一端可用 None 表示不设限。") from None

        m = np.ones(on.size, bool)
        if start is not None:
            m &= on >= start
        if end is not None:
            m &= on <= end
        return m


def ecc_pvalue(ecc, n):
    """各向同性零假设下, 2×2 散度矩阵特征值比的**精确**上尾概率。

    推导: 二维各向同性点云的散度矩阵服从 W_2(I, n−1)。其两个特征值之比 ``t``
    的密度, 在变量替换 ``u = (t−1)/(t+1)`` 之后是 ``ν·u·(1−u²)^(ν/2−1)``
    (``ν = n−2``), 于是上尾有闭式

        P(T > t) = (1 − u²)^(ν/2) = ( 4t / (t+1)² )^((n−2)/2)

    实测(200000 次 Monte Carlo, 各向同性)在 t = 2..8 上与公式吻合到 1% 以内,
    例如 n=20, t=4.90 时公式 0.57%、实测 0.60%。

    Parameters
    ----------
    ecc : float
        特征值比 λmax/λmin (``Emitter.eccentricity()``)。
    n : int
        定位点数。
    """
    n = int(n)
    if n < 3 or not np.isfinite(ecc) or ecc < 1.0:
        return 1.0
    q = 4.0 * float(ecc) / (float(ecc) + 1.0) ** 2
    return float(min(1.0, q ** ((n - 2) / 2.0)))


class EccentricityFilter(EmitterFilter):
    """用点云的各向异性剔除"两个 PSF 被拟合成一个"的假轨迹。

    两个分子同时亮、间距小于分辨率时, 单 PSF 拟合给的是强度加权质心; 相对亮度
    随时间起伏, 质心就在两点之间来回走 —— 整条点云沿连线被拉长。单一分子的点云
    则是各向同性的。

    判据用 :func:`ecc_pvalue` 的**精确零分布**, 而不是"若干倍的期望比值":
    后者在短轨迹上会误杀真实分子 (实测 n=8 时 11.9%、n=10 时 6.5%), 而精确 p 值
    在任何 n 下的假阳性都等于 ``p_crit``。

    注意: 统计量必须用**协方差**的特征值比。若先按轴归一化再算 (等价于用相关
    矩阵), 沿坐标轴方向的拉长会被完全抹掉 —— 两分子沿 x 排开时相关系数为 0,
    归一化后就是一团圆云。

    Parameters
    ----------
    p_crit : float
        p 值阈值, 小于它才**进入候选** (默认 1e-3)。p 值随样本量增长而变小,
        单独用它会误杀长轨迹 —— 见 ``ecc_min``。
    min_pts : int
        少于这么多定位点不做判定。短轨迹检验功效低, 只是"判不出来",
        不会被放行成假阳性。
    ecc_min : float | None
        **效应量下限**: 只有轴比也超过它才真的剔除, ``None`` = 只看 p 值
        (旧行为)。加这一条是因为 p 值回答"是否各向异性", 而我们要剔的是
        "拉长到物理上有意义" —— 两者在长轨迹上会脱节:

        实测 (σ_loc=0.15 px, n=828)::

            单分子 + 0.2 px 残余漂移 -> ecc 1.16 (p 显著与否都该保留)
            单分子 + 0.3 px 残余漂移 -> ecc 1.35, 纯 p 值拒收 73%
            两分子间距 0.5 px        -> ecc 1.94, 必须剔除

        纯 p 值在 n=828 时只要轴比 1.30 就判假, 于是 0.3 px 的残余漂移
        (远小于分辨率) 会把整条真轨迹杀掉; 红通道轨迹短 (n~250) 反而过得去,
        结果是"轨迹越长越容易被误杀"。``ecc_min=1.5`` 正好落在两者之间。
    """

    def __init__(self, p_crit=1e-3, min_pts=8, ecc_min=1.5):
        self.p_crit = float(p_crit)
        self.min_pts = int(min_pts)
        self.ecc_min = None if ecc_min is None else float(ecc_min)

    def keep(self, emitters):
        m = np.ones(len(emitters), bool)
        for k, e in enumerate(emitters):
            if e.n_locs < self.min_pts:
                continue
            ecc = e.eccentricity()
            if self.ecc_min is not None and ecc <= self.ecc_min:
                continue                       # 效应量不够, 不剔
            if ecc_pvalue(ecc, e.n_locs) < self.p_crit:
                m[k] = False
        return m