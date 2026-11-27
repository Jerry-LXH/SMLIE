"""变点检测的纯算法 —— 只依赖 numpy, 可以单独 import 做实验。

从旧的 ``StateAnalyzer`` 框架里搬出来, 现在是
:class:`smlm.trace.states.ChangepointStates` 的实现细节。三步:

1. :func:`binary_seg`    自上而下二分, 用 BIC 式判据决定要不要继续切;
2. :func:`refine_cp`     在邻域内重定位每个切点 (二分的切点常偏 1~2 帧);
3. :func:`agglomerative` 按段均值把段聚成**亮度等级** (同一分子的各 ON 段亮度相近)。

为什么需要 ``vst``
------------------
泊松噪声下方差 ∝ 均值, 而 RSS 判据隐含"整条 trace 噪声恒定"。亮段的大抖动会支配
RSS, 暗段的小变化就看不见了。方差稳定变换 (:func:`forward`) 先把噪声拉平再分段,
最后用 :func:`inverse` 把段均值变换回强度单位。

实测 (600 帧 6 个真实切点, 60 次重复, "exact" = 全找到且无多余)::

    对比度   raw exact    Anscombe exact
     1.5x      0.55          0.63
     2.2x      0.80          0.88
     3.8x      0.85          0.98
    10.0x      0.97          0.85        <- 亮分子反而 raw 更好

所以 VST 是个 trade-off, 默认关闭; 暗分子 (亮/暗 < 5x) 才建议开。
"""
from __future__ import annotations

import numpy as np

VST_CHOICES = (None, 'sqrt', 'anscombe')


# ---------------------------------------------------------------------------
#  基础量
# ---------------------------------------------------------------------------
def seg_rss(a, b, cum, cum2):
    """闭区间 ``[a, b]`` 的残差平方和 (用累加和 O(1) 算)。"""
    n = b - a + 1
    s = cum[b + 1] - cum[a]
    s2 = cum2[b + 1] - cum2[a]
    return max(s2 - s * s / n, 0.0)


def robust_noise(y):
    """基于逐帧差分的 MAD 噪声估计 (对阶跃鲁棒)。"""
    y = np.asarray(y, float).reshape(-1)
    if y.size < 2:
        return max(float(np.std(y)) if y.size else 0.0, 1e-6)
    dy = np.diff(y)
    n = np.median(np.abs(dy - np.median(dy))) / 0.6745 / np.sqrt(2)
    if n < 1e-6:
        n = float(np.std(y)) * 0.1
    return max(float(n), 1e-9)


# ---------------------------------------------------------------------------
#  方差稳定变换
# ---------------------------------------------------------------------------
def forward(y, vst=None):
    """强度 -> 变换域。``None`` 恒等; ``'sqrt'`` / ``'anscombe'`` 拉平泊松噪声。"""
    y = np.asarray(y, float).reshape(-1)
    if vst is None:
        return y
    if vst == 'anscombe':
        return 2.0 * np.sqrt(np.maximum(y, 0.0) + 3.0 / 8.0)
    if vst == 'sqrt':
        return np.sqrt(np.maximum(y, 0.0))
    raise ValueError(f'vst 只支持 {VST_CHOICES}, 收到 {vst!r}')


def inverse(z, vst=None):
    """变换域 -> 强度 (反向映射, 保证段均值能换回原始单位)。"""
    z = np.asarray(z, float)
    if vst is None:
        return z
    if vst == 'anscombe':
        return (z / 2.0) ** 2 - 3.0 / 8.0
    if vst == 'sqrt':
        return z ** 2
    raise ValueError(f'vst 只支持 {VST_CHOICES}, 收到 {vst!r}')


# ---------------------------------------------------------------------------
#  三步算法
# ---------------------------------------------------------------------------
def binary_seg(cum, cum2, a, b, penalty, min_seg):
    """自上而下二分分割, 返回切点列表 (切点 = 后一段的起始帧)。

    判据: ``n·log(rss0 / rss1) > penalty·log(n)``。左边是似然比, 右边是对数长度
    的惩罚 —— 段越长, 要求的分割收益越高。
    """
    n = b - a + 1
    if n < 2 * int(min_seg):
        return []
    rss0 = seg_rss(a, b, cum, cum2)
    if rss0 < 1e-20:
        return []
    best_gain, best_t = -np.inf, -1
    for t in range(a + int(min_seg) - 1, b - int(min_seg) + 1):
        rss1 = max(seg_rss(a, t, cum, cum2) + seg_rss(t + 1, b, cum, cum2),
                   1e-20)
        gain = n * np.log(rss0 / rss1)
        if gain > best_gain:
            best_gain, best_t = gain, t
    if best_gain > penalty * np.log(n):
        return (binary_seg(cum, cum2, a, best_t, penalty, min_seg)
                + [best_t + 1]
                + binary_seg(cum, cum2, best_t + 1, b, penalty, min_seg))
    return []


def refine_cp(cum, cum2, bounds, min_seg, n_iter=5, search_half=25):
    """在 ``±search_half`` 内重定位每个切点, 直到不再移动。"""
    bounds = list(bounds)
    if len(bounds) <= 2:
        return bounds
    min_seg = int(min_seg)
    for _ in range(int(n_iter)):
        moved = False
        for i in range(1, len(bounds) - 1):
            a, end = bounds[i - 1], bounds[i + 1] - 1
            lo0, hi0 = a + min_seg, end - min_seg + 1
            if lo0 > hi0:
                continue
            cur = bounds[i]
            lo = max(lo0, cur - int(search_half))
            hi = min(hi0, cur + int(search_half))
            best_rss, best_t = np.inf, cur
            for t in range(lo, hi + 1):
                r = (seg_rss(a, t - 1, cum, cum2)
                     + seg_rss(t, end, cum, cum2))
                if r < best_rss:
                    best_rss, best_t = r, t
            if best_t != bounds[i]:
                bounds[i], moved = best_t, True
        if not moved:
            break
    return bounds


def agglomerative(seg_mean, seg_len, threshold):
    """按均值差把段合并成等级, 返回**按均值升序**的 0-based 标签。

    反复找出均值最接近的两段, 距离小于 ``threshold`` 就并成一个等级。
    """
    seg_mean = np.asarray(seg_mean, float).reshape(-1)
    seg_len = np.asarray(seg_len, float).reshape(-1)
    labels = np.arange(seg_mean.size)

    def _centres(lab):
        u = np.unique(lab)
        return u, np.array([np.average(seg_mean[lab == j], weights=seg_len[lab == j])
                            for j in u])

    while True:
        u, cmean = _centres(labels)
        if u.size <= 1:
            break
        d = np.abs(cmean[:, None] - cmean[None, :])
        np.fill_diagonal(d, np.inf)
        i, j = np.unravel_index(int(np.argmin(d)), d.shape)
        if d[i, j] > threshold:
            break
        labels[labels == u[j]] = u[i]

    u, cmean = _centres(labels)
    out = np.empty_like(labels)
    for rank, pos in enumerate(np.argsort(cmean, kind='stable')):
        out[labels == u[pos]] = rank
    return out.astype(int)


# ---------------------------------------------------------------------------
#  一步到位
# ---------------------------------------------------------------------------
def segment(y, penalty=3.0, min_seg=5, refine_iter=5, merge_thr=2.0,
            vst=None):
    """把一条强度序列分成若干**亮度等级**。

    Parameters
    ----------
    y : array_like
        强度序列 (通常 ``trace['intensity']``)。
    penalty : float
        分割惩罚, 越大 -> 切点越少 -> 等级越少。
    min_seg : int
        最短段长 (帧)。天然抑制单帧 glitch, 所以不需要另外的 min_on / min_off。
    refine_iter : int
        切点精修的迭代轮数。
    merge_thr : float
        段均值合并阈值 (× 稳健噪声)。低于它就认为是同一个亮度等级。
        比较的是各段的**中位数**, 所以单帧 outlier 不会造出假等级。
    vst : None | 'sqrt' | 'anscombe'
        方差稳定变换, 默认不做。

    Returns
    -------
    levels : np.ndarray, (T,) int8
        逐帧亮度等级, ``0`` = 最暗级, 数值越大越亮。
    bounds : np.ndarray, (k + 1,) int
        各段的起止帧 (含起点, 不含终点)。
    """
    y = np.asarray(y, float).reshape(-1)
    n = y.size
    if n == 0:
        return np.zeros(0, np.int8), np.array([0], int)

    min_seg = int(min_seg)
    if n < 2 * min_seg:
        return np.zeros(n, np.int8), np.array([0, n], int)

    z = forward(y, vst)
    cum = np.concatenate(([0.0], np.cumsum(z)))
    cum2 = np.concatenate(([0.0], np.cumsum(z * z)))

    cp = binary_seg(cum, cum2, 0, n - 1, float(penalty), min_seg)
    bounds = refine_cp(cum, cum2, sorted({0, *cp, n}), min_seg,
                       refine_iter)
    bounds = np.asarray(sorted(set(bounds)), int)

    start, stop = bounds[:-1], bounds[1:]
    # 合并判据用**中位数**而不是均值: 单帧 outlier (热像素/宇宙线) 能把一个长段的
    # 均值抬高几个 σ, 凭空造出一个"新等级"; 中位数不受它影响。
    seg_mean = np.array([np.median(y[a:b]) for a, b in zip(start, stop)])
    seg_len = stop - start
    labels = agglomerative(seg_mean, seg_len, float(merge_thr) * robust_noise(y))

    levels = np.empty(n, np.int8)
    for a, b, lab in zip(start, stop, labels):
        levels[a:b] = lab
    return levels, bounds
