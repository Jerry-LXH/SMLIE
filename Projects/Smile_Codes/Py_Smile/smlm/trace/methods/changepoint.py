from ..base import StateAnalyzer, register_analyzer
import numpy as np

# ======================================================================
# Changepoint 实现
# ======================================================================
def _seg_rss(a, b, cumY, cumY2):
    n = b - a + 1
    s = cumY[b + 1] - cumY[a]
    s2 = cumY2[b + 1] - cumY2[a]
    return max(s2 - s * s / n, 0.0)


def _binary_seg(cumY, cumY2, a, b, penalty, min_seg):
    n = b - a + 1
    if n < 2 * min_seg:
        return []
    rss0 = _seg_rss(a, b, cumY, cumY2)
    if rss0 < 1e-20:
        return []
    best_gain, best_t = -np.inf, -1
    for t in range(a + min_seg - 1, b - min_seg + 1):
        rss1 = max(_seg_rss(a, t, cumY, cumY2) + _seg_rss(t + 1, b, cumY, cumY2), 1e-20)
        gain = n * np.log(rss0 / rss1)
        if gain > best_gain:
            best_gain, best_t = gain, t
    if best_gain > penalty * np.log(n):
        return (_binary_seg(cumY, cumY2, a, best_t, penalty, min_seg)
                + [best_t + 1]
                + _binary_seg(cumY, cumY2, best_t + 1, b, penalty, min_seg))
    return []


def _refine_cp(cumY, cumY2, bounds, min_seg, n_iter, search_half=25):
    bounds = list(bounds)
    if len(bounds) <= 2:
        return bounds
    for _ in range(n_iter):
        changed = False
        for i in range(1, len(bounds) - 1):
            a, b_end = bounds[i - 1], bounds[i + 1] - 1
            lo0, hi0 = a + min_seg, b_end - min_seg + 1
            if lo0 > hi0:
                continue
            curr = bounds[i]
            lo, hi = max(lo0, curr - search_half), min(hi0, curr + search_half)
            best_rss, best_t = np.inf, curr
            for t in range(lo, hi + 1):
                r = _seg_rss(a, t - 1, cumY, cumY2) + _seg_rss(t, b_end, cumY, cumY2)
                if r < best_rss:
                    best_rss, best_t = r, t
            if best_t != bounds[i]:
                bounds[i], changed = best_t, True
        if not changed:
            break
    return bounds


def _agglomerative(seg_mean, seg_len, threshold):
    labels = np.arange(len(seg_mean))
    while True:
        u = np.unique(labels)
        if u.size <= 1:
            break
        cmean = np.array([np.sum(seg_mean[labels == j] * seg_len[labels == j])
                          / np.sum(seg_len[labels == j]) for j in u])
        min_d, ma, mb = np.inf, -1, -1
        for j in range(u.size):
            for k in range(j + 1, u.size):
                d = abs(cmean[j] - cmean[k])
                if d < min_d:
                    min_d, ma, mb = d, u[j], u[k]
        if min_d > threshold:
            break
        labels[labels == mb] = ma
    return labels

@register_analyzer
class ChangepointAnalyzer(StateAnalyzer):
    method_name = "changepoint"

    def param_specs(self):
        o = super(ChangepointAnalyzer, self).param_specs()   # 显式写全, 不用 __class__ 单元
        o.update(
            penalty     = (2.0, "分割惩罚, 越大→切点越少→态越少"),
            min_seg_len = (5,   "最短段长 (帧)"),
            merge_thr   = (2.0, "合并阈值 (× noise)"),
            refine_iter = (5,   "切点位置迭代优化次数"),
        )
        return o

    def _fit_sequence(self, y):
        o = self.opts
        F = y.size
        cumY = np.concatenate([[0.0], np.cumsum(y)])
        cumY2 = np.concatenate([[0.0], np.cumsum(y ** 2)])

        cp = _binary_seg(cumY, cumY2, 0, F - 1, o["penalty"], o["min_seg_len"])
        bounds = sorted(set([0, *cp, F]))
        bounds = _refine_cp(cumY, cumY2, bounds, o["min_seg_len"], o["refine_iter"])

        nseg = len(bounds) - 1
        seg_mean = np.zeros(nseg)
        seg_len = np.zeros(nseg)
        seg_var = np.zeros(nseg)
        for s in range(nseg):
            a, b = bounds[s], bounds[s + 1] - 1
            seg_len[s] = b - a + 1
            seg_mean[s] = (cumY[b + 1] - cumY[a]) / seg_len[s]
            seg_var[s] = _seg_rss(a, b, cumY, cumY2) / max(seg_len[s] - 1, 1)

        noise = np.sqrt(np.sum(seg_var * seg_len) / np.sum(seg_len))
        if noise < 1e-6:
            noise = np.std(y) * 0.1
        labels = _agglomerative(seg_mean, seg_len, o["merge_thr"] * noise)

        seq = np.zeros(F, dtype=int)
        for s in range(nseg):
            seq[bounds[s]:bounds[s + 1]] = labels[s] + 1     # >0
        return seq