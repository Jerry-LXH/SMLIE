# smlm/drift/cloud.py
import numpy as np
from scipy.spatial import cKDTree
from .base import DriftEstimator


class CloudDriftEstimator(DriftEstimator):
    """基于点云最近邻匹配的 XY 漂移估计 (NP-Cloud 简化版)。

    参数
    ----
    frames_per_seg          : 每段帧数, 推荐 10; 闪烁多时 15-30。
    uncertainty             : 定位不确定度 (pixel), 搜索半径 r = r_scale*uncertainty。
                              None 时尝试用 locs.lp 的中位数 / pixel_size 估计。
    matched_pairs_threshold : 最少匹配对数, 不足则跳过该段。
    do_stage2               : 是否执行 Stage 2 (RR-NP Cloud, 全局增强参考)。
    r_scale                 : 搜索半径系数, 默认 4.5 (=1.5*3, 与原实现一致)。
    resample_factor         : Stage 2 参考点云重采样倍数。
    max_iter                : 每段最大迭代次数。
    verbose                 : 打印进度报告。
    """

    def __init__(self, frames_per_seg=10, uncertainty=1.0,
                 matched_pairs_threshold=30, do_stage2=False,
                 r_scale=1.5 * 3, resample_factor=12, max_iter=20,
                 verbose=True):
        self.frames_per_seg = int(frames_per_seg)
        self.uncertainty = float(uncertainty)
        self.matched_pairs_threshold = int(matched_pairs_threshold)
        self.do_stage2 = bool(do_stage2)
        self.r_scale = float(r_scale)
        self.resample_factor = int(resample_factor)
        self.max_iter = int(max_iter)
        self.verbose = verbose

    def _log(self, *a):
        if self.verbose:
            print(*a)

    # -------- 主流程 --------
    def _estimate(self, pts, frame, n_frames):
        fps = self.frames_per_seg
        r = self.r_scale * self.uncertainty
        min_locs = 2 * self.matched_pairs_threshold
        thr = self.matched_pairs_threshold
        end_seg = int(np.ceil(n_frames / fps))

        seg_of = frame // fps                       # 每个定位所属段(0-based)

        # ============ Stage 1: 以第 0 段为参考 ============
        self._log("=== Stage 1: first segment as reference ===")
        Q = pts[seg_of == 0]
        self._log(f"Reference seg 0 (frames 0-{fps-1}): {len(Q)} locs")

        delta_sum = np.zeros((n_frames, 2))          # 列 = (dy, dx)
        sum_trace = np.zeros(2)
        skip_count = 0

        for seg in range(1, end_seg):
            start = seg * fps
            end = min((seg + 1) * fps, n_frames)     # 独占上界
            n_seg = end - start
            mask_r = (frame >= start) & (frame < end)
            src = pts[mask_r]

            skipped = False
            matched = 0
            for it in range(self.max_iter):
                R = src - sum_trace
                if len(R) < min_locs:
                    self._log(f"[Warn] seg {seg} locs {len(R)}<{min_locs}, skip.")
                    skipped = True
                    break

                tree = cKDTree(R)
                nnDist, nnIdx = tree.query(Q, k=1)
                m = nnDist <= r
                matched = int(m.sum())
                if matched < thr:
                    self._log(f"[Warn] seg {seg} pairs {matched}<{thr}, skip.")
                    skipped = True
                    break

                dxy = R[nnIdx[m]] - Q[m]             # 残余漂移
                if it == 0:
                    judge = self.uncertainty / np.sqrt(matched)
                delta = dxy.mean(axis=0)
                sum_trace = sum_trace + delta
                if np.hypot(*delta) < judge:
                    self._log(f"  seg {seg:2d} (f{start}-{end-1}): "
                              f"{it} iters, {matched} pairs")
                    break

            if skipped:
                skip_count += 1
                delta_sum[start:end] = delta_sum[start - 1]
                continue

            frame_drift = (sum_trace - delta_sum[start - 1]) / n_seg
            for f in range(start, end):
                delta_sum[f] = delta_sum[f - 1] + frame_drift

        self._log(f"--- Stage 1: {end_seg-1} segs, {skip_count} skipped, "
                  f"total drift {np.hypot(*delta_sum[-1]):.2f} px ---")

        if not self.do_stage2:
            return delta_sum

        # ============ Stage 2: 全局增强参考 (RR-NP Cloud) ============
        self._log("=== Stage 2: full loc set as reference ===")
        delta_s1 = delta_sum.copy()
        corrected = pts - delta_s1[frame]            # 用 S1 结果预校正全部定位

        # 从各段均匀重采样构建 Q2
        n_target = self.resample_factor * len(Q)
        n_per_seg = int(np.ceil(n_target / end_seg))
        rng = np.random.default_rng(0)
        chunks = []
        for seg in range(end_seg):
            sel = corrected[seg_of == seg]
            if len(sel) > n_per_seg:
                idx = rng.choice(len(sel), n_per_seg, replace=False)
                sel = sel[idx]
            chunks.append(sel)
        Q2 = np.vstack(chunks) if chunks else corrected
        tree2 = cKDTree(Q2)
        self._log(f"Enhanced ref Q2: {len(Q2)} locs")

        seg_drift = np.zeros((end_seg, 2))
        seg_valid = np.ones(end_seg, bool)

        for seg in range(end_seg):
            start = seg * fps
            end = min((seg + 1) * fps, n_frames)
            mid = (start + end - 1) // 2
            mask_r = (frame >= start) & (frame < end)
            src = pts[mask_r]
            trace2 = delta_s1[mid].copy()            # 以 S1 中心帧为初值

            skipped = False
            for it in range(self.max_iter):
                R = src - trace2
                if len(R) < min_locs:
                    skipped = True
                    break
                nnDist, nnIdx = tree2.query(R, k=1)
                m = nnDist <= r
                matched = int(m.sum())
                if matched < thr:
                    skipped = True
                    break
                dxy = R[m] - Q2[nnIdx[m]]
                if it == 0:
                    judge = self.uncertainty / np.sqrt(matched)
                delta = dxy.mean(axis=0)
                trace2 = trace2 + delta
                if np.hypot(*delta) < judge:
                    break

            if skipped:
                seg_valid[seg] = False
                seg_drift[seg] = delta_s1[mid]
            else:
                seg_drift[seg] = trace2

        # 段中心帧位置线性插值到每帧
        seg_mid = np.array([( (s*fps) + (min((s+1)*fps, n_frames)-1) ) / 2
                            for s in range(end_seg)])
        if seg_valid.sum() >= 2:
            for col in (0, 1):
                delta_sum[:, col] = np.interp(
                    np.arange(n_frames),
                    seg_mid[seg_valid], seg_drift[seg_valid, col])
            delta_sum -= delta_sum[0]                # 归一化: 首帧漂移置零
        elif seg_valid.sum() == 1:
            delta_sum = np.repeat(seg_drift[seg_valid], n_frames, axis=0)
            self._log("[Warn] 仅 1 个有效段, 使用常量漂移。")
        else:
            delta_sum = delta_s1
            self._log("[Warn] 无有效段, Stage 2 回退到 Stage 1。")

        self._log(f"--- Stage 2: {seg_valid.sum()}/{end_seg} valid, "
                  f"total drift {np.hypot(*delta_sum[-1]):.2f} px ---")
        return delta_sum