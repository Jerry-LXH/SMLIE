from .base import TraceExtractor
from ..datatypes import Trace
import numpy as np


class ApertureTraceExtractor(TraceExtractor):
    def __init__(self, edge=3, do_bg_subtract=True, bg_edge=None,
                 bg_method='median', excl_radius=None,
                 signal_shape='circle',
                 sigma_clip_n=3.0, sigma_clip_iter=3, verbose=True):
        self.edge = int(edge)
        self.do_bg = bool(do_bg_subtract)
        self.bg_edge = int(bg_edge) if bg_edge is not None else self.edge + 3
        self.bg_method = bg_method.lower()
        self.excl_radius = int(excl_radius) if excl_radius is not None else self.edge
        self.signal_shape = signal_shape.lower()
        self.sc_n = float(sigma_clip_n)
        self.sc_iter = int(sigma_clip_iter)
        self.verbose = verbose
        if self.bg_edge <= self.edge:
            raise ValueError("bg_edge must be larger than edge.")
        if self.signal_shape not in ('circle', 'square'):
            raise ValueError("signal_shape must be 'circle' or 'square'.")

    # ---------- 主入口 ----------
    def extract(self, emitters, movie, n_frames=None) -> list[Trace]:
        movie = np.asarray(movie)
        F, rows, cols = movie.shape
        if n_frames is None:
            n_frames = F
        frames_axis = np.arange(n_frames)

        # Emitter.x / .y 是整条轨迹的数组; 用质心作为孔径中心。
        # centroid = [x_mean, y_mean], 这里需要 [row=y, col=x]。
        centers = np.array([[e.centroid[1], e.centroid[0]] for e in emitters],
                           dtype=float)                       # [row(y), col(x)]
        # Emitter 无 id 字段, 用枚举索引作为 trace 标识。
        ids = list(range(len(emitters)))
        N = len(emitters)

        if self.verbose:
            self._print_config(N, rows, cols, F)

        traces = []
        valid_counts = []

        for i in range(N):
            ry, cx = centers[i]
            r, c = int(round(ry)), int(round(cx))

            # ---- 信号积分 (圆/方孔径, 边界裁剪) ----
            raw, area = self._integrate_signal(movie, r, c, rows, cols)

            if not self.do_bg:
                traces.append(Trace(ids[i], frames_axis, raw.copy(), raw,
                                    None, np.zeros(n_frames), cx, ry,
                                    {'edge': self.edge, 'shape': self.signal_shape,
                                     'area': area, 'do_bg': False}))
                continue

            # ---- 背景环 + 邻居屏蔽 ----
            bg_stack, n_valid = self._background_pixels(
                movie, r, c, i, centers, rows, cols)
            valid_counts.append(n_valid)

            bg_pp, bg_std = self._estimate_background(bg_stack)  # (F,), (F,)

            intensity = raw - bg_pp * area
            snr = np.where((intensity > 0) & (bg_std > 0),
                           intensity / np.sqrt(intensity + area * bg_std**2),
                           0.0)

            traces.append(Trace(
                ids[i], frames_axis, intensity, raw, bg_pp, snr, cx, ry,
                meta={'edge': self.edge, 'bg_edge': self.bg_edge,
                      'shape': self.signal_shape, 'method': self.bg_method,
                      'area': area, 'n_bg_pixels': int(n_valid),
                      'excl_radius': self.excl_radius, 'do_bg': True}))

        if self.verbose:
            print(f"Done. {N} traces x {n_frames} frames extracted.")
            if valid_counts:
                vc = np.array(valid_counts)
                print(f"  Background pixels per emitter: "
                      f"min={vc.min()}, median={int(np.median(vc))}, max={vc.max()}")
        return traces

    # ---------- 信号积分 ----------
    def _integrate_signal(self, movie, r, c, rows, cols):
        rr, cc, mask = self._aperture_mask(r, c, self.edge, rows, cols,
                                           self.signal_shape)
        if mask.sum() == 0:
            F = movie.shape[0]
            return np.zeros(F), 0
        sub = movie[:, rr[0]:rr[1]+1, cc[0]:cc[1]+1].astype(float)  # (F,h,w)
        raw = sub[:, mask].sum(axis=1)                              # (F,)
        return raw, int(mask.sum())

    # ---------- 背景像素栈 ----------
    def _background_pixels(self, movie, r, c, idx, centers, rows, cols):
        r1, r2 = max(0, r - self.bg_edge), min(rows - 1, r + self.bg_edge)
        c1, c2 = max(0, c - self.bg_edge), min(cols - 1, c + self.bg_edge)
        rr = np.arange(r1, r2 + 1)[:, None]
        cc = np.arange(c1, c2 + 1)[None, :]
        dist = np.sqrt((rr - r)**2 + (cc - c)**2)

        # 圆环: 信号半径外、bg_edge 半径内
        valid = (dist > self.edge) & (dist <= self.bg_edge)

        # 剔除所有邻居 (圆形排除)
        for m in range(len(centers)):
            if m == idx:
                continue
            rn, cn = int(round(centers[m, 0])), int(round(centers[m, 1]))
            # 快速跳过远处邻居
            if (rn < r1 - self.excl_radius or rn > r2 + self.excl_radius or
                    cn < c1 - self.excl_radius or cn > c2 + self.excl_radius):
                continue
            dn = np.sqrt((rr - rn)**2 + (cc - cn)**2)
            valid &= (dn > self.excl_radius)

        sub = movie[:, r1:r2+1, c1:c2+1].astype(float)  # (F,h,w)
        bg_stack = sub[:, valid]                         # (F, P)
        return bg_stack, int(valid.sum())

    # ---------- 稳健背景估计 (逐帧 sigma-clip, 向量化) ----------
    def _estimate_background(self, bg_stack):
        F = bg_stack.shape[0]
        if bg_stack.shape[1] == 0:
            return np.zeros(F), np.zeros(F)

        clipped = self._sigma_clip_rows(bg_stack)         # nan 表示被剔除
        n_keep = np.isfinite(clipped).sum(axis=1)

        # 主估计
        if self.bg_method == 'median':
            bg_pp = np.nanmedian(clipped, axis=1)
        elif self.bg_method == 'mean':
            bg_pp = np.nanmean(clipped, axis=1)
        elif self.bg_method == 'trimmean':
            bg_pp = self._trimmean_rows(clipped, 20)
        else:
            raise ValueError(f"Unknown bg_method: {self.bg_method}")
        bg_std = np.nanstd(clipped, axis=1)

        # fallback: 有效像素 <5 的帧改用未裁剪中位数
        fb = n_keep < 5
        if fb.any():
            bg_pp[fb] = np.nanmedian(bg_stack[fb], axis=1)
            bg_std[fb] = np.nanstd(bg_stack[fb], axis=1)

        return np.nan_to_num(bg_pp), np.nan_to_num(bg_std)

    # ---------- 工具 ----------
    @staticmethod
    def _aperture_mask(r, c, radius, rows, cols, shape):
        r1, r2 = max(0, r - radius), min(rows - 1, r + radius)
        c1, c2 = max(0, c - radius), min(cols - 1, c + radius)
        rr = np.arange(r1, r2 + 1)[:, None]
        cc = np.arange(c1, c2 + 1)[None, :]
        if shape == 'circle':
            mask = (rr - r)**2 + (cc - c)**2 <= radius**2 + 1e-9
        else:  # square
            mask = np.ones((r2 - r1 + 1, c2 - c1 + 1), dtype=bool)
        return (r1, r2), (c1, c2), mask

    def _sigma_clip_rows(self, arr):
        """对 (F,P) 逐行独立 sigma-clip, 剔除位置置 nan。"""
        x = np.where(np.isfinite(arr), arr, np.nan).astype(float)
        for _ in range(self.sc_iter):
            med = np.nanmedian(x, axis=1, keepdims=True)
            std = np.nanstd(x, axis=1, keepdims=True)
            std = np.where(std < 1e-10, np.inf, std)  # std≈0 时不再剔除
            keep = np.abs(x - med) <= self.sc_n * std
            new = np.where(keep, x, np.nan)
            if np.array_equal(np.isnan(new), np.isnan(x)):
                break
            x = new
        return x

    @staticmethod
    def _trimmean_rows(arr, percent):
        """逐行截尾均值 (总共去掉 percent%%, 两端各半)。"""
        p = percent / 2.0
        lo = np.nanpercentile(arr, p, axis=1, keepdims=True)
        hi = np.nanpercentile(arr, 100 - p, axis=1, keepdims=True)
        m = (arr >= lo) & (arr <= hi)
        xm = np.where(m, arr, np.nan)
        return np.nanmean(xm, axis=1)

    def _print_config(self, N, rows, cols, F):
        print("=" * 40)
        print("  ApertureTraceExtractor: Config")
        print("=" * 40)
        print(f"  Emitters:        {N}")
        print(f"  Image size:      {rows} x {cols} x {F} frames")
        print(f"  Signal ROI:      {self.signal_shape}, radius/edge = {self.edge}")
        if self.do_bg:
            print(f"  Background:      ON (annulus, bg_edge = {self.bg_edge})")
            print(f"  Method:          '{self.bg_method}', "
                  f"excl_radius = {self.excl_radius}")
            print(f"  Sigma-clip:      {self.sc_n}-sigma, {self.sc_iter} iter")
        else:
            print("  Background:      OFF (raw integration only)")
        print("=" * 40 + "\n")