from .base import TraceExtractor
from ..datatypes import Trace
import numpy as np
from scipy.spatial import cKDTree


class ApertureTraceExtractor(TraceExtractor):
    def __init__(self, edge=3, do_bg_subtract=True, bg_edge=None,
                 bg_method='median', excl_radius=None,
                 signal_shape='circle',
                 bg_ref='frames', bg_tail=200,
                 model_locs=None, model_radius=None, self_radius=1.5,
                 model_weight=1.0,
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

        # ---------- 背景: 两种模式 ----------
        # 'frames' : 逐帧环形背景 (默认)。每帧单独用圆环像素 + 邻居掩膜估计背景,
        #            ``excl_radius`` 控制每个邻居被挖掉多大一块。
        #            ⚠ ``excl_radius`` 不要给得太大: 挖得太大等于把"邻居的外周"
        #            也从背景里删掉, 而外周正是最能代表本分子位置上背景水平的
        #            那部分光。一般取 1~2 px (只避开邻居的峰值核心) 即可。
        # 'tail'   : 同一套几何/掩膜, 但只取最后 ``bg_tail`` 帧的估计并聚合成
        #            一个常数, 整条 trace 都扣它 (背景噪声最低, 但要求背景不随时间变)。
        self.bg_ref = str(bg_ref).lower()
        if self.bg_ref not in ('tail', 'frames'):
            raise ValueError("bg_ref must be 'tail' or 'frames'.")
        self.bg_tail = int(bg_tail)
        if self.bg_tail < 5:
            raise ValueError("bg_tail must be >= 5 frames.")

        # ---------- 邻居的模型扣除 (可选) ----------
        # model_locs: 带 frame / x / y / sigma / photons 的对象(如 Localizations),
        #   坐标必须和 movie 同一坐标系。每帧把落在信号圆内的**其他**分子按
        #   高斯模型减掉, 这样信号圆不必再屏蔽邻居(屏蔽=乘性丢光, 见下面说明)。
        self.model_locs = model_locs
        self.model_radius = (float(model_radius) if model_radius is not None
                             else float(self.bg_edge) + 4.0)
        self.self_radius = float(self_radius)
        self.model_weight = float(model_weight)
        if model_locs is not None:
            self._m_fr = np.asarray(model_locs.frame).astype(int)
            self._m_x = np.asarray(model_locs.x, float)
            self._m_y = np.asarray(model_locs.y, float)
            self._m_sig = np.clip(np.asarray(model_locs.sigma, float), 0.5, 6.0)
            self._m_ph = np.asarray(model_locs.photons, float)
            self._m_tree = cKDTree(np.column_stack([self._m_x, self._m_y]))
        else:
            self._m_tree = None

    # ---------- 主入口 ----------
    def extract(self, emitters, movie, n_frames=None,
                mask_centers=None) -> list[Trace]:
        """提取孔径 trace。

        Parameters
        ----------
        emitters : sequence of Emitter
            孔径中心, 取自 ``emitter.centroid`` (``[x, y]``)。
        movie : (F, rows, cols) ndarray
        n_frames : int, optional
            取前多少帧, 默认全部。
        mask_centers : None | (M, 2) ndarray | sequence of Emitter
            **背景环做邻居屏蔽时用的中心点**, 顺序为 ``[row(y), col(x)]``。
            默认 ``None`` = 用 ``emitters`` 自己的中心(向后兼容)。

            双通道场景下取绿通道 trace 时必须传**绿通道自己的定位点**:
            孔径中心来自红通道连出来的 emitter, 它完全不知道绿通道里还有哪些
            亮着的供体分子。不传的话这些邻居会漏进背景环 -> 背景被高估 ->
            trace 被系统性低估, 表现为 FRET 效率偏高、
            ``1 - I_D^pre/I_D^post`` 甚至大于 1。
        """
        movie = np.asarray(movie)
        F, rows, cols = movie.shape
        if n_frames is None:
            n_frames = F
        frames_axis = np.arange(n_frames)

        # Emitter.x / .y 是整条轨迹的数组; 用质心作为孔径中心。
        # centroid = [x_mean, y_mean], 这里需要 [row=y, col=x]。
        centers = self._as_centers(emitters)            # [row(y), col(x)]
        # Emitter 无 id 字段, 用枚举索引作为 trace 标识。
        ids = list(range(len(emitters)))
        N = len(emitters)

        # ---- 邻居屏蔽中心: 建一次 KD 树, 每个 emitter 只取附近的点 ----
        self_mask = mask_centers is None
        mc = centers if self_mask else self._as_centers(mask_centers)
        tree = cKDTree(mc) if len(mc) else None
        query_r = self.bg_edge + self.excl_radius + 1.0

        if self.verbose:
            self._print_config(N, rows, cols, F)
            if not self_mask:
                print(f"  Neighbour mask:  {len(mc)} external centres "
                      f"(mask_centers given)")

        traces = []
        valid_counts = []

        for i in range(N):
            ry, cx = centers[i]
            r, c = int(round(ry)), int(round(cx))

            # ---- 信号积分 (圆/方孔径, 边界裁剪) ----
            raw, area = self._integrate_signal(movie, r, c, rows, cols)

            # ---- 邻居的模型扣除 (需要 model_locs) ----
            nb = None
            if self._m_tree is not None:
                off = self._disc_offsets(r, c, rows, cols, cx, ry)
                nb = self._model_neighbour_photons(cx, ry, off, n_frames)

            if not self.do_bg:
                if nb is not None:
                    raw = raw - nb
                traces.append(Trace(ids[i], frames_axis, raw.copy(), raw,
                                    None, np.zeros(n_frames), cx, ry,
                                    {'edge': self.edge, 'shape': self.signal_shape,
                                     'area': area, 'do_bg': False}))
                continue

            # ---- 背景环 + 邻居屏蔽 (两种模式共用) ----
            nb_rc = np.zeros((0, 2))
            if tree is not None:
                near = tree.query_ball_point((ry, cx), query_r)
                if near:
                    nb_rc = mc[np.asarray(near, dtype=int)]
                    # 去掉"自己": 与孔径中心重合的那一个中心点
                    # (外部中心表里没有身份信息, 只能按距离认; 1.5 px 以内都当作自己,
                    #  因为该范围已经在信号孔径 edge=3 之内, 不屏蔽反而会吃掉内圈背景环)
                    d_self = np.hypot(nb_rc[:, 0] - ry, nb_rc[:, 1] - cx)
                    j = int(np.argmin(d_self))
                    if self_mask or d_self[j] <= 1.5:
                        nb_rc = np.delete(nb_rc, j, axis=0)

            bg_stack, n_valid = self._background_pixels(
                movie, ry, cx, rows, cols, nb_rc)
            valid_counts.append(n_valid)

            # ---- 模式 A: 把环形背景聚合成一个常数 (只用最后 bg_tail 帧) ----
            if self.bg_ref == 'tail':
                bg_c, sd_c = self._bg_tail(bg_stack, n_frames)
                bg_pp = np.full(n_frames, bg_c, float)
                bg_std = np.full(n_frames, sd_c, float)
                intensity = raw - bg_c * area
                if nb is not None:
                    nb = np.minimum(nb, np.maximum(intensity, 0.0))
                    intensity = intensity - nb
                _den = np.maximum(intensity + area * bg_std ** 2, 1e-12)
                snr = np.where((intensity > 0) & (bg_std > 0),
                               intensity / np.sqrt(_den), 0.0)
                traces.append(Trace(
                    ids[i], frames_axis, intensity, raw, bg_pp, snr, cx, ry,
                    meta={'edge': self.edge, 'bg_edge': self.bg_edge,
                          'shape': self.signal_shape, 'method': self.bg_method,
                          'bg_ref': 'tail', 'bg_tail': self.bg_tail,
                          'area': area, 'n_bg_pixels': int(n_valid),
                          'excl_radius': self.excl_radius,
                          'bg_pp': float(bg_c), 'bg_disc': float(bg_c) * area,
                          'nb_sub_median': float(np.median(nb)) if nb is not None else 0.0,
                          'nb_sub_max': float(nb.max()) if nb is not None else 0.0}))
                continue

            # ---- 模式 B: 逐帧环形背景 ----
            bg_pp, bg_std = self._estimate_background(bg_stack)  # (F,), (F,)

            intensity = raw - bg_pp * area
            if nb is not None:
                # 邻居模型**不可能**把净信号扣成负数: 扣减量上限 = 实测的净超出部分。
                # 重叠时两个独立的单分子拟合会把公共通量各记一遍, 邻居的 photons 偏大,
                # 而在供体被淬灭、圆内几乎没有信号的时候, 这种过度扣除会直接把读数
                # 压到负值 —— 而这一段正是 FRET 公式要用的 I_D^pre。
                nb = np.minimum(nb, np.maximum(intensity, 0.0))
                intensity = intensity - nb
            # 分母开方前先夹到正数: 背景过扣时 intensity 可能很负, 直接开方会出
            # "invalid value encountered in sqrt" 警告(np.where 两边都会求值)
            _den = np.maximum(intensity + area * bg_std ** 2, 1e-12)
            snr = np.where((intensity > 0) & (bg_std > 0),
                           intensity / np.sqrt(_den), 0.0)

            traces.append(Trace(
                ids[i], frames_axis, intensity, raw, bg_pp, snr, cx, ry,
                meta={'edge': self.edge, 'bg_edge': self.bg_edge,
                      'shape': self.signal_shape, 'method': self.bg_method,
                      'area': area, 'n_bg_pixels': int(n_valid),
                      'excl_radius': self.excl_radius, 'do_bg': True,
                      'nb_sub_median': float(np.median(nb)) if nb is not None else 0.0,
                      'nb_sub_max': float(nb.max()) if nb is not None else 0.0}))

        if self.verbose:
            print(f"Done. {N} traces x {n_frames} frames extracted.")
            if valid_counts:
                vc = np.array(valid_counts)
                print(f"  Background pixels per emitter: "
                      f"min={vc.min()}, median={int(np.median(vc))}, max={vc.max()}")
        return traces

    # ---------- 信号积分 ----------
    def _disc_offsets(self, r, c, rows, cols, cx, cy):
        """信号孔径内每个像素相对**浮点**中心的偏移, 返回 (Npix, 2) = [dx, dy]。"""
        (r1, r2), (c1, c2), mask = self._aperture_mask(r, c, self.edge, rows, cols,
                                                       self.signal_shape)
        if mask.sum() == 0:
            return np.zeros((0, 2))
        rrg, ccg = np.meshgrid(np.arange(r1, r2 + 1), np.arange(c1, c2 + 1),
                               indexing='ij')
        return np.column_stack([(ccg - cx)[mask], (rrg - cy)[mask]])

    def _model_neighbour_photons(self, cx, cy, off, n_frames):
        """每帧落在信号圆内的**邻居**模型光子数, 返回 (n_frames,)。

        邻居 = 模型定位点里离孔径中心 > ``self_radius`` 的那些; 自己那个不扣,
        否则会把自己的信号也减掉。每个邻居按它自己的 (x, y, sigma, photons)
        高斯模型算"落在本圆内的通量占比", 再按帧累加。
        """
        ii = np.asarray(self._m_tree.query_ball_point((cx, cy), self.model_radius), int)
        if ii.size:
            ii = ii[np.hypot(self._m_x[ii] - cx, self._m_y[ii] - cy) > self.self_radius]
        if ii.size == 0 or off.shape[0] == 0:
            return np.zeros(n_frames)
        dx = off[:, 0][:, None] - (self._m_x[ii] - cx)[None, :]
        dy = off[:, 1][:, None] - (self._m_y[ii] - cy)[None, :]
        sig = self._m_sig[ii]
        w = np.exp(-(dx ** 2 + dy ** 2) / (2.0 * sig[None, :] ** 2)).sum(axis=0)
        w = w / (2.0 * np.pi * sig ** 2)          # 该邻居落在圆内的通量占比
        out = np.bincount(self._m_fr[ii], weights=self._m_ph[ii] * w * self.model_weight,
                          minlength=n_frames)
        return out[:n_frames]

    def _integrate_signal(self, movie, r, c, rows, cols):
        rr, cc, mask = self._aperture_mask(r, c, self.edge, rows, cols,
                                           self.signal_shape)
        if mask.sum() == 0:
            F = movie.shape[0]
            return np.zeros(F), 0
        sub = movie[:, rr[0]:rr[1]+1, cc[0]:cc[1]+1].astype(float)  # (F,h,w)
        raw = sub[:, mask].sum(axis=1)                              # (F,)
        return raw, int(mask.sum())

    # ---------- 模式 A: 环形背景聚合成常数 (最后 bg_tail 帧) ----------
    def _bg_tail(self, bg_stack, n_frames):
        """把**环形背景**在最后 ``bg_tail`` 帧上聚合成一个每像素常数。

        用的是和逐帧模式完全相同的几何与掩膜: 信号半径外的圆环、叶掉邻居
        (``excl_radius``)、逐帧 sigma-clip、``bg_method``; 只是不看整条 trace,
        只看最后 200 帧, 然后把这 200 帧的结果取中位数。

        为什么这样做:
          * 分子早就漂白了, 2000 帧的末尾两通道都是暗的, 这一段就是纯背景;
          * 200 帧 x ~100 个环内像素一起平均, 估计噪声比逐帧降一个量级;
          * 背景一旦是常数, ``I = raw − const`` 的噪声就只剩散粒噪声,
            供体被淬灭时窗口均值不会再因为背景的 1% 波动而变成负数。

        返回 ``(bg_pp, bg_std)``: 每像素背景 / 每像素噪声, 都是标量。
        """
        n_frames = int(n_frames)
        t0 = max(0, n_frames - self.bg_tail)
        tail = np.asarray(bg_stack[t0:n_frames], float)
        if tail.size == 0 or tail.shape[1] == 0:
            return 0.0, 0.0
        bg_pp_t, bg_std_t = self._estimate_background(tail)
        bg_pp_t = np.asarray(bg_pp_t, float)
        bg_std_t = np.asarray(bg_std_t, float)
        return (float(np.median(bg_pp_t)) if bg_pp_t.size else 0.0,
                float(np.median(bg_std_t)) if bg_std_t.size else 0.0)

    # ---------- 背景像素栈 ----------
    def _background_pixels(self, movie, ry, cx, rows, cols, nb_rc):
        """取背景环像素。

        ry / cx 用浮点(不再取整), 避免半个像素的环宽误差;
        nb_rc 是**已剔除自己**的邻居中心, 每个都挖掉 excl_radius 半径的圆。
        """
        r, c = int(round(ry)), int(round(cx))
        r1, r2 = max(0, r - self.bg_edge), min(rows - 1, r + self.bg_edge)
        c1, c2 = max(0, c - self.bg_edge), min(cols - 1, c + self.bg_edge)
        rr = np.arange(r1, r2 + 1)[:, None]
        cc = np.arange(c1, c2 + 1)[None, :]
        dist = np.sqrt((rr - ry)**2 + (cc - cx)**2)

        # 圆环: 信号半径外、bg_edge 半径内
        valid = (dist > self.edge) & (dist <= self.bg_edge)

        # 剔除邻居 (圆形排除)。nb_rc 来自 KD 树预筛, 数量已经不多了。
        for rn, cn in np.asarray(nb_rc, dtype=float):
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
    def _as_centers(obj):
        """把 Emitter 序列或 (M, 2) 数组统一成 ``[row(y), col(x)]`` 的浮点数组。"""
        arr = np.asarray(obj)
        if arr.ndim == 2 and arr.shape[1] == 2 and arr.dtype != object:
            return arr.astype(float)
        return np.array([[e.centroid[1], e.centroid[0]] for e in obj], dtype=float)

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
            if self.bg_ref == 'tail':
                print(f"  Aggregation:     one constant per emitter "
                      f"(last {self.bg_tail} frames)")
        else:
            print("  Background:      OFF (raw integration only)")
        print("=" * 40 + "\n")