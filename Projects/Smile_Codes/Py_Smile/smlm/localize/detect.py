# smlm/localize/detect.py
import time
import numpy as np
from scipy.ndimage import uniform_filter, maximum_filter, convolve

from .base import Detector
from ..datatypes import Candidates   # 注意: types.py 已改名 datatypes.py


class LocalMaxDetector(Detector):
    """自适应阈值局部极大值检测 (参考 MATLAB findMaxima)。

    对每帧: 盒滤波估局部背景 -> 泊松噪声模型定阈值 -> 7x7 局部极大
    -> 边缘屏蔽 -> 邻居计数去孤立噪点。
    """

    def __init__(self, k_sigma=4.0, edge=5, f_em=np.sqrt(2),
                 bg_window=15, max_size=7, min_neighbors=2, verbose=True):
        self.k_sigma = k_sigma          # 信噪比下界, 峰值需 > bg + k_sigma*sigma; 常用 4-6
        self.edge = edge                # 边缘屏蔽宽度 (给后续积分留边界)
        self.f_em = f_em                # EMCCD 过噪因子; 非 EMCCD 设 1
        self.bg_window = bg_window      # 背景盒滤波窗口
        self.max_size = max_size        # 局部极大邻域 (对应 se7)
        self.min_neighbors = min_neighbors  # >=T 的 3x3 邻居数下界, 去单像素噪点
        self.verbose = verbose

    def detect(self, stack):
        data = stack.data
        n_frames, rows, cols = data.shape

        # 边缘 mask 每帧相同, 只算一次
        edge_mask = np.zeros((rows, cols), dtype=bool)
        e = self.edge
        edge_mask[e:rows - e, e:cols - e] = True

        kernel3 = np.ones((3, 3))
        fs, ys, xs = [], [], []

        t0 = time.perf_counter()
        for f in range(n_frames):
            img = data[f].astype(np.float32, copy=False)

            # --- 局部背景 & 散粒噪声阈值 ---
            bg = uniform_filter(img, size=self.bg_window)
            sigma = self.f_em * np.sqrt(np.abs(bg))
            T = bg + self.k_sigma * sigma

            # --- 局部极大 ---
            local_max = img == maximum_filter(img, size=self.max_size)

            # --- 超阈值 & 邻居计数 ---
            above_T = img >= T
            neighbor_count = convolve(above_T.astype(np.float32), kernel3,
                                      mode='constant')

            # --- 综合 mask ---
            mask = (local_max & above_T & edge_mask &
                    (neighbor_count >= self.min_neighbors))

            r, c = np.nonzero(mask)
            if r.size:
                fs.append(np.full(r.size, f, dtype=np.int32))
                ys.append(r.astype(np.int32))
                xs.append(c.astype(np.int32))

            if self.verbose and (f + 1) % 50 == 0:
                print(f"processing frame {f + 1} / {n_frames}")

        if fs:
            frame = np.concatenate(fs)
            y = np.concatenate(ys)
            x = np.concatenate(xs)
        else:
            frame = np.empty(0, np.int32)
            y = np.empty(0, np.int32)
            x = np.empty(0, np.int32)

        if self.verbose:
            dt = time.perf_counter() - t0
            print(f"Total = {dt:.3f} s, per-frame = {dt / n_frames * 1000:.3f} ms")
            print(f"Detection complete, {frame.size / n_frames:.2f} detects/frame.")

        return Candidates(frame=frame, y=y, x=x)