# smlm/drift/base.py
"""漂移**估计**侧基类。

`DriftCorrector` (校正侧, 图像/定位) 已移到 `smlm/imgprocess/temporal.py`,
本模块只专注漂移估计。
"""
from __future__ import annotations
from abc import ABC
import numpy as np
from ..datatypes import Localizations, DriftTrace


class DriftEstimator(ABC):
    """漂移估计器基类。

    输入类型**不再限定为 Localizations**, 由子类决定:
      - 点云型 (Localizations): 默认走 `estimate()` -> `estimate_locs()`
        解析成 (pts, frame) 后交给 `_estimate()`;
      - 图像型 (ImageStack, 如后续的基于互相关的估计器): 直接重写
        `estimate()` 接受图像栈即可, 无需实现 `_estimate()`。

    输出统一为 `DriftTrace`: 累积漂移数组 (n_frames, 2), 列顺序 (dy, dx)
    即 (row, col), 与常见 SMLM 约定一致。
    """

    def estimate(self, data, n_frames=None) -> DriftTrace:
        """统一入口。

        data: 默认按 Localizations (点云) 处理; 图像型子类重写本方法
              以接受 ImageStack 等其它输入。
        """
        return self.estimate_locs(data, n_frames)

    def estimate_locs(self, locs: Localizations, n_frames=None) -> DriftTrace:
        """Localizations 专用入口: 解析成点集后调用 _estimate。"""
        if len(locs) == 0:
            raise ValueError("Localizations 为空, 无法估计漂移。")

        frame = np.asarray(locs.frame).astype(int)
        if n_frames is None:
            n_frames = int(frame.max()) + 1        # 假定帧号 0-based

        # 统一成 (row=y, col=x) 点集, 与内部算法约定一致
        pts = np.column_stack([np.asarray(locs.y, float),
                               np.asarray(locs.x, float)])

        delta = self._estimate(pts, frame, n_frames)   # (n_frames, 2) = (dy, dx)
        delta = np.asarray(delta, float)

        return DriftTrace(dx=delta[:, 1].copy(),
                          dy=delta[:, 0].copy(),
                          metadata={'estimator': type(self).__name__})

    def _estimate(self, pts: np.ndarray, frame: np.ndarray,
                  n_frames: int) -> np.ndarray:
        """点云型核心算法。

        pts: (N,2)=(row,col); frame: (N,); 返回 (n_frames,2)=(drow,dcol)。

        图像型子类可忽略本方法 (直接重写 estimate())。
        """
        raise NotImplementedError