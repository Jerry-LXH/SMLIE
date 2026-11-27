# smlm/fret/datatypes.py
"""FRET 子包自己的数据类型。

放这里而不是 ``smlm/datatypes.py``: 这些类型只服务 FRET, 别的子包不会用到,
放全局反而会稀释 "datatypes 是通用数据" 这条约定。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


# ---------------------------------------------------------------------------
#  事件
# ---------------------------------------------------------------------------
@dataclass
class StepEvent:
    """一次「受体漂白 + 供体去淬灭」候选台阶。

    单位都是 photon/frame。``IA_*`` 是受体(acceptor)通道, ``ID_*`` 是供体(donor)。
    """
    emitter_id: int
    frame: int                      # 台阶帧 j
    lag: int                        # 供体一侧的滞后帧数
    IA_pre: float
    IA_post: float
    ID_pre: float
    ID_post: float                  # 供体亮态电平 (抗闪烁: 长窗内亮帧的中位数)
    dIA: float                      # 红降(含 α·ΔI_D 的反向贡献)
    dID: float                      # 绿升(抗闪烁口径, 用亮态电平)
    dID_long: float                 # 绿升(长窗**纯中位数**口径, 仅供诊断: 闪烁会让它偏低)
    n_on_d: int = -1                # 长窗内处于亮态的帧数
    sig_a: float = np.nan           # 受体通道白噪声 σ
    sig_d: float = np.nan           # 供体通道白噪声 σ
    score: float = np.nan
    flags: dict = field(default_factory=dict)

    # 兼容 dict 风格读取, 方便按名字批量判据
    def __getitem__(self, k):
        if k in self.flags:
            return self.flags[k]
        return getattr(self, k)

    def get(self, k, default=None):
        try:
            return self[k]
        except (KeyError, AttributeError):
            return default

    @property
    def gamma_ratio(self):
        """逐分子 ΔI_A/ΔI_D (未加回 α)。"""
        return self.dIA / self.dID if self.dID else np.nan

    def __repr__(self):
        return (f'StepEvent(#{self.emitter_id}, j={self.frame}, lag={self.lag}, '
                f'ΔI_A={self.dIA:.0f}, ΔI_D={self.dID:.0f})')


# 判据名字(顺序即报表顺序)
CRITERIA = ('pre_a', 'off_a', 'red', 'grn', 'hold', 'post_d', 'drop_d')
CRITERIA_STRICT = CRITERIA[:6]      # 不含绝对幅度下限


def passed(ev: StepEvent, keys=CRITERIA) -> bool:
    return all(bool(ev.flags.get(k, False)) for k in keys)


# ---------------------------------------------------------------------------
#  FRET 结果
# ---------------------------------------------------------------------------
@dataclass
class FretResult:
    """一个分子的 FRET 结果(一行 CSV)。"""
    emitter_id: int
    frame: int
    lag: int = -1
    x: float = np.nan
    y: float = np.nan
    n_pre: int = -1
    n_on: int = -1
    IA_pre: float = np.nan
    ID_pre: float = np.nan
    IA_post: float = np.nan
    ID_post: float = np.nan
    S_A: float = np.nan
    E: float = np.nan
    E_don: float = np.nan
    snr_A: float = np.nan
    t_bleach: float = np.nan

    def as_row(self):
        return {'idx': self.emitter_id, 't_bleach_s': self.t_bleach,
                'frame_j': self.frame, 'lag': self.lag,
                'x_px': self.x, 'y_px': self.y,
                'n_pre': self.n_pre, 'n_on': self.n_on,
                'I_A_pre': self.IA_pre, 'I_D_pre': self.ID_pre,
                'I_A_post': self.IA_post, 'I_D_post': self.ID_post,
                'S_A': self.S_A, 'E': self.E, 'E_don': self.E_don,
                'snr_A': self.snr_A}


# ---------------------------------------------------------------------------
#  复核
# ---------------------------------------------------------------------------
@dataclass
class ReviewDecision:
    """人工复核的一条结论。"""
    emitter_id: int
    keep: bool
    note: str = ''
