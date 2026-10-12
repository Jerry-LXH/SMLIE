# smlm/localize/loc/doublehelix.py
"""3D 双螺旋 (Double Helix / DH-PSF) 定位 —— 预留骨架。

双螺旋 PSF 的 z 信息编码在两个旁瓣的旋转角度/间距上, 模型显著不同于高斯
(需要 DH 相位掩模 PSF 或双高斯旁瓣模型)。本模块先给出类骨架占位:

  - 接口与其它 Localizer 一致: localize(stack, cand) -> Localizations;
  - 产物同样复用 Localizations.z / z_lp (单位 nm);
  - 接入具体 DH-PSF 模型时, 参照 astigmatic.py 实现即可 (共享
    common.extract_rois / common.run_pool / fisher.fisher_from_jacobian)。
"""
from __future__ import annotations
from ..base import Localizer


class DoubleHelixLocalizer(Localizer):
    """3D 双螺旋定位器 (预留, 尚未实现)。"""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "DoubleHelixLocalizer 为 3D 双螺旋方案预留, 尚未实现。"
            "请先接入 DH-PSF 模型后使用 (可参考 astigmatic.py 的实现模式)。"
        )

    def localize(self, stack, cand):
        raise NotImplementedError(
            "DoubleHelixLocalizer 为 3D 双螺旋方案预留, 尚未实现。"
        )
