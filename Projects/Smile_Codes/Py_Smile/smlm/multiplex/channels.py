# smlm/multiplex/channels.py
"""把人工取的框变成 ``{通道: 窗口}``, 再从一幅双视角图里裁出 ``MultiImageStack``。

只处理**已经知道该取哪一块**之后的机械步骤:

* ``analysis_crops(calib, ...)``: 从 ``DualViewCalibration.roi`` 里取出标定时用的
  窗口, 加/减漂移余量, 红通道由"绿 + split_x"派生。
* ``multi_from_crops(stack, crops, reg=...)``: 按框整数裁切, 拼成
  :class:`~smlm.datatypes.MultiImageStack`; 粗偏移在整数裁切时就被吸收, 不插值。

新流程(标定与分析分离)走 ``MultichannelROI`` + ``bundle``, 见
``smTest_dual.ipynb`` / ``smTest_dual_calibration.ipynb``; 这里保留的是
"手上只有一个旧标定结果 + 一堆 sif"时(批量处理)还要用的最小集合。
"""
from __future__ import annotations

import numpy as np

from ..datatypes import MultiImageStack
from ..imgprocess import Crop

DEFAULT_SPLIT_X = 256


def windows_from_crops(crops, reference='g', moving='r'):
    """把人工指定的框转成 ``{通道: (y0, x0, h, w)}``。

    ``crops = {'g': dict(y=(y0, y1), x=(x0, x1)), 'r': {...}}``
    (``MultiImageStack.crop()`` 记进 ``metadata['analysis_crops']`` 的就是这个格式)。
    通道数不限: ``crops`` 里有几个通道就处理几个(``reference`` 必须在其中)。
    """
    if reference not in crops:
        raise ValueError(f'crops 缺少参考通道 {reference!r}: {list(crops)}')
    order = [c for c in (reference, moving) if c in crops]
    order += [c for c in crops if c not in order]
    out = {}
    for ch in order:
        c = crops[ch]
        y0, y1 = (int(v) for v in c['y'])
        x0, x1 = (int(v) for v in c['x'])
        if y1 <= y0 or x1 <= x0:
            raise ValueError(f'通道 {ch!r} 的框无效: y={c["y"]} x={c["x"]}')
        out[ch] = (y0, x0, y1 - y0, x1 - x0)
    return out


def multi_from_crops(stack, crops, *, reference='g', reg=None, margin=0,
                     pixel_size=None, roi_dict=None, verbose=True):
    """按人工指定的框切出**多通道**并堆成 ``MultiImageStack``。

    ``crops`` 与 :func:`analysis_crops` 是同一套格式::

        {'g': dict(y=(y0, y1), x=(x0, x1)), 'r': {...}}

    Parameters
    ----------
    reg : dict, optional
        ``{通道: DualViewCalibration}``(参考通道不写)。存进
        ``MultiImageStack.reg``, 之后 ``multi.geometry()`` 就能还原全部坐标换算。
    roi_dict : dict, optional
        记进 ``metadata['roi']``(通常就是 ``calib.roi``), 供几何对象还原 margin。
    """
    w = windows_from_crops(crops, reference)
    shapes = {(h, ww) for _, _, h, ww in w.values()}
    if len(shapes) != 1:
        raise ValueError(f'各通道框尺寸不一致 {sorted(shapes)}; '
                         f'MultiImageStack 要求同形')
    channels = list(w)
    data = np.stack([np.asarray(Crop(y=(y0, y0 + h), x=(x0, x0 + ww))(stack).data)
                     for ch, (y0, x0, h, ww) in w.items()], axis=0)
    offs = np.array([(y0, x0) for y0, x0, _, _ in w.values()], dtype=int)
    meta = {**dict(stack.metadata), 'roi': dict(roi_dict or {}),
            'margin': int(margin)}
    if verbose:
        for ch, (y0, x0, h, ww) in w.items():
            print(f'  {ch}: detector y{y0}..{y0 + h} x{x0}..{x0 + ww} ({h}x{ww})')
    return MultiImageStack(
        data=data, channels=channels, reference=channels.index(reference),
        reg=dict(reg or {}), offsets=offs,
        pixel_size=pixel_size if pixel_size is not None else stack.pixel_size,
        metadata=meta)


def mirror_crop(crop, split_x, axis='x'):
    """把一个通道的框搬到另一个通道: 探测器坐标 + ``split_x``。"""
    y0, y1 = (int(v) for v in crop['y'])
    x0, x1 = (int(v) for v in crop['x'])
    s = int(split_x)
    if axis == 'x':
        return dict(y=(y0, y1), x=(x0 + s, x1 + s))
    return dict(y=(y0 + s, y1 + s), x=(x0, x1))


def analysis_crops(calib, x_trim=0, y_trim=0, split_x=None):
    """从**标定取框**派生分析取框: 绿 = 标定绿框向内收 trim, 红 = 绿 + split_x。

    x 方向直接沿用标定取框(那组框已经让开了鬼影带), 所以分析窗口一定落在标定
    覆盖范围内, 且两通道覆盖同一块物理视场 —— 公共可用区就是窗口本身。
    想在标定取框内再往里收, 就调大 ``x_trim``(两侧对称, 对齐关系不变)。
    """
    roi = calib.roi or {}
    sx = int(roi.get('split_x', split_x or DEFAULT_SPLIT_X))
    gx = roi.get('green_x')
    rows = roi.get('ref_rows', roi.get('row_range'))
    if gx is None or rows is None:
        raise ValueError('calib.roi 缺 green_x / ref_rows, 无法派生分析取框')
    y0, y1 = int(rows[0]) + int(y_trim), int(rows[1]) - int(y_trim)
    x0, x1 = int(gx[0]) + int(x_trim), int(gx[1]) - int(x_trim)
    if y1 <= y0 or x1 <= x0:
        raise ValueError(f'分析取框为空: y={(y0, y1)} x={(x0, x1)}')
    g = dict(y=(y0, y1), x=(x0, x1))
    return g, mirror_crop(g, sx)
