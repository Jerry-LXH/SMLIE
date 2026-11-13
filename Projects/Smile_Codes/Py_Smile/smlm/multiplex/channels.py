# smlm/multiplex/channels.py
"""双视角几何: 切通道、定义"同一物理窗口"、粗偏移估计、边缘余量检查。

背景
----
探测器上一幅图被分成两半: 左侧是参考通道(绿/供体), 右侧是被变换的通道(红/受体),
两者在探测器坐标下相差约 ``split_x`` 个像素。下游分析需要两通道落在**同一物理
窗口**里, 否则坐标系会差一个 256 px 的大平移, 得让多项式去拟合它 —— 数值条件差,
而且粗偏移会隐藏在系数里, 出错时很难发现。

所以这里做两件事:
  1. ``cut_pair(same_window=True)``: 以参考通道为准定义窗口, 另一通道取同一物理
     区域(平移 ``coarse``)。粗偏移在**整数裁切**时就被吸收了, 不涉及任何插值。
  2. ``refine_coarse_offset()``: 用两通道的平均图做局域互相关, 把粗偏移修正到
     ±1 px, 剩下的亚像素残差交给配准模型(只在坐标层)。

边缘余量
--------
窗口四周要让出 ``margin`` px 给漂移留空间。给少了, 漂移会把真实的分子推出窗口
边界 —— 这种损失是**静默**的(只是少点), 所以 ``check_drift_margin()`` 会明确
报出来, 只在超过时 warning。
"""
from __future__ import annotations

from dataclasses import dataclass
import warnings
import numpy as np

from ..datatypes import ImageStack, MultiImageStack
from ..imgprocess import Crop

DEFAULT_SPLIT_X = 256
DEFAULT_ROW_RANGE = (0, 512)
# 分界线鬼影带的**估计**宽度, 仅用于画图提示; 默认**不**自动裁 ——
# 取框由人看着图决定(见 viz.plot_cut_windows 与 notebook 里的 GREEN_CROP / RED_CROP)。
SEAM_HINT = 14
DEFAULT_SEAM = 0


# ---------------------------------------------------------------------------
#  ROI 描述 (与旧 *_dualcal.npz 里的 roi 字典保持同一套键)
# ---------------------------------------------------------------------------
@dataclass
class DualViewROI:
    """双视角几何。列范围按**探测器坐标**存, 与旧版 npz 一致。"""
    split_x: int = DEFAULT_SPLIT_X
    row_range: tuple = DEFAULT_ROW_RANGE
    n_cols: int = 512

    @property
    def green_x(self):
        return (0, int(self.split_x))

    @property
    def red_x(self):
        return (int(self.split_x), int(self.n_cols))

    def as_dict(self, green_x=None, red_x=None, **extra):
        """转成 dict(存进 npz 用)。

        ``green_x`` / ``red_x`` 可以覆盖, 用来记录**实际用于标定的**窗口
        (即 ``seam_windows()`` 裁过 seam 的那一对)。``split_x`` / ``row_range``
        仍描述探测器的双视角布局。
        """
        d = {'split_x': int(self.split_x),
             'row_range': [int(v) for v in self.row_range],
             'green_x': [int(v) for v in (green_x or self.green_x)],
             'red_x': [int(v) for v in (red_x or self.red_x)]}
        d.update(extra)
        return d

    @classmethod
    def from_dict(cls, d, n_cols=None):
        """从 ``calib.roi`` 还原**探测器布局**(两个半幅, 不含取框)。

        注意: ``roi['green_x']`` / ``roi['red_x']`` 存的是**标定用的取框**,
        不是探测器半幅。半幅始终是 ``(0, split_x)`` 和 ``(split_x, 2*split_x)``。
        """
        sx = int(d.get('split_x', DEFAULT_SPLIT_X))
        if n_cols is None:
            nc = int(d.get('n_cols', 2 * sx))
        else:
            nc = int(n_cols)
        return cls(split_x=sx,
                   row_range=tuple(int(v) for v in
                                   d.get('row_range', DEFAULT_ROW_RANGE)),
                   n_cols=nc)

    @classmethod
    def from_stack(cls, stack, split_x=DEFAULT_SPLIT_X, row_range=None):
        _, h, w = stack.data.shape
        return cls(split_x=int(split_x),
                   row_range=tuple(row_range) if row_range else (0, h),
                   n_cols=int(w))

    def check(self, stack):
        """几何自洽性检查: 宽度必须是 2*split_x, 高度要够 row_range。"""
        _, h, w = stack.data.shape
        if w != 2 * int(self.split_x):
            warnings.warn(f'视场宽 {w} != 2*split_x = {2 * self.split_x};'
                          ' 双视角布局可能不对(或做过额外裁剪)', RuntimeWarning,
                          stacklevel=2)
        if h < int(self.row_range[1]):
            warnings.warn(f'视场高 {h} < row_range 上限 {self.row_range[1]};'
                          ' 配准结果会外推', RuntimeWarning, stacklevel=2)
        return self


# ---------------------------------------------------------------------------
#  切通道
# ---------------------------------------------------------------------------
def windows_from_crops(crops, reference='g', moving='r'):
    """把人工指定的框转成 ``{通道: (y0, x0, h, w)}``。

    ``crops = {'g': dict(y=(y0, y1), x=(x0, x1)), 'r': {...}}``
    与 ``coords.geometry_from_crops`` 用**同一种**格式, 两边可以照抄。
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


def split_from_crops(stack, crops, *, reference='g', moving='r', verbose=False):
    """按人工指定的两个框切出两通道。

    取框是**人的判断**, 不应该是算出来的: 分界线那里有一条鬼影/串扰亮带
    (通道图像里的分子可能是另一个通道的尾巴), 但具体到哪一列结束、哪一列开始,
    不同光路、不同焦面都不一样。用 ``viz.plot_cut_windows`` 看一眼再决定。

    Returns
    -------
    (stack_g, stack_r, windows)
        ``windows`` 是 ``{ch: (y0, x0, h, w)}`` (探测器坐标), 可以直接记进 ``calib.roi``。
    """
    w = windows_from_crops(crops, reference, moving)
    stacks = {}
    for ch, (y0, x0, h, ww) in w.items():
        stacks[ch] = Crop(y=(y0, y0 + h), x=(x0, x0 + ww))(stack)
    if verbose:
        for ch, (y0, x0, h, ww) in w.items():
            print(f'  {ch}: detector y{y0}..{y0 + h} x{x0}..{x0 + ww} '
                  f'({h}x{ww})')
    return stacks[reference], stacks[moving], w


def multi_from_crops(stack, crops, *, reference='g', reg=None, margin=0,
                     pixel_size=None, roi_dict=None, verbose=True):
    """按人工指定的框切出**多通道**并堆成 ``MultiImageStack``。

    ``crops`` 与 ``analysis_crops`` / ``geometry_from_crops`` 是同一套格式::

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


def crops_from_ref(crop_ref, split_x, reference='g', moving='r'):
    """人工只给**参考通道**一个框, 另一个通道自动 = 该框 + ``split_x``。

    物理依据: 同一个物点在探测器上出现在 ``f(p)`` 和 ``split_x + f(p)`` 两处,
    所以两通道的窗口只差一个 ``split_x`` 的平移。这样带来的好处:

      * 两个窗口覆盖**同一块物理视场** —— 公共可用区就是框本身, 不会因为
        两个框取得不一样大而在边缘互相丢弃;
      * 两通道 local 坐标几乎相等::

            green_local = det - x0
            red_local   = (det + split_x) - (x0 + split_x) = green_local

        于是配对位移≈0, 多项式系数接近恒等, 出了问题一眼能看出来;
      * 只需要人决定**一个**框, 少一半出错机会。

    注意: 因为两窗口相差 ``split_x``, 而鬼影带骑在分界线上, 所以框的左右界
    要同时避开带子的两侧 —— 绿框右界 ≤ split_x - 半带宽, 红框左界 ≥
    split_x + 半带宽。``viz.plot_cut_windows`` 会把这条带画出来。
    """
    return {reference: dict(crop_ref),
            moving: mirror_crop(crop_ref, split_x)}


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


def seam_windows(roi, seam=0):
    """两个通道的切分窗口, 各自在**分界线那一侧**让出 ``seam`` px。

    为什么必须让: 探测器上两条半幅交界处(本机实测列 246~266)有一圈亮带,
    比平台高 20~70%, 那是两个通道的鬼影/串扰区。落在里面的"分子"可能是另一个
    通道的尾巴 , 放进标定会污染配对。

    注意这会改变两通道 local 坐标的**相对关系**::

        green_local = p                       (绿窗口原点仍在 0)
        red_local   = p - seam                (红窗口原点右移 seam)
        =>  red_local = green_local - seam

    这只是个常数平移, 多项式的常数项会吸收它; 但 ``calib.roi`` 必须记录**裁剪后**
    的窗口, 否则 ``coords.DualViewGeometry`` 的偏移记账会错 seam px。
    """
    s = int(seam)
    gx0, gx1 = (int(v) for v in roi.green_x)
    rx0, rx1 = (int(v) for v in roi.red_x)
    if s > 0:
        if gx1 - s <= gx0 or rx1 - s <= rx0:
            raise ValueError(f'seam={s} px 太大, 半幅宽度只有 {gx1 - gx0} px')
        gx1 -= s        # 绿: 让步 —— 丢掉靠近分界线的右侧
        rx0 += s        # 红: 让步 —— 丢掉靠近分界线的左侧
    return (gx0, gx1), (rx0, rx1)


def split_dual_view(stack, roi=None, *, split_x=DEFAULT_SPLIT_X,
                    row_range=None, seam=0, verbose=False):
    """把整幅双视角图切成两个通道。

    Parameters
    ----------
    seam : int
        两通道在**分界线那侧**各让出多少 px, 避开双视角的重叠/鬼影带
        (见 ``seam_windows``)。标定和分析应当用**同一个小值**, 否则分析窗口会
        比标定覆盖范围更靠外, 变成外推。

    Returns
    -------
    (stack_g, stack_r, roi, offsets)
        offsets 是各通道窗口左上角的**探测器坐标** (row, col);
        用 ``seam_windows`` 可以拿到实际用到的列范围, 记进 ``calib.roi``。
    """
    roi = roi or DualViewROI.from_stack(stack, split_x, row_range)
    roi.check(stack)
    ry = roi.row_range
    gx, rx = seam_windows(roi, seam)
    stack_g = Crop(y=tuple(ry), x=gx)(stack)
    stack_r = Crop(y=tuple(ry), x=rx)(stack)
    offsets = np.array([[int(ry[0]), int(gx[0])],
                        [int(ry[0]), int(rx[0])]], dtype=int)
    if verbose:
        note = f'  (seam={seam} px 已让出分界带)' if seam else ''
        print(f'  green {tuple(stack_g.shape)} {gx}   '
              f'red {tuple(stack_r.shape)} {rx}{note}')
    return stack_g, stack_r, roi, offsets


def cut_pair(stack, roi=None, *, split_x=DEFAULT_SPLIT_X, row_range=None,
             coarse=(0, DEFAULT_SPLIT_X), margin=0, shape=None,
             verbose=True):
    """以参考通道(绿)为准, 裁出一对**同一物理窗口**的通道图。

    Parameters
    ----------
    coarse : (dy, dx)
        被变换通道(红)的窗口相对参考窗口在**探测器坐标**下的偏移。
        典型值 (0, split_x): 因为双视角本来就把红放在右半幅。
    margin : int
        窗口四周让出的边距(px), 给漂移留空间。**手动设置**, 漂移超出会被
        ``check_drift_margin`` 报出来。
    shape : (h, w) 可选
        目标窗口尺寸。给定时以参考通道的可用区域中心取这个尺寸(不够则报错)。

    Returns
    -------
    (data_g, data_r, offsets)
        data_* 是 (T,h,w) 的 ndarray(视图), offsets 是 (2,2) 的探测器坐标,
        offsets[0] / offsets[1] 分别是绿 / 红窗口的左上角。
        两通道形状**完全相同**, 因此可以直接堆成 MultiImageStack。
    """
    roi = roi or DualViewROI.from_stack(stack, split_x, row_range)
    roi.check(stack)
    _, h, w = stack.data.shape
    dy, dx = int(coarse[0]), int(coarse[1])
    m = int(margin)

    # 参考(绿)通道的可用区域
    g_r0, g_r1 = int(roi.row_range[0]) + m, int(roi.row_range[1]) - m
    g_c0, g_c1 = int(roi.green_x[0]) + m, int(roi.green_x[1]) - m

    if shape is not None:
        hh, ww = int(shape[0]), int(shape[1])
        cr = (g_r0 + g_r1) // 2
        cc = (g_c0 + g_c1) // 2
        g_r0, g_r1 = cr - hh // 2, cr - hh // 2 + hh
        g_c0, g_c1 = cc - ww // 2, cc - ww // 2 + ww

    # 红通道窗口 = 同一物理区域 + coarse
    r_r0, r_r1 = g_r0 + dy, g_r1 + dy
    r_c0, r_c1 = g_c0 + dx, g_c1 + dx

    # 越界就把两通道一起往回收(保持物理对应关系, 也保持形状一致)
    lo_r = int(roi.row_range[0])
    hi_r = int(roi.row_range[1])
    shrink = []
    if r_r0 < lo_r:
        shrink.append(lo_r - r_r0)
    if r_c0 < int(roi.red_x[0]):
        shrink.append(int(roi.red_x[0]) - r_c0)
    if r_r1 > hi_r:
        shrink.append(r_r1 - hi_r)
    if r_c1 > int(roi.red_x[1]):
        shrink.append(r_c1 - int(roi.red_x[1]))
    if shrink:
        s = max(shrink)
        g_r0 += s; g_r1 -= s; g_c0 += s; g_c1 -= s
        r_r0 += s; r_r1 -= s; r_c0 += s; r_c1 -= s
        warnings.warn(f'窗口越界, 两通道各向内收缩 {s} px 以保持同一物理视场',
                      RuntimeWarning, stacklevel=2)

    if g_r1 <= g_r0 or g_c1 <= g_c0:
        raise ValueError(f'窗口过小: 行 [{g_r0},{g_r1}) 列 [{g_c0},{g_c1});'
                         f' margin={margin} 可能给得太大了')

    data_g = stack.data[:, g_r0:g_r1, g_c0:g_c1]
    data_r = stack.data[:, r_r0:r_r1, r_c0:r_c1]
    offsets = np.array([[g_r0, g_c0], [r_r0, r_c0]], dtype=int)

    if verbose:
        print(f'  同一物理窗口: green rows[{g_r0},{g_r1}) cols[{g_c0},{g_c1})  '
              f'red rows[{r_r0},{r_r1}) cols[{r_c0},{r_c1})')
        print(f'  -> {data_g.shape} (两通道同形, 粗偏移 {(dy, dx)} 已吸收)')
    return data_g, data_r, offsets


def _half_bounds(roi, ch, reference, n_rows, n_cols):
    """该通道在探测器上允许的 (row0, row1, col0, col1)。"""
    if roi is not None:
        rows = tuple(int(v) for v in roi.row_range)
        if ch == reference:
            return (rows[0], rows[1], 0, int(roi.split_x))
        return (rows[0], rows[1], int(roi.split_x), int(roi.n_cols))
    return (0, int(n_rows), 0, int(n_cols))


def cut_channels(stack, roi=None, offsets=None, *, reference='g', margin=0,
                 shape=None, channels=None, verbose=True):
    """按 ``{通道: (dy, dx)}`` 切出**同一物理窗口**的多通道。

    ``offsets[ch]`` 是相对参考通道窗口原点的位移。两通道时
    ``{'g': (0, 0), 'r': (0, 256)}`` 就等价于旧的双视角切分。

    Returns
    -------
    (data, windows)
        ``data`` 是 ``{通道: (T,H,W) 视图}``, ``windows`` 是
        ``{通道: (y0, x0, h, w)}``; 所有通道的 h/w **完全相同**。
    """
    roi = roi or DualViewROI.from_stack(stack)
    channels = list(channels or (offsets or {}).keys())
    if reference not in channels:
        raise ValueError(f'参考通道 {reference!r} 不在 {channels} 里')
    offs = {reference: (0, 0)}
    offs.update({c: (int(offsets[c][0]), int(offsets[c][1]))
                 for c in channels if c != reference and c in (offsets or {})})
    for c in channels:
        offs.setdefault(c, (0, 0))

    _, n_rows, n_cols = stack.data.shape
    m = int(margin)
    ry0, ry1, rx0, rx1 = _half_bounds(roi, reference, reference, n_rows, n_cols)
    ry0, ry1 = ry0 + m, ry1 - m
    rx0, rx1 = rx0 + m, rx1 - m
    if shape is not None:
        hh, ww = int(shape[0]), int(shape[1])
        cr, cc = (ry0 + ry1) // 2, (rx0 + rx1) // 2
        ry0, ry1 = cr - hh // 2, cr - hh // 2 + hh
        rx0, rx1 = cc - ww // 2, cc - ww // 2 + ww

    def _win(ch):
        dy, dx = offs[ch]
        return (ry0 + dy, rx0 + dx, ry1 + dy, rx1 + dx)

    # 越界就所有通道一起往回收, 保持同一物理视场与相同形状
    shrink = 0
    for ch in channels:
        y0, x0, y1, x1 = _win(ch)
        br0, br1, bc0, bc1 = _half_bounds(roi, ch, reference, n_rows, n_cols)
        shrink = max(shrink, br0 - y0, x0 - bc0, y1 - br1, x1 - bc1)
    if shrink > 0:
        warnings.warn(f'窗口越界, 所有通道向内收缩 {shrink} px 以保持同一物理视场',
                      RuntimeWarning, stacklevel=2)
        ry0 += shrink; ry1 -= shrink; rx0 += shrink; rx1 -= shrink
    if ry1 <= ry0 or rx1 <= rx0:
        raise ValueError(f'窗口过小: 行 [{ry0},{ry1}) 列 [{rx0},{rx1}); '
                         f'margin={margin} 可能给得太大')

    data, windows = {}, {}
    for ch in channels:
        y0, x0, y1, x1 = _win(ch)
        data[ch] = stack.data[:, y0:y1, x0:x1]
        windows[ch] = (int(y0), int(x0), int(y1 - y0), int(x1 - x0))
    if verbose:
        for ch in channels:
            y0, x0, h, w = windows[ch]
            print(f'  {ch}: detector y{y0}..{y0 + h} x{x0}..{x0 + w} ({h}x{w})')
    return data, windows


def build_multi(stack, roi=None, *, channels=None, reference='g',
                offsets=None, coarse=None, margin=0, shape=None, reg=None,
                pixel_size=None, verbose=True):
    """切窗口 -> 堆成 ``MultiImageStack`` (所有通道已对齐到同一物理窗口)。

    Parameters
    ----------
    offsets : dict, optional
        ``{通道: (dy, dx)}``, 相对参考通道窗口原点的位移(探测器坐标)。
        默认按双视角布局 ``{'g': (0, 0), 'r': (0, split_x)}``。
    channels : sequence, optional
        通道名; 默认 ``list(offsets)``。
    reg : dict, optional
        ``{通道: DualViewCalibration}``(参考通道不写), 存进 ``.reg``。
    coarse : tuple, optional
        旧参数 ``(dy, dx)``, 等价于 ``offsets[non_ref]``。
    """
    roi = roi or DualViewROI.from_stack(stack)
    if offsets is None:
        dy, dx = (int(coarse[0]), int(coarse[1])) if coarse is not None \
            else (0, int(roi.split_x))
        others = [c for c in (channels or ['r']) if c != reference]
        offsets = {reference: (0, 0), **{c: (dy, dx) for c in others}}
    channels = list(channels or offsets)
    data_d, windows = cut_channels(stack, roi, offsets, reference=reference,
                                   margin=margin, shape=shape,
                                   channels=channels, verbose=verbose)
    data = np.stack([data_d[c] for c in channels], axis=0)
    offs = np.array([windows[c][:2] for c in channels], dtype=int)
    meta = {**dict(stack.metadata),
            'roi': roi.as_dict(margin=int(margin),
                               offsets={c: list(map(int, v))
                                        for c, v in offsets.items()})}
    return MultiImageStack(
        data=data, channels=channels, reference=channels.index(reference),
        reg=dict(reg or {}), offsets=offs,
        pixel_size=pixel_size if pixel_size is not None else stack.pixel_size,
        metadata=meta)


# ---------------------------------------------------------------------------
#  粗偏移精修 (只修正几 px 的残差)
# ---------------------------------------------------------------------------
def refine_coarse_offset(data_g, data_r, search=8, patch=None, min_var=1e-6):
    """用两通道的平均图做整数 NCC 搜索, 修正 ``coarse`` 的残余偏移。

    这里用**整个内部区域**做互相关(而不是中央 block): 稀疏小球图的中央可能是
    空白, 只取中央 block 会因为没有结构而静默返回 (0, 0)。

    Parameters
    ----------
    data_g, data_r : (T,H,W) 或 (H,W)
        已经切到"同一物理窗口"的两通道图像。理论上只差几个像素。
    search : int
        在 ±search px 内整数搜索。
    patch : int, optional
        只取内部区域中央的 patch×patch 做匹配。None = 用整个内部区域。

    Returns
    -------
    (dy, dx)
        需要**加**到 ``coarse`` 上的整数修正量, 也等于红窗口原点要移动的量。
        推导: 设 ``red_origin - green_origin = coarse + e``。按 ``coarse`` 裁切后
        ``local_red = local_green + e``, 即 ``R[x + e] = G[x]``; 而下面的搜索
        正是找 ``d`` 使 ``R[x + d] = G[x]`` 的 NCC 最大, 故 ``e = d``。
    """
    def _mean(a):
        a = np.asarray(a, float)
        return a.mean(axis=0) if a.ndim == 3 else a

    G = _mean(data_g)
    R = _mean(data_r)
    if G.shape != R.shape:
        raise ValueError(f'两通道形状不一致: {G.shape} vs {R.shape}')
    h, w = G.shape
    m = int(search)
    ih, iw = h - 2 * m, w - 2 * m
    if ih < 8 or iw < 8:
        raise ValueError(f'窗口太小({h}x{w}), 无法在 ±{m} px 内搜索')

    p_h, p_w = (ih, iw) if patch is None else (min(int(patch), ih),
                                               min(int(patch), iw))
    y0 = m + (ih - p_h) // 2
    x0 = m + (iw - p_w) // 2

    ref = G[y0:y0 + p_h, x0:x0 + p_w]
    ref0 = ref - ref.mean()
    ref_norm = float(np.sqrt((ref0 ** 2).sum()))
    if ref_norm <= 0:
        raise ValueError('参考通道图像无对比度(全平), 无法估计粗偏移')

    best = (-np.inf, 0, 0)                    # NCC 越大越好
    for dy in range(-m, m + 1):
        for dx in range(-m, m + 1):
            sub = R[y0 + dy: y0 + dy + p_h, x0 + dx: x0 + dx + p_w]
            s0 = sub - sub.mean()
            den = ref_norm * float(np.sqrt((s0 ** 2).sum()))
            if den <= min_var:
                continue
            ncc = float((ref0 * s0).sum() / den)
            if ncc > best[0]:
                best = (ncc, dy, dx)
    if not np.isfinite(best[0]):
        raise ValueError('NCC 搜索没有有效候选(图像可能全平)')
    return int(best[1]), int(best[2])


# ---------------------------------------------------------------------------
#  边缘余量检查
# ---------------------------------------------------------------------------
def check_drift_margin(dx, dy, margin, *, name='', strict=False, frame0=None):
    """检查漂移是否超出窗口边缘余量 ``margin``。

    漂移超出余量意味着有些帧的真实信号被推出了窗口 —— 这是**静默**的损失,
    所以这里主动报出来。

    Parameters
    ----------
    dx, dy : 长度 n_frames 的累积漂移(px), 与 DriftTrace 的约定一致。
    margin : 手动设置的边距(px)。
    strict : True 时超限直接抛异常, False 只 warning。

    Returns
    -------
    dict(ok, margin, need, excess_dy, excess_dx, name)
        ``need`` 是"刚好不越界"所需的最小边距, 可以直接拿去改 margin。
    """
    dx = np.asarray(dx, float)
    dy = np.asarray(dy, float)
    i0 = 0 if frame0 is None else int(frame0)
    ex_dy = float(np.max(np.abs(dy - dy[i0]))) if len(dy) else 0.0
    ex_dx = float(np.max(np.abs(dx - dx[i0]))) if len(dx) else 0.0
    need = int(np.ceil(max(ex_dy, ex_dx))) + 1
    ok = need <= int(margin)
    info = {'ok': bool(ok), 'name': name, 'margin': int(margin), 'need': need,
            'drift_dy': ex_dy, 'drift_dx': ex_dx,
            'excess': max(0, need - int(margin))}
    if not ok:
        msg = (f'{name or "channel"}: 漂移 ±({ex_dy:.2f}, {ex_dx:.2f}) px 超出'
               f'边缘余量 margin={margin} px, 建议 margin >= {need} px')
        if strict:
            raise RuntimeError(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
    return info
