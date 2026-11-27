# smlm/multiplex/coords.py
"""坐标层: 分析窗口 local 坐标 <-> 参考通道坐标, 以及公共可用区。

为什么不放进 imgprocess
-----------------------
``imgprocess`` 里全是会**重采样像素**的 ``Transform``。本模块一个像素都不碰 ——
它只搬运坐标。混进图像变换模块会让人以为这里做了插值, 从而怀疑 PSF 被改形。

为什么图像只做漂移、变换只在坐标层
----------------------------------
如果对图像做空间变化的形变重采样, PSF 的有效插值核会随位置变化 -> 同一种分子在
视场不同位置得到不同的有效 PSF -> 逐分子 kappa 抖动。这种抖动**无法**被一个全局
gamma 吸收, 会直接broaden FRET 分布。而坐标是点, 变换它没有任何插值代价。

两套窗口 (关键)
---------------
标定是在**标定窗口**的 local 坐标里拟合的(通常是两通道各用满半幅); 分析时往往还要
再裁一次(避开分界线 / 留漂移余量), 于是有了**分析窗口**。两者位置不同, 所以::

    ch分析local  --+ (分析原点 - 标定原点)-->  ch标定local
    ch标定local  --calib.transform-------->    ref标定local
    ref标定local --- (ref分析原点 - 标定原点)--> ref分析local

漏掉任何一步都会带进一个**静默**的系统偏移(几 px ~ 几十 px)。notebook 里手写的
``OFF_G`` / ``OFF_R`` / ``CORR_OFF_*`` 就是在补这几步; 本模块把它封装掉, 并给出
``report()`` / ``check()`` 做自检。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import warnings
import numpy as np

from .transform import DualViewCalibration


def _as_window(w):
    """窗口 -> (y0, x0, h, w)。接受 (y0,x0,h,w) 或含 y0/y1/x0/x1 的 dict。"""
    if isinstance(w, dict):
        if 'y1' in w and 'x1' in w:
            return (int(w['y0']), int(w['x0']),
                    int(w['y1']) - int(w['y0']), int(w['x1']) - int(w['x0']))
        return (int(w['y0']), int(w['x0']), int(w['h']), int(w['w']))
    y0, x0, h, ww = w
    return (int(y0), int(x0), int(h), int(ww))


def _perimeter(y0, x0, h, w, n=64):
    """矩形窗口四边上均匀取点(探测器坐标), 用于精确求映射后的包围盒。"""
    t = np.linspace(0.0, 1.0, int(n))
    return np.vstack([
        np.column_stack([x0 + t * w, np.full_like(t, y0)]),
        np.column_stack([x0 + t * w, np.full_like(t, y0 + h)]),
        np.column_stack([np.full_like(t, x0), y0 + t * h]),
        np.column_stack([np.full_like(t, x0 + w), y0 + t * h]),
    ])


def calib_windows_from_calibs(calibs, reference):
    """从 ``{通道: calib}`` 还原各通道的**标定窗口** ``{通道: (y0,x0,h,w)}``。

    每个 calib 都是“某个通道 -> reference”的映射, 它自己的 ``roi`` 里记录了
    该通道与参考通道在手写取框时的列范围(``red_x`` / ``green_x``)和行范围。
    多通道时逐个取, 不再假定只有一个 moving。
    """
    out = {}
    for ch, cal in (calibs or {}).items():
        if cal is None or ch == reference or not getattr(cal, 'roi', None):
            continue
        roi = cal.roi
        if roi.get('type') == 'MultichannelROI':
            windows = roi.get('windows')
            if isinstance(windows, dict):
                out.update({name: _as_window(window)
                            for name, window in windows.items()})
                continue
        if 'red_x' not in roi:
            continue
        x0, x1 = (int(v) for v in roi['red_x'])
        rows = roi.get('mov_rows', roi.get('row_range'))
        if rows is None:
            continue
        y0, y1 = int(rows[0]), int(rows[1])
        out[ch] = (y0, x0, y1 - y0, x1 - x0)

    # 参考通道的标定窗口也必须记下来: 漏了它 calib_delta(reference) 就退化成 0,
    # 分析窗口相对标定窗口的行/列偏移会整个丢掉(几十 px 的静默错位)。
    if reference not in out:
        for cal in (calibs or {}).values():
            roi = getattr(cal, 'roi', None) or {}
            if 'green_x' not in roi:
                continue
            x0, x1 = (int(v) for v in roi['green_x'])
            rows = roi.get('ref_rows', roi.get('row_range'))
            if rows is None:
                continue
            y0, y1 = int(rows[0]), int(rows[1])
            out[reference] = (y0, x0, y1 - y0, x1 - x0)
            break
    return out


@dataclass
class DualViewGeometry:
    """双视角几何: 分析窗口 + 标定窗口 + 通道间配准, 提供坐标互转。

    Parameters
    ----------
    windows : dict
        ``{通道名: (y0, x0, h, w)}``, **分析**窗口, 探测器坐标。
    calib : DualViewCalibration | dict | None
        通道 -> 参考通道 的配准。允许三种写法:
          * 单个 ``DualViewCalibration``(2 通道, 向后兼容);
          * ``{通道: DualViewCalibration}``(任意通道数, 参考通道不在字典里);
          * ``None``(只做窗口搬运, 不做跨通道映射)。
        规范化后统一存在 ``self.calibs`` 里, ``self.moving`` 仅在恰好 1 个
        非参考通道时才非 None。
    calib_windows : dict, optional
        标定窗口(探测器坐标)。默认从各 ``calib.roi`` 推。
    reference : str
        参考(不动)通道名, 默认 'g'。
    check_coverage : bool
        分析窗口是否落在标定覆盖范围内。超出 = 外推, 结果不可信 -> warning。
    """
    windows: dict
    calib: object = None
    calib_windows: dict = None
    reference: str = 'g'
    roi: dict = None
    margin: int = 0
    check_coverage: bool = True
    strict: bool = False
    _warned: list = field(default_factory=list, repr=False)

    def __post_init__(self):
        self.windows = {str(k): _as_window(v) for k, v in self.windows.items()}
        if self.reference not in self.windows:
            raise ValueError(f'参考通道 {self.reference!r} 不在 windows 里: '
                             f'{list(self.windows)}')

        # calib 允许三种写法: 单个变换(2 通道) / {通道: 变换}(N 通道) / None
        if self.calib is None:
            self.calibs = {}
        elif isinstance(self.calib, dict):
            self.calibs = {str(k): v for k, v in self.calib.items()
                           if v is not None}
        else:
            name = getattr(self.calib, 'moving', None)
            if name not in self.windows:
                name = next((c for c in self.windows
                             if c != self.reference), None)
            if name is None:
                raise ValueError('给了 calib 但没有非参考通道')
            self.calibs = {name: self.calib}
        unknown = [c for c in self.calibs if c not in self.windows]
        if unknown:
            raise ValueError(f'calib 里的通道 {unknown} 不在 windows 里: '
                             f'{list(self.windows)}')
        if self.reference in self.calibs:
            raise ValueError(f'参考通道 {self.reference!r} 不该有配准(恒等)')

        moving = [c for c in self.windows if c != self.reference]
        self.moving = moving[0] if len(moving) == 1 else None

        if self.calib_windows is not None:
            self.calib_windows = {str(k): _as_window(v)
                                  for k, v in self.calib_windows.items()}
        elif self.calibs:
            self.calib_windows = calib_windows_from_calibs(self.calibs, self.reference)
        if self.calibs and not self.calib_windows:
            warnings.warn('标定里没有 roi, 无法还原标定窗口; 假定分析窗口 = 标定窗口',
                          RuntimeWarning, stacklevel=2)
            self.calib_windows = dict(self.windows)

        for fn in ('_check_shape', '_check_coverage'):
            msg = getattr(self, fn)()
            if not msg:
                continue
            if self.strict:
                raise ValueError(msg)
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            self._warned.append(msg)

    # ---------- 检查 ----------
    @property
    def ref_window(self):
        return self.windows[self.reference]

    def _check_shape(self):
        rg = self.ref_window
        bad = [c for c, w in self.windows.items()
               if c != self.reference and (w[2] != rg[2] or w[3] != rg[3])]
        if not bad:
            return None
        desc = ', '.join(f'{c}: {self.windows[c][2]}x{self.windows[c][3]}'
                         for c in bad)
        return (f'通道窗口形状与参考 {self.reference} ({rg[2]}x{rg[3]}) 不同: '
                f'{desc}; 不能直接堆成 MultiImageStack; 公共可用区按交集算')

    def _check_coverage(self):
        """分析窗口(换算到标定 local 后)是否超出标定点的覆盖范围。

        为什么对 order<=1 不算问题: 仿射是**精确可外推**的 —— 超出跳标定点范围
        也只是把一个平面延伸出去, 不会像高阶多项式那样在边缘乱飞。所以低阶时
        只给一条提示, 不报警; 高阶时才是真风险。
        """
        if not self.calibs or not self.calib_windows:
            return None
        bad = []
        for ch, w in self.windows.items():
            cal = self.calibs.get(ch)
            if cal is None or ch not in self.calib_windows:
                continue
            cov = cal.coverage()
            if not cov:
                continue
            d = self.calib_delta(ch)
            xs = (d[0], d[0] + w[3])
            ys = (d[1], d[1] + w[2])
            over_x = max(0.0, cov['x'][0] - xs[0], xs[1] - cov['x'][1])
            over_y = max(0.0, cov['y'][0] - ys[0], ys[1] - cov['y'][1])
            if over_x > 1 or over_y > 1:
                bad.append((ch, over_x, over_y))
        if not bad:
            return None
        desc = '; '.join(f'{ch} 超出 x {ox:.1f} / y {oy:.1f} px'
                         for ch, ox, oy in bad)
        worst = max(max(ox, oy) for _, ox, oy in bad)
        order = max(int(getattr(c, 'poly_order', 1)) for c in self.calibs.values())
        if order <= 1:
            return (f'分析窗口比标定点覆盖大 {worst:.1f} px; 当前是 order '
                    f'{order}(仿射), 外推精确, 可忽略。   {desc}')
        return f'[!] 分析窗口超出标定覆盖范围(高阶多项式外推不可信): {desc}'

    def check(self):
        """列出所有已知问题(不抛异常)。"""
        return [m for m in (self._check_shape(), self._check_coverage()) if m]

    # ---------- 第一层: 探测器 <-> 分析窗口 local ----------
    def local_to_detector(self, xy, ch=None):
        ch = ch or self.reference
        y0, x0, _, _ = self.windows[ch]
        return np.asarray(xy, float).reshape(-1, 2) + np.array([x0, y0], float)

    def detector_to_local(self, xy, ch=None):
        ch = ch or self.reference
        y0, x0, _, _ = self.windows[ch]
        return np.asarray(xy, float).reshape(-1, 2) - np.array([x0, y0], float)

    def window_bounds(self, ch=None, frame='analysis'):
        """窗口边界 (x0, x1, y0, y1)。"""
        ch = ch or self.reference
        y0, x0, h, w = self.windows[ch]
        if frame == 'detector':
            return (float(x0), float(x0 + w), float(y0), float(y0 + h))
        return (0.0, float(w), 0.0, float(h))

    def window_rect(self, ch=None):
        """``(y0, x0, h, w)`` 分析窗口(探测器坐标)。"""
        return self.windows[ch or self.reference]

    # ---------- 第二层: 分析 local <-> 标定 local ----------
    def calib_delta(self, ch):
        """分析窗口 local 与标定窗口 local 的平移差, (x, y) 顺序。

        即 notebook 里的 ``OFF_G`` / ``OFF_R``。
        """
        if not self.calib_windows or ch not in self.calib_windows:
            return np.zeros(2)
        w, cw = self.windows[ch], self.calib_windows[ch]
        return np.array([float(w[1] - cw[1]), float(w[0] - cw[0])])

    def to_calib_frame(self, xy, ch=None):
        """分析窗口 local -> 标定窗口 local(纯平移)。"""
        ch = ch or self.reference
        return np.asarray(xy, float).reshape(-1, 2) + self.calib_delta(ch)

    def from_calib_frame(self, xy, ch=None):
        """标定窗口 local -> 分析窗口 local(纯平移)。"""
        ch = ch or self.reference
        return np.asarray(xy, float).reshape(-1, 2) - self.calib_delta(ch)

    # ---------- 第三层: 通道 -> 参考通道 ----------
    def to_reference_calib(self, xy, ch=None):
        """ch 分析 local -> **参考通道标定 local**(与标定系数同坐标系)。"""
        ch = ch or self.moving or self.reference
        p = self.to_calib_frame(xy, ch)
        c = self.calibs.get(ch)
        if c is None:
            return p
        return np.asarray(c.transform(p), float).reshape(-1, 2)

    def to_reference(self, xy, ch=None):
        """ch 分析 local -> **参考通道分析 local**(可直接和裁切后的图像比对)。"""
        ch = ch or self.moving or self.reference
        return self.from_calib_frame(self.to_reference_calib(xy, ch), self.reference)

    def to_channel(self, xy, ch=None, n_iter=10, tol=1e-9):
        """``to_reference`` 的反向: 参考分析 local -> ch 分析 local。"""
        ch = ch or self.moving or self.reference
        p = self.to_calib_frame(xy, ch=self.reference)
        if ch == self.reference or ch not in self.calibs:
            return self.from_calib_frame(p, ch)
        return self.from_calib_frame(self._inv_poly(p, ch, n_iter, tol), ch)

    def _inv_poly(self, xy, ch=None, n_iter=10, tol=1e-9):
        """``apply_poly`` 的逆。

        注意: 必须反演**多项式**本身, 不能拿 ``affine`` 代替 —— 两者虽然都拟合同一批
        点, 但一次多项式与仿射是两个独立的最小二乘结果(只差浮点噪声), 混用会让
        往返误差落在 1e-1 px 量级, 而不是 1e-10。
        """
        c = self.calibs.get(ch or self.moving)
        if c is None:
            return np.asarray(xy, float).reshape(-1, 2)
        xy = np.asarray(xy, float).reshape(-1, 2)
        if c.poly_order <= 1:
            # order 1 是线性的: u = a0 + a1*nx + a2*ny, 解析求逆
            M = np.array([[c.poly_a[1], c.poly_a[2]],
                          [c.poly_b[1], c.poly_b[2]]], float)
            off = np.array([c.poly_a[0], c.poly_b[0]], float)
            n = (xy - off) @ np.linalg.inv(M)
            return np.asarray(c.center, float) + float(c.scale) * n
        p = xy.copy()                           # 高阶: 多项式接近恒等 -> xy 是好初值
        for _ in range(int(n_iter)):
            r = xy - c.apply_poly(p)
            if np.max(np.abs(r)) < tol:
                break
            J = c.jacobian_at(p)
            det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
            det = np.where(np.abs(det) < 1e-12, 1e-12, det)
            p = p + np.column_stack([
                (J[:, 1, 1] * r[:, 0] - J[:, 0, 1] * r[:, 1]) / det,
                (-J[:, 1, 0] * r[:, 0] + J[:, 0, 0] * r[:, 1]) / det])
        return p

    def apply_locs(self, locs, ch=None, frame='analysis'):
        """ch 的 Localizations -> 参考坐标系。

        ``frame='analysis'``: 落在参考**分析窗口**坐标(与裁切后图像一致);
        ``frame='calib'``   : 落在参考**标定**坐标。
        """
        from dataclasses import replace
        xy = np.column_stack([locs.x, locs.y])
        out = (self.to_reference(xy, ch) if frame == 'analysis'
               else self.to_reference_calib(xy, ch))
        return replace(locs, x=out[:, 0], y=out[:, 1])

    __call__ = to_reference

    # ---------- 公共可用区 ----------
    def valid_box(self, n=64):
        """所有通道**同时**覆盖的区域, 用参考通道**分析**坐标表示。"""
        x0, x1, y0, y1 = self.window_bounds(self.reference, 'analysis')
        box_x, box_y = (x0, x1), (y0, y1)
        for ch in self.windows:
            if ch == self.reference:
                continue
            my0, mx0, mh, mw = self.windows[ch]
            corners = _perimeter(my0, mx0, mh, mw, n=n)
            p = self.to_reference(self.detector_to_local(corners, ch), ch)
            box_x = (max(box_x[0], float(p[:, 0].min())),
                     min(box_x[1], float(p[:, 0].max())))
            box_y = (max(box_y[0], float(p[:, 1].min())),
                     min(box_y[1], float(p[:, 1].max())))
        return {'x': box_x, 'y': box_y}

    def in_valid(self, xy, ch=None):
        """判断**参考分析坐标**里的点是否落在公共可用区内。"""
        vb = self.valid_box()
        p = self.to_reference(xy, ch) if ch else np.asarray(xy, float).reshape(-1, 2)
        x, y = p.T
        return ((x >= vb['x'][0]) & (x <= vb['x'][1]) &
                (y >= vb['y'][0]) & (y <= vb['y'][1]))

    # ---------- 诊断 ----------
    def report(self):
        print(f'  reference = {self.reference!r}, channels = {list(self.windows)}, '
              f'配准到参考: {list(self.calibs) or "无"}')
        for ch, w in self.windows.items():
            tag = '  <- 参考' if ch == self.reference else ''
            print(f'  {ch}: 分析窗口 detector y{w[0]}..{w[0] + w[2]} '
                  f'x{w[1]}..{w[1] + w[3]}  ({w[2]}x{w[3]}){tag}')
            cw = (self.calib_windows or {}).get(ch)
            if cw and (w[0] != cw[0] or w[1] != cw[1]):
                d = self.calib_delta(ch)
                print(f'      标定窗口 y{cw[0]} x{cw[1]}  ->  '
                      f'分析坐标 + ({d[0]:+.0f}, {d[1]:+.0f}) = 标定坐标')
        for ch, cal in self.calibs.items():
            cov = cal.coverage()
            print(f'  calib[{ch}]: order {cal.poly_order}, '
                  f'coarse_offset {tuple(cal.coarse_offset)}, '
                  f'coverage x{cov.get("x")} y{cov.get("y")}')
        vb = self.valid_box()
        print(f"  公共可用区(参考分析坐标): x {vb['x'][0]:.1f}~{vb['x'][1]:.1f}, "
              f"y {vb['y'][0]:.1f}~{vb['y'][1]:.1f}")
        for m in self._warned:
            print(f'  [!] {m}')

    def __repr__(self):
        chs = ','.join(self.calibs) or '-'
        return (f'DualViewGeometry(ref={self.reference!r}, '
                f'channels={list(self.windows)}, calib={chs})')


# ---------------------------------------------------------------------------
#  通道间坐标映射 (单点工具, 供漂移 / trace 等复用)
# ---------------------------------------------------------------------------
def map_channel_xy(xy, geom, source, target):
    """Map ``[x, y]`` points between two channel analysis windows.

    Both channels must be the geometry reference or carry a registration to it;
    the mapping is ``to_reference`` followed by ``to_channel``.
    """
    points = np.asarray(xy, dtype=float).reshape(-1, 2)
    for channel in (source, target):
        if channel not in geom.windows:
            raise KeyError(f'unknown channel {channel!r}; '
                           f'available: {list(geom.windows)}')
    for channel in {source, target}:
        if channel != geom.reference and channel not in geom.calibs:
            raise ValueError(f'channel {channel!r} has no registration to '
                             f'geometry reference {geom.reference!r}')
    if source == target:
        return points.copy()
    return geom.to_channel(geom.to_reference(points, source), target)


# ---------------------------------------------------------------------------
#  工厂
# ---------------------------------------------------------------------------