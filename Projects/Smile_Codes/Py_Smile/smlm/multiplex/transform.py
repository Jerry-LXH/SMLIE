# smlm/multiplex/transform.py
"""通道之间的 2D 配准模型 (mov -> ref), 以及它的落盘/载入。

坐标约定
--------
一切输入输出都是**各通道自己的 local 坐标**(窗口左上角为 0), 与 ImageStack 一致。
xy 的分量顺序是 (x, y) = (列, 行)。

粗偏移的处理
------------
双视角把同一个物点放到探测器的 ``f(p)`` 和 ``split_x + f(p)`` 两处。两个通道
**各自从自己的半幅起点裁切**, 于是::

    green_local = f(p) - 0                 = f(p)
    red_local   = (split_x + f(p)) - split_x = f(p)

**两通道 local 坐标本来就重合** —— 256 px 的大平移在裁切那一刻就被吃掉了。
所以多项式只需要拟合几十 nm 量级的残差形变, 天然接近恒等, 数值条件很好。
实测标定点的配对位移中位数是 (0.13, -1.43) px, 也直接印证了这一点。

``coarse_offset`` = (dy, dx) = "mov 窗口原点 - ref 窗口原点"(探测器坐标), 它
**只是元数据**:
  - 报告/诊断: 让人一眼看出两个窗口在探测器上的相对位置;
  - 传给 ``channels.cut_pair(coarse=...)``, 保证分析窗口与标定窗口同一布局;
  - 供 ``coords.DualViewGeometry`` 做布局一致性检查。
它**不进入** ``apply_poly`` 的数学 —— 进了反而会错 256 px。

向后兼容
--------
旧版 ``*_dualcal.npz`` 没有 ``coarse_offset`` 键, ``load()`` 时默认 (0, 0) ——
于是 ``mov + coarse == mov``, 变换与原版逐点完全相同(已验证 0.000e+00 px)。
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import json
import numpy as np

NPZ_VERSION = 2


@dataclass
class DualViewCalibration:
    """通道间 2D 配准: mov 通道 local 坐标 -> ref 通道 local 坐标。

    poly_* : 归一化坐标下的多项式系数  (u = F@poly_a, v = F@poly_b)
    affine : (3,2),  ``[x y 1] @ affine -> (u, v)``

    设计要点:
      - 多项式项表由本类的 ``_terms()`` 自己生成, 不依赖拟合侧的代码 ——
        单独载入 npz 就能用, 也不会因为两边定义不同而错位。
      - ``coarse_offset`` 只是**元数据**, 不参与数学(见模块 docstring)。
    """
    poly_a: np.ndarray
    poly_b: np.ndarray
    center: np.ndarray
    scale: float
    affine: np.ndarray
    poly_order: int = 3
    roi: dict = field(default_factory=dict)
    metrics: dict = field(default_factory=dict)
    coarse_offset: tuple = (0, 0)        # (dy, dx) 探测器坐标下的粗偏移
    reference: str = 'g'                 # 参考(不动)的通道名
    moving: str = 'r'                    # 被变换的通道名
    version: int = NPZ_VERSION

    # ---------- 项表 / 设计矩阵 ----------
    @staticmethod
    def _terms(order):
        """完整 order 阶二维多项式的项 (i = x 幂, j = y 幂), i + j <= order。

        顺序 = 总幂次递增、同阶按 x 幂递减; 低阶是高阶的前缀, 所以自动降阶时
        系数一一对应。**改了这里必须同步改 fit.poly_terms()**, 否则旧 npz 错位。
        """
        return [(i, d - i) for d in range(int(order) + 1) for i in range(d, -1, -1)]

    @property
    def coarse_xy(self):
        """粗偏移转成 (x, y) 顺序, 便于和坐标数组比较(仅诊断用)。"""
        dy, dx = self.coarse_offset
        return np.array([float(dx), float(dy)])

    def _feat(self, xy_shifted):
        n = (np.asarray(xy_shifted, float) - np.asarray(self.center, float)) \
            / float(self.scale)
        x, y = n[:, 0], n[:, 1]
        return np.column_stack([x ** i * y ** j
                                for i, j in self._terms(self.poly_order)])

    # ---------- 变换 ----------
    def apply_poly(self, xy):
        F = self._feat(np.asarray(xy, float).reshape(-1, 2))
        return np.column_stack([F @ self.poly_a, F @ self.poly_b])

    def apply_affine(self, xy):
        p = np.asarray(xy, float).reshape(-1, 2)
        return np.column_stack([p, np.ones(len(p))]) @ self.affine

    def transform(self, xy, method='poly'):
        if method == 'poly':
            return self.apply_poly(xy)
        if method == 'affine':
            return self.apply_affine(xy)
        raise ValueError("method 只支持 'poly' / 'affine'")

    def apply_locs(self, locs, method='poly'):
        """把 mov 通道的 Localizations 搬到 ref 通道坐标系(不改原对象)。"""
        xy = self.transform(np.column_stack([locs.x, locs.y]), method)
        return replace(locs, x=xy[:, 0], y=xy[:, 1])

    __call__ = transform

    # ---------- 局部分析 ----------
    def jacobian_at(self, xy, h=0.5):
        """在 xy 处的雅可比 J = d(u,v)/d(x,y), 返回 (M,2,2) 或 (2,2)。

        用途: 把一个通道的**位移**(如漂移)换算到另一个通道 —— 全局平移在
        两通道里的像元位移只差这个线性化矩阵。
        """
        p = np.asarray(xy, float).reshape(-1, 2)
        single = np.asarray(xy, float).ndim == 1
        out = np.empty((len(p), 2, 2))
        for k in range(2):
            e = np.zeros(2)
            e[k] = float(h)
            out[:, :, k] = (self.apply_poly(p + e) - self.apply_poly(p - e)) / (2.0 * h)
        return out[0] if single else out

    def coverage(self):
        """标定点覆盖的矩形范围 (来自 fit 时写入的 metrics, 为 mov 输入坐标)。"""
        return self.metrics.get('coverage', {})

    # ---------- 诊断 ----------
    def report(self):
        print(f'  {self.moving} -> {self.reference}   '
              f'(order={self.poly_order}, version={self.version})')
        print(f'  coarse_offset (dy,dx) = {tuple(self.coarse_offset)}  '
              f'(窗口在探测器上的位置差, 已由裁切吸收, 不参与变换)')
        print(f'  roi = {self.roi}')
        for k, v in self.metrics.items():
            print(f'  {k:20s} = {v}')

    @property
    def T33(self):
        """MATLAB 版 3x3 齐次矩阵(行向量右乘), 便于和旧代码逐项对照。"""
        A = np.hstack([np.asarray(self.affine, float), np.zeros((3, 1))])
        A[2, 2] = 1.0
        return A

    def __repr__(self):
        return (f'DualViewCalibration({self.moving}->{self.reference}, '
                f'order={self.poly_order}, coarse={tuple(self.coarse_offset)})')

    # ---------- 持久化 ----------
    def save(self, path):
        np.savez(path,
                 poly_a=self.poly_a, poly_b=self.poly_b,
                 center=self.center, scale=float(self.scale),
                 affine=self.affine, poly_order=int(self.poly_order),
                 coarse_offset=np.asarray(self.coarse_offset, float),
                 reference=str(self.reference), moving=str(self.moving),
                 version=int(self.version),
                 roi=json.dumps(self.roi, ensure_ascii=False),
                 metrics=json.dumps(self.metrics, ensure_ascii=False))
        return path

    @classmethod
    def load(cls, path):
        z = np.load(path, allow_pickle=False)
        kw = dict(poly_a=z['poly_a'], poly_b=z['poly_b'],
                  center=z['center'], scale=float(z['scale']),
                  affine=z['affine'], poly_order=int(z['poly_order']),
                  roi=json.loads(str(z['roi'])),
                  metrics=json.loads(str(z['metrics'])))
        # v2 才有的字段: 旧文件缺了也能读(此时 coarse=(0,0), 变换与原版逐点相同)
        if 'coarse_offset' in z:
            kw['coarse_offset'] = tuple(float(v) for v in z['coarse_offset'])
        if 'reference' in z:
            kw['reference'] = str(z['reference'])
        if 'moving' in z:
            kw['moving'] = str(z['moving'])
        kw['version'] = int(z['version']) if 'version' in z else 1
        return cls(**kw)
