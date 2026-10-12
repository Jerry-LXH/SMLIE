from __future__ import annotations
import numpy as np
from .base import PSFModel


class CalibratedPSF(PSFModel):
    """基于实测 PSF 数据的经验 PSF (spline 压缩表示)。

    适用场景: 已测到某个 z 扫描的 PSF 体积, 不想走"相位反演→Zernike"路径
    (那条路本质等价于 CircularAperturePSF, 直接传拟合出的 zernike 系数即可),
    而是以实测数据为真值。为减少数据量并支持任意亚像素/任意 z 求值,
    用**粗网格的三次张量积 B-spline** 表示整个体积。

    数据模型 (约定与 CircularAperturePSF 一致):
      data      : (n_z, n_y, n_x) 实测 PSF 强度体积, 单位光子; 内部每层归一化为
                  求和=1 的单光子核, 负值截断到 0;
      pixel_size: 实测像元尺寸 (nm/px);
      z_axis    : (n_z,) 各层轴向位置 (nm), 严格单调递增;
      focus     : 焦面 z (nm), render 时 δz = z - focus。

    Spline 表示 (核心, 数据压缩):
      PSF(z,y,x) = Σ_{i,j,k} c_{ijk} · B_i(z)·B_j(y)·B_k(x)
      在 (z, y, x) 上每隔 (z_control_step, control_step, control_step) 个点取
      控制网格 (n_cz, n_cy, n_cx), 系数数 << 原始 n_z·ny·nx (如 100×64×64 → 
      25×16×16 ≈ 64 倍压缩), 且可在任意 (z,y,x) 连续求值(亚像素、任意 z)。
      (实现可用 scipy.interpolate 的张量积样条 / RegularGridInterpolator(cubic))

    接口:
      CalibratedPSF(data, pixel_size, z_axis, ...)   # 从原始体积建 spline
      CalibratedPSF.from_file(path, ...)             # 从磁盘加载(体积或已存 spline)
      .save(path)                                    # 存压缩 spline 表示
      .psf_at(z, spacing=None)                       # 单分子在绝对 z 处的核(不依赖渲染)
      .render(emitters, photon_frame, fov)           # 与其它 PSF 一致
    """

    def __init__(self, data, pixel_size: float, z_axis=None, focus: float = 0.0,
                 n_channels: int = 1, spline_degree: int = 3,
                 control_step: int = 2, z_control_step: int = 1):
        data = np.asarray(data, dtype=float)
        if data.ndim != 3:
            raise ValueError(
                f"data 应为 (n_z, n_y, n_x) 的 3D 体积, 得到 shape {data.shape}")
        if float(pixel_size) <= 0:
            raise ValueError(f"pixel_size 必须为正, 得到 {pixel_size}")
        self.n_z, self.ny, self.nx = data.shape
        self.pixel_size = float(pixel_size)          # nm/px

        if z_axis is None:
            z_axis = np.arange(self.n_z)
        z_axis = np.asarray(z_axis, dtype=float)
        if z_axis.shape != (self.n_z,):
            raise ValueError(f"z_axis 长度应等于 {self.n_z}, 得到 {z_axis.shape}")
        if not np.all(np.diff(z_axis) > 0):
            raise ValueError("z_axis 必须严格单调递增")
        self.z_axis = z_axis                          # (n_z,) nm

        self.focus = float(focus)
        self.n_channels = int(n_channels)
        self.spline_degree = int(spline_degree)       # 默认 3 = 三次
        self.control_step = int(control_step)         # xy 控制点步长(像素)
        self.z_control_step = int(z_control_step)     # z 控制点步长(层)

        # 归一化: 每层求和=1 (单光子核), 负值截断到 0
        d = np.clip(data, 0.0, None)
        s = d.sum(axis=(1, 2), keepdims=True)
        self._raw = np.where(s > 0, d / np.where(s > 0, s, 1.0), 0.0)

        self._spline = None                           # 压缩表示 (拟合后填充)
        self._cache = {}

    # ------------------------------------------------------------------
    # Spline 拟合 / 求值
    # ------------------------------------------------------------------
    def _fit_spline(self):
        """构建粗网格三次张量积 B-spline 表示 (惰性, 首次求值/渲染时调用)。

        计划:
          1. 在 (z, y, x) 轴每隔 (z_control_step, control_step, control_step)
             个采样取一个控制节点, 节点数 n_cz/n_cy/n_cx 决定压缩率;
          2. 用三次张量积 B-spline 拟合 self._raw (最小二乘或插值),
             得到系数 c (n_cz, n_cy, n_cx);
          3. 存入 self._spline (含节点向量 + 系数 + 原始网格元数据)。
        数据量: (n_cz·n_cy·n_cx) << (n_z·ny·nx)。
        """
        raise NotImplementedError("CalibratedPSF._fit_spline 尚未实现")

    def _evaluate(self, z: float, dy: np.ndarray, dx: np.ndarray) -> np.ndarray:
        """在相对焦面偏移 dy/dx (nm, 查询网格) 处求值 z 平面的 PSF 核。

        计划:
          - 查询坐标 (nm) → 样条坐标 (px): (offset - y_center)/pixel_size;
          - 用三次张量积样条在 (z, y, x) 同时求值 (向量化);
          - z 越界夹到 z_axis 端面; 负值截断到 0。
        返回 (len(dy), len(dx))。
        """
        raise NotImplementedError("CalibratedPSF._evaluate 尚未实现")

    # ------------------------------------------------------------------
    # 调试/单测小工具 (与 CircularAperturePSF.psf_at 对齐)
    # ------------------------------------------------------------------
    def psf_at(self, z: float, spacing: float | None = None) -> np.ndarray:
        """返回单个分子在绝对 z (nm) 处的归一化强度 PSF 核 (不含渲染)。

        spacing: 输出核网格间距 (nm); None = 实测 pixel_size。
        形状 (n_y, n_x), 求和 ≈ 1, 峰值在中心。
        坐标轴: (np.arange(n) - n // 2) * spacing。
        """
        raise NotImplementedError("CalibratedPSF.psf_at 尚未实现")

    # ------------------------------------------------------------------
    # 磁盘 IO
    # ------------------------------------------------------------------
    @classmethod
    def from_file(cls, path, spline_only: bool = False, **kwargs):
        """从磁盘加载:
          - spline_only=False: 从原始体积文件 (npz/tif, 含 z_axis/pixel_size) 载入,
            内部再建样条;
          - spline_only=True : 直接载入 save() 存出的压缩 spline 表示。
        kwargs 透传给 __init__ (如 focus/control_step 等)。
        """
        raise NotImplementedError("CalibratedPSF.from_file 尚未实现")

    def save(self, path):
        """把压缩 spline 表示存盘: 系数 + 节点/网格 + 元数据 (z_axis, pixel_size, ...)。"""
        raise NotImplementedError("CalibratedPSF.save 尚未实现")

    # ------------------------------------------------------------------
    # 渲染 (与其它 PSF 一致)
    # ------------------------------------------------------------------
    def render(self, emitters, photon_frame, fov):
        """把该帧光子按实测 PSF 铺到像面, 返回 (C, H, W) 理想光子图。

        计划: 对每个 emitter 按其 xyz[:,2] 在 spline 上**直接**求值 FOV 像素
        中心偏移处的核 (无需像 circular 的 z-grid 体积 + z 插值), 在局部窗口
        内加权累加; 核按 emitter 身份缓存。
        """
        raise NotImplementedError("CalibratedPSF.render 尚未实现")

    def clear_cache(self):
        self._cache.clear()
