from __future__ import annotations
import numpy as np
from .base import PSFModel


class CircularAperturePSF(PSFModel):
    """圆孔径傅立叶光学 PSF (光瞳 → |FT|² → PSF), 支持离焦(3D)。

    光瞳函数 P(ρ, θ) = A(ρ) · exp(i·2π·W(ρ, θ)):
      - A : 圆孔径 (归一化径向 ρ ≤ 1), 由 NA 决定截止频率 f_c = NA/λ;
      - W : 波前 (waves) = 静态 Zernike 像差 + 离焦项:
              W(ρ) = Σ_k z_k·Z_k(ρ,θ) + (NA²/2λ)·δz·ρ²,  δz = z - focus。

    实现路线 (逐 z 2D FFT 体积 + 插值, 与 GaussianPSF "核缓存 + 局部窗口"一致):
      1. [网格] 对一组 δz 平面, 每个平面在光瞳上叠加离焦相位后零填充 2D FFT,
         得到过采样强度 PSF 体积 (n_z, kh, kw), xy 间距 = kernel_spacing;
         体积只算一次并缓存。
      2. [插值] render 时按 emitter.xyz[:, 2] 在 z 上线性插值, 并在 xy 上做
         亚像素双线性插值 (同时完成到 fov.pixel_size 的重采样), 局部窗口截断叠加。

    单位: wavelength / xyz / pixel_size / focus 均为 nm; Zernike 系数统一存 'waves';
    z_range 为 emitter 的绝对 z 范围 (nm), 内部换算为相对 focus 的 δz。
    外部原地修改 emitters.xyz(漂移)后, 需 clear_cache() 重置。
    """

    def __init__(self, na: float, wavelength: float,
                 zernike=None, zernike_unit: str = 'waves',
                 pupil_size: int = 128, oversample: float = 2.0,
                 n_channels: int = 1, trunc: float | None = None,
                 focus: float = 0.0,
                 z_range: tuple[float, float] | None = None,
                 z_step: float | None = None):
        if not (0.0 < float(na) < 2.0):
            raise ValueError(f"NA 应在 (0, 2) 范围, 得到 {na}")
        if float(wavelength) <= 0:
            raise ValueError(f"wavelength 必须为正, 得到 {wavelength}")
        self.na = float(na)
        self.wavelength = float(wavelength)          # nm
        self.pupil_size = int(pupil_size)            # 光瞳网格尺寸(直径方向采样数)
        self.oversample = float(oversample)          # PSF 过采样倍数
        self.n_channels = int(n_channels)
        self.trunc = float(trunc) if trunc is not None else None
        self.focus = float(focus)                    # 焦平面 z (nm), δz = z - focus
        self.z_range = (float(z_range[0]), float(z_range[1])) if z_range else None
        self.z_step = float(z_step) if z_step is not None else None
        # Zernike 系数统一存为 "waves" 单位
        self.zernike = self._normalize_zernike(zernike)
        if zernike_unit == 'nm':
            self.zernike = {k: v / self.wavelength
                            for k, v in self.zernike.items()}
        elif zernike_unit != 'waves':
            raise ValueError(
                f"zernike_unit 只能是 'waves' 或 'nm', 得到 {zernike_unit!r}")
        self._cache = {}                             # PSF 体积/渲染缓存

    # ------------------------------------------------------------------
    # 由 NA / λ 决定的物理量
    # ------------------------------------------------------------------
    @property
    def f_cutoff(self) -> float:
        """截止频率 (cycles/nm) = NA / λ。"""
        return self.na / self.wavelength

    @property
    def nyquist_pixel(self) -> float:
        """衍射极限 PSF 的 Nyquist 采样间距 (nm) = λ / (2·NA)。"""
        return self.wavelength / (2.0 * self.na)

    @property
    def kernel_spacing(self) -> float:
        """过采样 PSF 核的网格间距 (nm) = λ / (2·NA·oversample)。"""
        return self.wavelength / (2.0 * self.na * self.oversample)

    @property
    def defocus_waves_per_nm(self) -> float:
        """离焦灵敏度: ρ² 项波前系数 (waves/nm) = NA² / (2λ)。

        离焦波前 W_defocus(ρ) = defocus_waves_per_nm · δz · ρ², δz = z - focus。
        """
        return self.na ** 2 / (2.0 * self.wavelength)

    # ------------------------------------------------------------------
    # Zernike 系数规范化: Noll 索引 1..15, 单位换算到 'waves'
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_zernike(zernike) -> dict:
        """把 Zernike 系数归一成 {Noll 索引: 系数} (保留输入单位, 限 1..15)。"""
        if zernike is None:
            return {}
        if isinstance(zernike, dict):
            coeffs = {int(k): float(v) for k, v in zernike.items() if float(v) != 0.0}
        else:                                       # 长度 ≤ 16 的数组, 索引即 Noll
            arr = np.asarray(zernike, dtype=float).ravel()
            coeffs = {i + 1: float(v) for i, v in enumerate(arr) if v != 0.0}
        for idx in coeffs:
            if not (1 <= idx <= 15):
                raise ValueError(
                    f"Zernike Noll 索引 {idx} 超出支持范围 1..15 (至多 15 阶)")
        return coeffs

    # ------------------------------------------------------------------
    # Zernike 多项式 (Noll 1..15, Noll 归一化)
    # ------------------------------------------------------------------
    @staticmethod
    def _zernike(noll: int, rho: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """返回归一化 Zernike 多项式 Z_j(ρ, θ) (Noll 索引 j, 1..15)。"""
        r, t, sq = rho, theta, np.sqrt
        if noll == 1:                                   # 平移 (piston)
            return np.ones_like(r)
        if noll == 2:                                   # x 倾斜 (tip)
            return 2.0 * r * np.cos(t)
        if noll == 3:                                   # y 倾斜 (tilt)
            return 2.0 * r * np.sin(t)
        if noll == 4:                                   # 离焦 (defocus)
            return sq(3.0) * (2.0 * r ** 2 - 1.0)
        if noll == 5:                                   # 斜像散
            return sq(6.0) * r ** 2 * np.sin(2.0 * t)
        if noll == 6:                                   # 垂直像散
            return sq(6.0) * r ** 2 * np.cos(2.0 * t)
        if noll == 7:                                   # 垂直彗差
            return sq(8.0) * (3.0 * r ** 3 - 2.0 * r) * np.sin(t)
        if noll == 8:                                   # 水平彗差
            return sq(8.0) * (3.0 * r ** 3 - 2.0 * r) * np.cos(t)
        if noll == 9:                                   # 垂直三叶
            return sq(8.0) * r ** 3 * np.sin(3.0 * t)
        if noll == 10:                                  # 斜三叶
            return sq(8.0) * r ** 3 * np.cos(3.0 * t)
        if noll == 11:                                  # 初级球差
            return sq(5.0) * (6.0 * r ** 4 - 6.0 * r ** 2 + 1.0)
        if noll == 12:
            return sq(10.0) * (4.0 * r ** 4 - 3.0 * r ** 2) * np.cos(2.0 * t)
        if noll == 13:
            return sq(10.0) * (4.0 * r ** 4 - 3.0 * r ** 2) * np.sin(2.0 * t)
        if noll == 14:
            return sq(10.0) * r ** 4 * np.cos(4.0 * t)
        if noll == 15:
            return sq(10.0) * r ** 4 * np.sin(4.0 * t)
        raise ValueError(f"不支持的 Noll 索引: {noll}")

    # ------------------------------------------------------------------
    # 光瞳 / z 网格 / PSF 体积
    # ------------------------------------------------------------------
    def _build_pupil(self, z_offset: float,
                     n_pad: int | None = None) -> np.ndarray:
        """构造轴向偏移 z_offset (δz, nm) 处的复光瞳函数。

        网格: 归一化频率坐标 (cycles/sample, 原点在角上, 与 np.fft.fft2 约定一致),
        孔径圆心在原点(周期性中心), 孔径半径 R = (pupil_size/2)/n_pad 个采样。
        返回 (n_pad, n_pad) 复光瞳, 可直接 fft2 -> fftshift 得到**居中**的 PSF
        (避免偶数网格半像素错位)。

        W(waves) = Σ z_k·Z_k(ρ,θ) + defocus_waves_per_nm·z_offset·ρ²,
        其中 ρ 为归一化径向 [0,1] (除以孔径半径 R)。
        """
        ps = self.pupil_size
        if n_pad is None:
            n_pad = max(int(round(ps * self.oversample)), ps)
        fx = np.fft.fftfreq(n_pad)                      # [-0.5, 0.5), 原点在 0
        X, Y = np.meshgrid(fx, fx)
        rho_cyc = np.hypot(X, Y)                        # 频率半径 (cycles/sample)
        R = (ps / 2.0) / n_pad                          # 孔径半径 (cycles/sample)
        A = (rho_cyc <= R).astype(float)
        rho = np.clip(rho_cyc / R, 0.0, 1.0)            # 归一化径向 [0,1] (Zernike 用)
        theta = np.arctan2(Y, X)

        W = np.zeros_like(rho)                          # 波前, 单位 waves
        for noll, coeff in self.zernike.items():
            W += coeff * self._zernike(noll, rho, theta)
        W += self.defocus_waves_per_nm * float(z_offset) * rho ** 2
        return A * np.exp(1j * 2.0 * np.pi * W)

    def _resolve_z_grid(self, emitters) -> np.ndarray:
        """确定 PSF 体积的 δz 网格 (nm, 相对 focus)。

        优先用 init 的 z_range/z_step (z_range 为绝对 z, 内部减 focus);
        未给 z_range 时从 emitters.xyz[:, 2] 范围加 1 步边距自动生成;
        未给 z_step 时默认取离焦引起 ~λ/4 (0.25 waves) 波前变化的 δz 间隔。
        """
        if self.z_step is not None:
            dz_step = float(self.z_step)
        else:
            # ρ=1 处波前变化 0.25 waves 对应的 δz
            dz_step = 0.25 / self.defocus_waves_per_nm

        if self.z_range is not None:
            lo = self.z_range[0] - self.focus
            hi = self.z_range[1] - self.focus
        else:
            z = np.asarray(emitters.xyz[:, 2], dtype=float)
            if z.size == 0:
                return np.array([0.0])
            lo = float(z.min()) - self.focus - dz_step
            hi = float(z.max()) - self.focus + dz_step

        if hi <= lo:                                    # 单平面 (全部同 z)
            return np.array([0.5 * (lo + hi)])
        n = int(np.ceil((hi - lo) / dz_step)) + 1
        return np.linspace(lo, hi, n)

    def _build_psf_volume(self, z_grid: np.ndarray) -> np.ndarray:
        """计算过采样强度 PSF 体积: 逐 z 平面 |FFT(pupil)|²。

        光瞳在 n_pad = round(pupil_size·oversample) 网格上按"角点原点"约定构建,
        fft2 后 fftshift 即得**居中**的 PSF (无偶数网格半像素错位), 取模平方,
        每平面归一化使核求和 = 1。返回 (n_z, n_pad, n_pad),
        xy 间距 = kernel_spacing。
        """
        ps = self.pupil_size
        n_pad = max(int(round(ps * self.oversample)), ps)
        vol = np.empty((len(z_grid), n_pad, n_pad), dtype=float)
        for i, dz in enumerate(z_grid):
            pupil = self._build_pupil(float(dz), n_pad)     # (n_pad, n_pad) 复数
            field = np.fft.fftshift(np.fft.fft2(pupil))
            psf = np.abs(field) ** 2
            s = psf.sum()
            if s > 0:
                psf /= s                                    # 单光子核求和 = 1
            vol[i] = psf
        return vol

    # ------------------------------------------------------------------
    # 调试/单测小工具: 不依赖任何渲染, 直接取单个分子某 z 处的 PSF 核
    # ------------------------------------------------------------------
    def psf_at(self, z: float) -> np.ndarray:
        """返回单个分子在绝对 z (nm) 处的归一化强度 PSF 核 (不含渲染)。

        等价于 _build_psf_volume 的单平面结果:
          - 形状 (n_pad, n_pad), n_pad = round(pupil_size·oversample);
          - 网格间距 = kernel_spacing (nm), 峰值在阵列中心;
          - 核求和 = 1 (单光子), 可直接作为通量权重。

        调参/单测时非常方便: z = focus 处应看到衍射极限 Airy 环
        (第一暗环 r ≈ 0.61·λ/NA), z 偏离 focus 时出现离焦形变/环状扩散。

        坐标轴(单位 nm)可用下式构造 (原点取 fftshift 的 DC 索引 n//2):
            n = kernel.shape[0]
            r_axis = (np.arange(n) - n // 2) * self.kernel_spacing
        """
        dz = float(z) - self.focus
        vol = self._build_psf_volume(np.array([dz]))
        return vol[0]

    # ------------------------------------------------------------------
    # 插值工具
    # ------------------------------------------------------------------
    @staticmethod
    def _interp_z(vol: np.ndarray, z_grid: np.ndarray, dz: float) -> np.ndarray:
        """在 δz 上对 PSF 体积做线性插值, 越界则夹到端面。"""
        n = len(z_grid)
        if n == 1 or dz <= z_grid[0]:
            return vol[0]
        if dz >= z_grid[-1]:
            return vol[-1]
        j = int(np.searchsorted(z_grid, dz)) - 1
        j = min(max(j, 0), n - 2)
        t = (dz - z_grid[j]) / (z_grid[j + 1] - z_grid[j])
        return (1.0 - t) * vol[j] + t * vol[j + 1]

    @staticmethod
    def _bilinear(K: np.ndarray, kc: np.ndarray,
                  qy: np.ndarray, qx: np.ndarray) -> np.ndarray:
        """在核 K (间距均匀, 坐标 kc) 上对查询网格 (qy × qx) 做双线性插值。

        返回 (len(qy), len(qx)); 落在核外的点取 0。
        """
        ny, nx = K.shape
        h = kc[1] - kc[0]
        fy = (qy - kc[0]) / h
        fx = (qx - kc[0]) / h
        iy = np.floor(fy).astype(int)
        ix = np.floor(fx).astype(int)
        wy = (fy - iy)[:, None]
        wx = (fx - ix)[None, :]
        vy = ((iy >= 0) & (iy < ny - 1))[:, None]
        vx = ((ix >= 0) & (ix < nx - 1))[None, :]
        iyc = np.clip(iy, 0, ny - 2)
        ixc = np.clip(ix, 0, nx - 2)
        K00 = K[iyc[:, None], ixc[None, :]]
        K01 = K[iyc[:, None], ixc[None, :] + 1]
        K10 = K[iyc[:, None] + 1, ixc[None, :]]
        K11 = K[iyc[:, None] + 1, ixc[None, :] + 1]
        out = (K00 * (1 - wy) * (1 - wx) + K01 * (1 - wy) * wx
               + K10 * wy * (1 - wx) + K11 * wy * wx)
        return np.where(vy & vx, out, 0.0)

    # ------------------------------------------------------------------
    # 渲染
    # ------------------------------------------------------------------
    def _build_render_kernels(self, emitters, fov, z_grid, vol):
        """为每个 emitter 预计算像面局部核及其像素切片 (对齐 GaussianPSF 接口)。

        对每个分子: 按其 z 插值出核体积切片, 再在像素中心处对 xy 做双线性插值
        (同时完成到 fov.pixel_size 的重采样与亚像素平移), 乘通量因子 (px/h)²。
        返回 (kernels, slices), slices[i] = (r0, r1, c0, c1)。
        """
        px = fov.pixel_size
        h = self.kernel_spacing
        n_pad = vol.shape[1]
        # 注意: fftshift 把 DC 放在索引 n_pad//2 (偶数网格), 原点取这里
        kc = (np.arange(n_pad) - n_pad // 2) * h            # 核样点偏移 (nm)
        half_ext = (n_pad // 2) * h                         # 核物理半幅 (nm)

        if self.trunc is not None:
            radius_nm = min(self.trunc * self.nyquist_pixel, half_ext)
        else:
            radius_nm = half_ext
        r_pix = max(int(np.ceil(radius_nm / px)), 2)
        flux = (px / h) ** 2                                # 间距 h→px 的通量补偿

        kernels, slices = [], []
        for i in range(emitters.n):
            x0, y0, z0 = (emitters.xyz[i, 0], emitters.xyz[i, 1],
                          emitters.xyz[i, 2])
            K = self._interp_z(vol, z_grid, z0 - self.focus)
            c = int(np.round(x0 / px - 0.5))                # 峰值所在像素列
            rw = int(np.round(y0 / px - 0.5))               # 峰值所在像素行
            c0, c1 = c - r_pix, c + r_pix + 1
            r0, r1 = rw - r_pix, rw + r_pix + 1
            dx = (np.arange(c0, c1) + 0.5) * px - x0         # 像素中心相对偏移 (nm)
            dy = (np.arange(r0, r1) + 0.5) * px - y0
            k = self._bilinear(K, kc, dy, dx) * flux
            kernels.append(k)
            slices.append((r0, r1, c0, c1))
        return kernels, slices

    def render(self, emitters, photon_frame, fov):
        """把该帧光子按 3D PSF 铺到像面, 返回 (C, H, W) 理想光子图。

        1. 解析 δz 网格并(按需)构建/复用缓存的过采样 PSF 体积;
        2. 逐 emitter 在 z 上插值 + xy 双线性插值得到像面局部核 (缓存于 xyz 身份);
        3. 在局部窗口内以 photon_frame 加权累加 (参考 GaussianPSF.render)。
        """
        img = np.zeros((self.n_channels, fov.height, fov.width))
        if emitters.n == 0:
            return img

        # 1. PSF 体积: 仅当 δz 网格变化时重算
        z_grid = self._resolve_z_grid(emitters)
        vc = self._cache.get('_volume')
        if (vc is None or vc[0].shape != z_grid.shape
                or not np.allclose(vc[0], z_grid)):
            vol = self._build_psf_volume(z_grid)
            self._cache['_volume'] = (z_grid, vol)
        else:
            z_grid, vol = vc

        # 2. 每分子局部核: 按 emitters.xyz 身份缓存 (与 GaussianPSF 一致)
        key = ('render', id(emitters))
        cached = self._cache.get(key)
        if cached is None or cached[0] is not emitters.xyz:
            kernels, slices = self._build_render_kernels(
                emitters, fov, z_grid, vol)
            cached = (emitters.xyz, kernels, slices)
            self._cache[key] = cached
        _, kernels, slices = cached

        # 3. 局部窗口加权累加
        H, W = fov.height, fov.width
        for i in range(emitters.n):
            r0, r1, c0, c1 = slices[i]
            R0, C0 = max(r0, 0), max(c0, 0)
            R1, C1 = min(r1, H), min(c1, W)
            if R1 <= R0 or C1 <= C0:
                continue
            k = kernels[i][R0 - r0:R1 - r0, C0 - c0:C1 - c0]
            for c in range(self.n_channels):
                img[c, R0:R1, C0:C1] += photon_frame[i, c] * k
        return img

    def clear_cache(self):
        self._cache.clear()