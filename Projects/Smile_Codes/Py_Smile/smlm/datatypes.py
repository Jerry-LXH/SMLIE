# smlm/types.py
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

@dataclass
class RawDataset:
    """所有 loader 的统一输出:原始 ADU 图像 + metadata。

    约定与 ImageStack 一致:data 为 (T, H, W), 帧在**第一维**。
    """
    data: np.ndarray                 # (T, H, W) 原始 ADU 计数, float32
    exposure_time: float = None
    em_gain: float = None
    metadata: dict = field(default_factory=dict)

    @property
    def n_frames(self):
        return self.data.shape[0]


@dataclass
class ImageStack:
    """光子单位的图像栈,携带下游标定信息。camera 层的输出建议包成它。"""
    data: np.ndarray                 # (T, H, W) float32, 光子数
    pixel_size: float = None         # nm/pixel, 定位转 nm 时要用
    metadata: dict = field(default_factory=dict)

    @property
    def n_frames(self): return self.data.shape[0]
    @property
    def shape(self): return self.data.shape


@dataclass
class MultiImageStack:
    """多通道图像栈:各通道已裁到**同一物理窗口**,共用参考通道坐标系。

    为什么不直接给 ImageStack 加一个通道维:
        ``localize``(``extract_rois`` 按 ``data[f, ...]`` 取帧)、``imgprocess``
        (Crop / Pipeline / DriftCorrector 全按 ``(T,H,W)``)、``viz`` 都硬编码了
        "(T,H,W), 第 0 维是帧"。加一维会让 ``n_frames`` 静默变成通道数, 属于
        最坏的一类破坏。所以单独一个类。

    与单通道流程的衔接:
        ``channel(ch)`` 返回 ``data[i]`` 的**视图**包成的 ImageStack,
        因此 detect / localize / Crop / DriftCorrector / viz 全部可直接复用,
        定位输出仍是标准 ``Localizations``。

    坐标约定:
        - ``data`` : (C, T, H, W), 第 0 维是**通道**。
        - ``offsets[c] = (row, col)``: 该通道窗口左上角在**探测器坐标**里的位置。
          窗口是按参考通道定义的同一物理区域, 所以两通道 offsets 之差就是
          通道间的粗偏移(≈ split_x)。
        - 各通道自己的 local 坐标仍从 0 开始, 与 ImageStack 完全一致。
        - 通道之间的**精细**对应关系存在 ``reg``: ``{ch: DualViewCalibration}``,
          参考通道是恒等映射, 因此不在字典里。
    """
    data: np.ndarray                                  # (C, T, H, W) 光子数
    channels: list[str] = field(default_factory=lambda: ['g', 'r'])
    reference: int = 0                                # 参考通道索引, 恒等
    reg: dict = field(default_factory=dict)           # {ch: 变换对象}, 参考通道无
    offsets: np.ndarray = None                        # (C,2) 窗口左上角 (row, col)
    pixel_size: float = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.data = np.asarray(self.data)
        if self.data.ndim != 4:
            raise ValueError('MultiImageStack.data 必须是 (C,T,H,W), '
                             f'得到 {self.data.shape}')
        if len(self.channels) != self.data.shape[0]:
            raise ValueError(f'通道名 {len(self.channels)} 个, 但 data 的 C = '
                             f'{self.data.shape[0]}')
        if not 0 <= int(self.reference) < len(self.channels):
            raise ValueError(f'reference = {self.reference} 越界')
        self.reference = int(self.reference)
        if self.offsets is not None:
            self.offsets = np.asarray(self.offsets, dtype=int).reshape(-1, 2)
            if self.offsets.shape[0] != len(self.channels):
                raise ValueError(f'offsets 需要 {len(self.channels)} 行, '
                                 f'得到 {self.offsets.shape[0]}')

    # ---------- 基本属性 ----------
    @property
    def n_channels(self): return len(self.channels)

    @property
    def n_frames(self): return self.data.shape[1]

    @property
    def shape(self): return self.data.shape           # (C, T, H, W)

    @property
    def ref_name(self): return self.channels[self.reference]

    # ---------- 通道访问 ----------
    def index(self, ch) -> int:
        """通道名(或索引) -> 通道索引。"""
        if isinstance(ch, (int, np.integer)):
            i = int(ch)
            if not 0 <= i < self.n_channels:
                raise IndexError(f'通道索引 {i} 越界')
            return i
        try:
            return list(self.channels).index(ch)
        except ValueError:
            raise KeyError(f'未知通道 {ch!r}, 可选 {self.channels}') from None

    def channel(self, ch) -> ImageStack:
        """取单个通道, 返回 ImageStack **视图**(不拷贝)。

        metadata 里补上通道名与该通道窗口在探测器上的偏移, 与 ``Crop``
        写在 ``crop_offset`` 里的语义一致, 下游坐标换算可直接复用。
        """
        i = self.index(ch)
        off = ((0, 0) if self.offsets is None
               else tuple(int(v) for v in self.offsets[i]))
        meta = dict(self.metadata)
        meta.update({'channel': self.channels[i],
                     'multi': True,
                     'crop_offset': off})
        return ImageStack(data=self.data[i], pixel_size=self.pixel_size,
                          metadata=meta)

    def channel_offset(self, ch) -> tuple[int, int]:
        """该通道窗口左上角在探测器坐标里的 (row, col)。"""
        if self.offsets is None:
            return (0, 0)
        return tuple(int(v) for v in self.offsets[self.index(ch)])

    # ---------- 坐标层: 通道 local -> 参考通道 local ----------
    def geometry(self, margin=None, check_coverage=True, strict=False):
        """构造 ``DualViewGeometry``(分析窗口 = offsets, 配准 = ``reg``)。

        这是本类与坐标层的**唯一**接口: 窗口原点、各通道配准、标定 roi 都汇总在
        这里, 下游(漂移估计 / 定位搬运 / 叠图)因此不用再手传一堆偏移量。
        """
        from .multiplex.coords import DualViewGeometry
        _, _, h, w = self.data.shape
        offs = (np.zeros((len(self.channels), 2), int) if self.offsets is None
                else np.asarray(self.offsets, int))
        windows = {ch: (int(offs[i, 0]), int(offs[i, 1]), int(h), int(w))
                   for i, ch in enumerate(self.channels)}
        roi_d = (self.metadata or {}).get('roi')
        roi_d = roi_d if isinstance(roi_d, dict) else None
        return DualViewGeometry(
            windows=windows, calib=dict(self.reg or {}),
            reference=self.ref_name, roi=roi_d,
            margin=int(margin if margin is not None
                       else (roi_d or {}).get('margin', 0)),
            check_coverage=check_coverage, strict=strict)

    def to_reference(self, xy, ch) -> np.ndarray:
        """该通道 local 坐标 -> **参考通道** local 坐标(含窗口与标定偏移)。

        注意这里只做**坐标**变换, 不动像素 —— 亚像素配准由坐标承载, 图像只做
        整体平移(漂移校正)。
        """
        if ch == self.ref_name:
            return np.asarray(xy, float).reshape(-1, 2).copy()
        return self.geometry(check_coverage=False).to_reference(xy, ch)

    # ---------- 派生 ----------
    def subset(self, channels) -> 'MultiImageStack':
        """取若干通道(保持给定顺序), reference 随之重新定位。"""
        idx = [self.index(c) for c in channels]
        names = [self.channels[i] for i in idx]
        if self.ref_name not in names:
            raise ValueError(f'通道子集必须包含参考通道 {self.ref_name!r}')
        offs = None if self.offsets is None else np.asarray(self.offsets)[idx]
        return MultiImageStack(
            data=self.data[idx], channels=names,
            reference=names.index(self.ref_name),
            reg={c: v for c, v in self.reg.items() if c in names},
            offsets=offs, pixel_size=self.pixel_size,
            metadata=dict(self.metadata))

    def with_data(self, data, offsets=None,
                  metadata=None) -> 'MultiImageStack':
        """换掉像素(例如漂移校正后), 通道/配准/元数据保持一致。"""
        return MultiImageStack(
            data=data, channels=list(self.channels), reference=self.reference,
            reg=dict(self.reg),
            offsets=self.offsets if offsets is None else offsets,
            pixel_size=self.pixel_size,
            metadata=dict(self.metadata) if metadata is None else dict(metadata))

    def __repr__(self):
        c, t, h, w = self.data.shape
        return (f'MultiImageStack(C={c}, T={t}, H={h}, W={w}, '
                f'channels={self.channels}, ref={self.ref_name!r})')


@dataclass
class Candidates:
    """检测输出:整数像素级候选点。"""
    frame: np.ndarray                # int
    y: np.ndarray                    # int, 像素
    x: np.ndarray                    # int

    def __len__(self): return len(self.x)


@dataclass
class Localizations:
    """定位输出:列式存储,每个数组长度 = 分子总数。SMLM 的通用做法。

    3D 约定:z / z_lp 为轴向坐标与轴向定位精度, 单位 nm(由 z 标定给出,
    如像散 σx/σy→z、双螺旋角度→z)。2D 定位时保持 None, 与现有下游完全兼容。
    """
    frame: np.ndarray                # int
    x: np.ndarray                    # 亚像素坐标 (pixel)
    y: np.ndarray
    photons: np.ndarray              # 单分子光子数
    bg: np.ndarray                   # 局部背景 (光子/pixel)
    sigma: np.ndarray = None         # PSF 宽度 (2D 各向同性; 各向异性方案放 extra)
    lp: np.ndarray = None            # 横向定位精度 (nm)
    z: np.ndarray = None             # 轴向位置 (nm), 3D 定位才有; 2D 为 None
    z_lp: np.ndarray = None          # 轴向定位精度 (nm), 3D 才有; 2D 为 None
    extra: dict = field(default_factory=dict)   # 拟合残差、CRLB、sigma_x/sigma_y 等

    def __len__(self): return len(self.x)

    def __getitem__(self, idx):
        """支持布尔掩码 / 整数数组 / 切片, 同步切所有列(含 extra 里的等长数组)。"""
        def take(v):
            return None if v is None else np.asarray(v)[idx]

        new_extra = {}
        for k, v in (self.extra or {}).items():
            arr = np.asarray(v)
            # 只对与分子数等长的数组做切片, 标量/CRLB 常数等原样保留
            if arr.ndim >= 1 and arr.shape[0] == len(self):
                new_extra[k] = arr[idx]
            else:
                new_extra[k] = v

        return Localizations(
            frame=take(self.frame), x=take(self.x), y=take(self.y),
            photons=take(self.photons), bg=take(self.bg),
            sigma=take(self.sigma), lp=take(self.lp),
            z=take(self.z), z_lp=take(self.z_lp),
            extra=new_extra,
        )

    def to_dataframe(self):
        import pandas as pd
        cols = {'frame': self.frame, 'x': self.x, 'y': self.y,
                'photons': self.photons, 'bg': self.bg}
        if self.sigma is not None: cols['sigma'] = self.sigma
        if self.lp is not None:    cols['lp'] = self.lp
        if self.z is not None:     cols['z'] = self.z
        if self.z_lp is not None:  cols['z_lp'] = self.z_lp
        return pd.DataFrame(cols)

@dataclass
class DriftTrace:
    """逐帧累积漂移(drift estimate 与 correct 之间的契约)。

    约定(关键):
      - dx / dy 长度均为 n_frames, 单位 pixel, 索引为 0-based 帧号。
      - dx 对应列方向(x / col), dy 对应行方向(y / row)。
      - 校正定义: 帧 t 的坐标减去 (dx[t], dy[t]) 即为校正后坐标;
        图像按 (-dy[t], -dx[t]) 平移以抵消漂移。
    """
    dx: np.ndarray                   # (n_frames,) 列(x)方向累积漂移, pixel
    dy: np.ndarray                   # (n_frames,) 行(y)方向累积漂移, pixel
    pixel_size: float = None         # nm/pixel, 便于换算成 nm
    metadata: dict = field(default_factory=dict)

    def __len__(self): return len(self.dx)

    @property
    def n_frames(self): return len(self.dx)

    @property
    def total_drift(self):
        """末帧相对首帧的总漂移模长 (pixel)。"""
        return float(np.hypot(self.dx[-1] - self.dx[0],
                              self.dy[-1] - self.dy[0]))

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame({'frame': np.arange(len(self.dx)),
                             'dx': self.dx, 'dy': self.dy})

    @classmethod
    def zeros(cls, n_frames, pixel_size=None, **metadata):
        """创建一条零漂移轨迹 (常用于测试基线 / 无漂移占位)。"""
        return cls(dx=np.zeros(n_frames, dtype=float),
                dy=np.zeros(n_frames, dtype=float),
                pixel_size=pixel_size,
                metadata={'kind': 'zero', **metadata})
    

@dataclass
class Emitter:
    """单个发光体的连接轨迹。photons/bg/sigma 等按 loc_idx 回查源 Localizations,
    避免冗余存储。闪烁/漂白/常亮统一用派生属性表达。"""
    frame:   np.ndarray          # 出现的帧 (升序)
    x:       np.ndarray
    y:       np.ndarray
    loc_idx: np.ndarray          # 指向源 Localizations 的行索引
    bleached: bool = False       # True=在片内观测到漂白; False=片尾仍常亮

    # ---- 时间属性 ----
    @property
    def on_frame(self) -> int:   return int(self.frame[0])
    @property
    def last_frame(self) -> int: return int(self.frame[-1])
    @property
    def n_locs(self) -> int:     return int(self.frame.size)
    @property
    def span(self) -> int:                      # on..last 闭区间长度
        return self.last_frame - self.on_frame + 1
    @property
    def n_dark(self) -> int:                     # 闪烁的暗帧数
        return self.span - self.n_locs
    @property
    def duty_cycle(self) -> float:               # 占空比: 常亮≈1, 闪烁<1
        return self.n_locs / self.span

    # ---- 空间属性 ----
    @property
    def xy(self) -> np.ndarray:
        return np.column_stack([self.x, self.y])
    @property
    def centroid(self) -> np.ndarray:
        return np.array([self.x.mean(), self.y.mean()])

    def covariance(self) -> np.ndarray:
        p = self.xy - self.centroid
        return (p.T @ p) / max(self.n_locs - 1, 1)

    def eccentricity(self) -> float:
        w = np.linalg.eigvalsh(self.covariance())   # 升序
        return float(w[-1] / max(w[0], 1e-12))

@dataclass
class EmitterStats:
    centroid: np.ndarray          # [N x 2]
    survival_sec: np.ndarray
    n_locs: np.ndarray
    duty_cycle: np.ndarray
    bleached: np.ndarray
    eccentricity: np.ndarray
    photons_per_sec: np.ndarray
    photons_sum: np.ndarray
    sigma_mean: np.ndarray
    bg_mean: np.ndarray


@dataclass
class Trace:
    emitter_id: int
    frames: np.ndarray            # 连续帧轴 (稠密, 含 OFF 帧)
    intensity: np.ndarray         # 分析主对象: 扣背景后信号 (或 raw)
    raw: np.ndarray               # 原始积分强度 (raw_series)
    bg: np.ndarray | None = None  # 每像素背景 (bg_series)
    snr: np.ndarray | None = None # 逐帧空间 SNR (snr_series)
    x: float = 0.0                # 提取中心 (像素坐标)
    y: float = 0.0
    meta: dict = field(default_factory=dict)  # edge / bg_edge / method / n_bg_pixels ...

    def __len__(self):
        return self.frames.size

@dataclass
class StateInfo:
    label: int
    mean_intensity: float
    std_intensity: float
    occupancy: float                     # 占空比 (相对 active 段)

@dataclass
class StateResult:
    emitter_id: int
    trace: np.ndarray                    # 被分析的信号 (通常 intensity)
    sequence: np.ndarray                 # 逐帧状态标签, 0 = OFF/背景
    transitions: np.ndarray              # (n,2): [frame_index, new_state]
    n_states: int
    state_info: list[StateInfo]
    fit_trace: np.ndarray                # 阶梯拟合
    background: float
    trace_bgsub: np.ndarray
    fit_trace_bgsub: np.ndarray
    bleach_frame: int | None             # 索引 (非 frame 值); None = 未漂白
    start_frame: int
    lifetime: int | None
    method: str
    meta: dict = field(default_factory=dict)