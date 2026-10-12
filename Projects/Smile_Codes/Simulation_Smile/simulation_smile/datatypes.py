"""核心数据类型与全局约定 (原 core.py)。

单位约定：长度 nm，时间 s，强度/亮度以“光子数”为单位。
随机数约定：一律使用外部传入的 np.random.Generator，模块内部绝不 seed。
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class EmitterSet:
    """发射体的静态/内禀属性。所有量都是 (N,) 或 (N, k) 的 numpy 数组。"""
    xyz: np.ndarray                       # (N, 3) 位置，单位 nm
    brightness: np.ndarray                # (N,) 内禀亮度（单位激发下的发射率，光子/s）
    channel: np.ndarray | None = None     # (N,) 光谱通道索引，双通道/FRET 预留
    orientation: np.ndarray | None = None # (N, 3) 偶极取向，暂留空
    size: np.ndarray | None = None        # (N,) 纳米粒子物理尺寸 nm，暂留空
    ids: np.ndarray | None = None         # (N,) 用于 ground truth 追踪

    def __post_init__(self):
        self.xyz = np.asarray(self.xyz, dtype=float)
        self.brightness = np.asarray(self.brightness, dtype=float)
        if self.channel is None:
            self.channel = np.zeros(self.n, dtype=int)
        if self.ids is None:
            self.ids = np.arange(self.n)

    @property
    def n(self) -> int:
        return self.xyz.shape[0]
    
    @property
    def config(self) -> dict:
        """自动收集所有非下划线开头的实例属性,做 JSON 友好转换。"""
        out = {}
        for k, v in vars(self).items():
            if k.startswith("_"):
                continue
            if isinstance(v, np.ndarray):
                v = v.tolist()
            out[k] = v
        return out


@dataclass
class PhotonStream:
    """光物理模拟输出：每个 emitter、每帧、每通道发射的光子数。"""
    counts: np.ndarray                    # (N, T, C)
    states: np.ndarray | None = None      # (N, T) 亮/暗/漂白状态，供 ML 标签


@dataclass
class SimulationResult:
    """整条管线的最终产物。ML 用途下 ground truth 与图像同等重要。"""
    images: np.ndarray                    # (T, C, H, W) 数字图像
    ground_truth: EmitterSet
    photon_stream: PhotonStream
    metadata: dict = field(default_factory=dict)


@dataclass
class FieldOfView:
    """视场定义：像素尺寸与像面网格。"""
    height: int                           # 像素行数
    width: int                            # 像素列数
    pixel_size: float                     # 每像素对应的物理尺寸 nm

    @property
    def shape(self) -> tuple[int, int]:
        return (self.height, self.width)

    def pixel_centers(self) -> tuple[np.ndarray, np.ndarray]:
        """返回像素中心的物理坐标网格 (X, Y)，单位 nm。原点在左上角。"""
        xs = (np.arange(self.width) + 0.5) * self.pixel_size
        ys = (np.arange(self.height) + 0.5) * self.pixel_size
        X, Y = np.meshgrid(xs, ys)        # (H, W)
        return X, Y
