from __future__ import annotations
import numpy as np
from .datatypes import SimulationResult, FieldOfView
from tqdm import tqdm

def _describe(obj):
    """记录组件类型;若组件自己暴露了 config 就一并存下。"""
    info = {"class": type(obj).__qualname__}
    if hasattr(obj, "config"):
        info["config"] = obj.config
    return info

class SimulationPipeline:
    def __init__(self, source, illumination, photophysics, psf, detector,
                 fov: FieldOfView, seed: int = 0):
        self.source = source
        self.illumination = illumination
        self.photophysics = photophysics
        self.psf = psf
        self.detector = detector
        self.fov = fov
        self.seed = seed

    def run(self, n_frames: int = 1, dt: float = 0.03) -> SimulationResult: # dt即曝光时间
        rng = np.random.default_rng(self.seed)

        emitters = self.source.sample(rng)                       # 分布
        stream = self.photophysics.simulate(                     # 光物理 + shot
            emitters, self.illumination, n_frames, dt, rng)      # 这一步已经以frames作为单位了

        images = []
        # 使用 tqdm 包装 range(n_frames)，并添加描述文字
        for t in tqdm(range(n_frames), desc="Generating frames", unit="frame"):
            photon_frame = stream.counts[:, t, :]                # (N, C)
            ideal = self.psf.render(emitters, photon_frame,      # 空间铺开(多项分配)
                                    self.fov, rng)
            digital = self.detector.capture(ideal, rng)          # 相机
            images.append(digital)
        images = np.stack(images)                                # (T, C, H, W)

        return SimulationResult(
            images=images,
            ground_truth=emitters,
            photon_stream=stream,
            metadata={
                "seed": self.seed,
                "dt": dt,
                "pixel_size": self.fov.pixel_size,
                "components": {
                    "source":       _describe(self.source),
                    "illumination": _describe(self.illumination),
                    "photophysics": _describe(self.photophysics),
                    "psf":          _describe(self.psf),
                    "detector":     _describe(self.detector),
                },
            },
        )