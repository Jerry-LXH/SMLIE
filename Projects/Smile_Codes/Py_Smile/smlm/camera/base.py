# smlm/camera/base.py
from abc import ABC, abstractmethod
import numpy as np
from ..datatypes import ImageStack, RawDataset   


class Camera(ABC):
    @abstractmethod
    def adu_to_photons(self, images: np.ndarray) -> np.ndarray:
        """纯数值核心:ADU 数组 → 光子数组。"""
        ...

    def process(self, raw, pixel_size=None) -> ImageStack:
        """桥:RawDataset(ADU) → ImageStack(光子)。所有子类共用。"""
        photons = self.adu_to_photons(raw.data)
        # pixel_size 优先用显式传入,其次尝试从 metadata 读
        if pixel_size is None:
            pixel_size = raw.metadata.get('pixel_size')
        meta = dict(raw.metadata)
        meta['camera'] = type(self).__name__
        return ImageStack(data=photons.astype(np.float32),
                          pixel_size=pixel_size,
                          metadata=meta)


