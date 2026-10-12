# smlm/io/sif.py
from .base import Loader, RawDataset
from ._readsif import readsif
import numpy as np


class SIFLoader(Loader):
    def load(self, path, interval=1) -> RawDataset:
        info = readsif(path)
        # _readsif 返回的 imageData 为 (H, W, T) 列优先;统一转成 (T, H, W)
        imgs = np.moveaxis(info['imageData'], -1, 0)   # (H, W, T) -> (T, H, W)
        if interval > 1:
            imgs = imgs[::interval, :, :]              # 与 TIFLoader 一致:按帧轴抽取
        # order='C' 保证输出是 C 连续的 (T, H, W), 便于下游逐帧处理
        return RawDataset(
            data=np.ascontiguousarray(imgs, dtype=np.float32),
            exposure_time=info.get('exposureTime'),
            em_gain=info.get('gainDAC'),
            metadata=info,
        )