# smlm/io/tif.py
import tifffile
import numpy as np
from .base import Loader, RawDataset


class TIFLoader(Loader):
    def load(self, path, interval=1) -> RawDataset:
        arr = tifffile.imread(path)          # (T, H, W)
        if arr.ndim == 2:
            arr = arr[None, :, :]            # 单帧也统一成 (1, H, W)
        if interval > 1:
            arr = arr[::interval, :, :]      # 按帧轴抽取
        meta = {}
        with tifffile.TiffFile(path) as tf:
            if tf.imagej_metadata:
                meta.update(tf.imagej_metadata)
        return RawDataset(data=np.ascontiguousarray(arr, dtype=np.float32), metadata=meta)