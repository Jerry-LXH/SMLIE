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

def save_tif(image, path, metadata=None, compress=None):
    """把 ImageStack / RawDataset / ndarray 存成 (T, H, W) 的 TIFF。

    - ImageStack / RawDataset: 取 ``.data``(帧在第一维), 并读出 ``pixel_size``;
    - ndarray: 直接写入(2D 自动补成单帧 (1, H, W));
    - pixel_size (nm/pixel) 会写入 ImageJ 分辨率标签(unit='um'), 便于下游读回。
    返回写出的路径。
    """
    if hasattr(image, 'data'):                  # ImageStack / RawDataset
        arr = np.asarray(image.data)
        px = getattr(image, 'pixel_size', None)
    else:                                        # 纯 ndarray
        arr = np.asarray(image)
        px = None
    if arr.ndim == 2:
        arr = arr[None, :, :]

    ij_meta = {'axes': 'TYX'}
    if metadata:
        ij_meta.update({k: v for k, v in metadata.items()
                        if isinstance(v, (str, int, float, bool))})

    kwargs = dict(imagej=True, metadata=ij_meta)
    if px:                                       # nm/pixel -> 像素/微米
        px_um = float(px) / 1000.0
        ij_meta['unit'] = 'um'
        kwargs['resolution'] = (1.0 / px_um, 1.0 / px_um)
    if compress is not None:
        kwargs['compression'] = compress

    tifffile.imwrite(path, arr, **kwargs)
    print(f"[io] saved tif {arr.shape} -> {path}")
    return path