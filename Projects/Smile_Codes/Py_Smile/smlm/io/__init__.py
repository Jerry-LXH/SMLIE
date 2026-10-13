# smlm/io/__init__.py
from .base import RawDataset, Loader
from .sif import SIFLoader
from .tif import TIFLoader, save_tif
from .cache import run_cached, save_locs, load_locs, result_path

_LOADERS = {'sif': SIFLoader, 'tif': TIFLoader, 'tiff': TIFLoader}


def load_raw(path, interval=1) -> RawDataset:
    """按扩展名分派读取。"""
    ext = str(path).lower().rsplit('.', 1)[-1]
    loader = _LOADERS.get(ext)
    if loader is None:
        raise ValueError(f'Unsupported format: .{ext}') 
    return loader().load(path, interval)


__all__ = ['RawDataset', 'Loader', 'SIFLoader', 'TIFLoader', 'load_raw','save_tif','run_cached', 'save_locs', 'load_locs', 'result_path']