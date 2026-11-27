# smlm/io/cache.py
import os
import json
import hashlib
import zipfile
import numpy as np
from ..datatypes import Localizations

# Localizations 的标准字段(顺序无所谓; 3D 定位的 z/z_lp 一并缓存)
_ARRAY_FIELDS = ['frame', 'x', 'y', 'photons', 'bg', 'sigma', 'lp', 'z', 'z_lp']


# ---------------------------------------------------------------------------
#  路径: 与源文件同目录同名, 换后缀
# ---------------------------------------------------------------------------
def result_path(src_path, suffix='.locs.npz'):
    base, _ = os.path.splitext(src_path)
    return base + suffix

def data_fingerprint(stack, n_probe=32):
    """对输入 ImageStack 算内容指纹: 覆盖一切上游改动(preprocess/camera)。

    只对 ``n_probe`` 帧(等间隔, 含首末帧)做全量哈希, 不全量哈希整个栈。
    原因: 全量要把跨步的裁切视图先 copy 成连续数组(2000x312x228 float64 = 1.1 GB),
    在 notebook 里实测单次 110 s —— 缓存**命中**也要付这个钱, 完全得不偿失。
    抽样仍能发现 gain / crop / 帧窗 / 相机参数 这类上游改动。
    """
    arr = np.asarray(stack.data)
    h = hashlib.blake2b(digest_size=16)
    h.update(str(arr.shape).encode())
    h.update(str(arr.dtype).encode())
    n = int(arr.shape[0]) if arr.ndim >= 1 else 0
    if n > 0:
        k = int(min(max(1, int(n_probe)), n))
        idx = np.unique(np.linspace(0, n - 1, k).round().astype(int))
        h.update(np.ascontiguousarray(arr[idx]).tobytes())
    else:
        h.update(np.ascontiguousarray(arr).tobytes())
    # 可读的上游元信息也一并纳入(便于排查, 也防纯 metadata 改动)
    h.update(json.dumps(describe(stack.metadata),
                        sort_keys=True, default=str).encode())
    return h.hexdigest()


def image_signature(stack):
    """Return a lightweight signature for a 3D image stack.

    The signature covers the full shape and dtype plus the first frame's
    content. It is intended for fast cache validation, not as a full-stack
    content fingerprint.
    """
    arr = np.asarray(stack.data)
    if arr.ndim != 3:
        raise ValueError(f'expected a 3D image stack, got shape {arr.shape}')
    if arr.shape[0] == 0:
        raise ValueError('cannot validate an empty image stack')

    first_frame = np.ascontiguousarray(arr[0])
    digest = hashlib.blake2b(digest_size=16)
    digest.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
    digest.update(arr.dtype.str.encode())
    digest.update(first_frame.tobytes())
    return {
        'shape': list(arr.shape),
        'dtype': arr.dtype.str,
        'first_frame_hash': digest.hexdigest(),
    }


# ---------------------------------------------------------------------------
#  把 method / detector / localizer 递归转成可 JSON 化的配置字典
# ---------------------------------------------------------------------------
def describe(obj, _depth=0):
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (list, tuple)):
        return [describe(v, _depth + 1) for v in obj]
    if isinstance(obj, dict):
        return {str(k): describe(v, _depth + 1) for k, v in obj.items()}
    if callable(obj):                              # 函数 / lambda
        return getattr(obj, '__name__', repr(obj))
    if _depth > 6:
        return repr(obj)
    if hasattr(obj, '__dict__'):                   # 普通对象: 抽公开属性
        d = {'__class__': type(obj).__name__}
        for k, v in vars(obj).items():
            if k.startswith('_'):
                continue
            d[k] = describe(v, _depth + 1)
        return d
    return repr(obj)


def config_hash(config):
    s = json.dumps(config, sort_keys=True, default=str)
    return hashlib.sha1(s.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
#  保存 / 加载 Localizations
# ---------------------------------------------------------------------------
def save_locs(locs, src_path, method=None, data_hash=None, extra_meta=None,
              compressed=True):
    path = result_path(src_path)
    arrays = {}
    meta = {'format_version': 1,
            'fields': [], 'extra_array': [], 'extra_scalar': {}}

    for name in _ARRAY_FIELDS:
        v = getattr(locs, name, None)
        if v is not None:
            arrays[name] = np.asarray(v)
            meta['fields'].append(name)

    # extra: 数组存进 npz, 标量存进 meta
    for k, v in (locs.extra or {}).items():
        if isinstance(v, np.ndarray) or (hasattr(v, '__len__')
                                         and not isinstance(v, (str, dict))):
            arrays['extra__' + k] = np.asarray(v)
            meta['extra_array'].append(k)
        else:
            meta['extra_scalar'][k] = v

    if method is not None:
        meta['config'] = describe(method)
        meta['config_hash'] = config_hash(meta['config'])

    if data_hash is not None:
        meta['data_hash'] = data_hash 
        
    if extra_meta:
        meta['meta'] = extra_meta

    arrays['__meta__'] = np.array(json.dumps(meta, default=str))
    if compressed:
        np.savez_compressed(path, **arrays)
    else:
        np.savez(path, **arrays)
    print(f"[cache] saved {len(locs)} locs -> {path}")
    return path


def load_locs(src_path, *, strict=False):
    """Load cached localizations, optionally requiring and validating the file.

    With ``strict=False``, a missing cache returns ``(None, None)`` for
    legacy callers. Strict mode raises on missing or malformed cache files.
    Set ``compressed=False`` in :func:`save_locs` to avoid decompression cost
    for caches intended for repeated direct loads.
    """
    path = result_path(src_path)
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(f'localization cache not found: {path}')
        return None, None

    try:
        with np.load(path, allow_pickle=False) as npz:
            if '__meta__' not in npz.files:
                raise ValueError('missing __meta__')
            meta = json.loads(str(npz['__meta__'].item()))
            fields = meta.get('fields')
            if not isinstance(fields, list):
                raise ValueError('invalid fields metadata')
            required = {'frame', 'x', 'y', 'photons', 'bg'}
            if not required.issubset(fields):
                raise ValueError('missing required localization fields')
            if strict and meta.get('format_version') != 1:
                raise ValueError('unsupported or missing cache format version')

            kw = {}
            lengths = set()
            for name in fields:
                if name not in npz.files:
                    raise ValueError(f'missing array {name!r}')
                value = np.asarray(npz[name])
                if value.ndim == 0:
                    raise ValueError(f'localization array {name!r} is scalar')
                lengths.add(value.shape[0])
                kw[name] = value
            if len(lengths) != 1:
                raise ValueError('localization arrays have different lengths')

            extra = dict(meta.get('extra_scalar', {}))
            for name in meta.get('extra_array', []):
                key = 'extra__' + name
                if key not in npz.files:
                    raise ValueError(f'missing extra array {name!r}')
                extra[name] = np.asarray(npz[key])
            kw['extra'] = extra
    except (KeyError, OSError, ValueError, zipfile.BadZipFile) as exc:
        if strict:
            raise ValueError(f'invalid localization cache {path}: {exc}') from exc
        raise

    return Localizations(**kw), meta


def save_cached_locs(locs, stack, src_path, *, method=None, filter_kw=None,
                     compressed=False):
    """Save filtered localizations with the image signature used by run_cached."""
    frame = np.asarray(locs.frame)
    if (frame.size and
            (np.any(frame < 0) or np.any(frame >= int(stack.shape[0])))):
        raise ValueError('localization frame indices outside the input stack')

    return save_locs(
        locs,
        src_path,
        method=method,
        extra_meta={
            'localization_stage': 'filtered',
            'image_signature': image_signature(stack),
            'filter_config': describe(filter_kw or {}),
            'pixel_size': (None if stack.pixel_size is None
                           else float(stack.pixel_size)),
        },
        compressed=compressed,
    )


# ---------------------------------------------------------------------------
# Strict cache reader. Callers explicitly choose cache read versus recompute.
# ---------------------------------------------------------------------------
def run_cached(stack, src_path, *, read_cache=False):
    """Load cached filtered localizations after validating the image stack.

    With ``read_cache=False`` this is a no-op and returns ``(None, None)``.
    When enabled, missing, legacy, malformed, or image-incompatible caches
    raise instead of silently triggering recomputation.
    """
    if not isinstance(read_cache, bool):
        raise TypeError('read_cache must be a bool')
    if not read_cache:
        return None, None
    if src_path is None:
        raise ValueError('src_path is required when read_cache=True')

    locs, meta = load_locs(src_path, strict=True)
    cache_info = meta.get('meta')
    if not isinstance(cache_info, dict):
        raise ValueError(
            f'localization cache {result_path(src_path)} is missing '
            'validation metadata')
    if cache_info.get('localization_stage') != 'filtered':
        raise ValueError(
            f'localization cache {result_path(src_path)} does not contain '
            'filtered localizations; recompute with read_cache=False')
    if cache_info.get('image_signature') != image_signature(stack):
        raise ValueError(
            f'localization cache {result_path(src_path)} does not match '
            'the input image shape, dtype, or first frame')
    if (locs.frame.size and
            (np.any(locs.frame < 0)
             or np.any(locs.frame >= int(stack.shape[0])))):
        raise ValueError(
            f'localization cache {result_path(src_path)} contains frame '
            'indices outside the input stack')
    return locs, meta