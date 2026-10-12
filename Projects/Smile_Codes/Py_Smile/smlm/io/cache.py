# smlm/io/cache.py
import os
import json
import hashlib
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

def data_fingerprint(stack):
    """对输入 ImageStack 算内容指纹: 覆盖一切上游改动(preprocess/camera)。"""
    arr = np.ascontiguousarray(stack.data)
    h = hashlib.blake2b(digest_size=16)
    h.update(str(arr.shape).encode())
    h.update(str(arr.dtype).encode())
    h.update(arr.tobytes())                      # 全量哈希, 最可靠
    # 可读的上游元信息也一并纳入(便于排查, 也防纯 metadata 改动)
    h.update(json.dumps(describe(stack.metadata),
                        sort_keys=True, default=str).encode())
    return h.hexdigest()

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
def save_locs(locs, src_path, method=None, data_hash=None, extra_meta=None):
    path = result_path(src_path)
    arrays = {}
    meta = {'fields': [], 'extra_array': [], 'extra_scalar': {}}

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
    np.savez_compressed(path, **arrays)
    print(f"[cache] saved {len(locs)} locs -> {path}")
    return path


def load_locs(src_path):
    """返回 (Localizations, meta) 或 (None, None)。"""
    path = result_path(src_path)
    if not os.path.exists(path):
        return None, None

    npz = np.load(path, allow_pickle=False)
    meta = json.loads(str(npz['__meta__']))

    kw = {f: npz[f] for f in meta['fields']}
    extra = dict(meta.get('extra_scalar', {}))
    for k in meta.get('extra_array', []):
        extra[k] = npz['extra__' + k]
    kw['extra'] = extra

    return Localizations(**kw), meta


# ---------------------------------------------------------------------------
#  带缓存的运行入口
# ---------------------------------------------------------------------------
def combined_key(method, stack):
    return {
        'config_hash': config_hash(describe(method)),
        'data_hash': data_fingerprint(stack),
    }


def run_cached(method, stack, src_path, force=False):
    key = combined_key(method, stack)

    if not force:
        locs, meta = load_locs(src_path)
        if locs is not None:
            cached = {'config_hash': meta.get('config_hash'),
                      'data_hash': meta.get('data_hash')}
            if cached == key:
                print(f"[cache] hit -> {result_path(src_path)} ({len(locs)} locs)")
                return locs
            # 告诉用户到底是哪一环变了
            if cached['config_hash'] != key['config_hash']:
                print("[cache] method config changed, recompute.")
            if cached['data_hash'] != key['data_hash']:
                print("[cache] input stack changed (preprocess/camera), recompute.")

    locs = method.run(stack)
    save_locs(locs, src_path, method=method, data_hash=key['data_hash'])
    return locs