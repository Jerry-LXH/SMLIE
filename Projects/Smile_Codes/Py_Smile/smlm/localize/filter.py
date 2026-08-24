import numpy as np

def filter_locs(locs, photons=(20, 1e7), bg=(0, 1e5),
                sigma=(0.5, 5.0), lp=(0.0, 1.0),
                z=None, z_lp=None,
                return_mask=False, verbose=True):
    """按物理量合理范围过滤, 对应 MATLAB filterBadLocs。

    每项传 (min, max); min/max 任一为 None 表示该侧不限制;
    整项传 None 表示跳过该字段(比如没有 sigma 时)。
    区间为开区间 (>lo, <hi), 与 MATLAB 版本一致。
    z / z_lp 为 3D 定位的可选过滤项, 2D 定位 (字段为 None) 时自动跳过。
    """
    n0 = len(locs)
    mask = np.ones(n0, dtype=bool)
    ranges = {'photons': photons, 'bg': bg, 'sigma': sigma, 'lp': lp,
              'z': z, 'z_lp': z_lp}

    for fld, rng in ranges.items():
        if rng is None:
            continue
        v = getattr(locs, fld, None)
        if v is None:                 # 该字段本身缺失, 跳过
            continue
        v = np.asarray(v, dtype=float)
        lo, hi = rng
        mask &= np.isfinite(v)        # 顺便剔除 NaN/Inf
        if lo is not None:
            mask &= v > lo
        if hi is not None:
            mask &= v < hi

    out = locs[mask]
    if verbose:
        n1 = len(out)
        print(f"[filter] removed {n0 - n1} / {n0} "
              f"({100 * (n0 - n1) / max(n0, 1):.2f}%), kept {n1}")
    return (out, mask) if return_mask else out