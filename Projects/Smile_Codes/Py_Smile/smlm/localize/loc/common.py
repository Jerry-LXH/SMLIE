# smlm/localize/loc/common.py
"""各定位方案共享的基础设施: ROI 提取 + 多进程 worker 池。

全局约定(所有方案一致):
  - stack.data: (T, H, W); y = row, x = col; 坐标为 0-based 像素中心。
  - 每个方案的 fit 函数签名:
        fit_single(roi, *fit_args) -> 相对 ROI 中心的一维估计数组, 失败返回 None。
  - 多进程 worker 通过 (fit_func, fit_args) 传递; fit_func 必须是模块级
    函数(可 pickle), 以便跨进程分发。
"""
from __future__ import annotations
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm


def extract_rois(stack, cand, edge):
    """从 (T, H, W) 栈中抠出 (2*edge+1) x (2*edge+1) 的 ROI, 过滤越界候选。

    返回 (rois: list[np.ndarray], meta: (n,3) int [frame, row, col])。
    """
    data = stack.data
    _, H, W = data.shape
    e = int(edge)
    frames = np.asarray(cand.frame, dtype=int)
    rows = np.asarray(cand.y, dtype=int)
    cols = np.asarray(cand.x, dtype=int)

    rois, meta = [], []
    for f, rr, cc in zip(frames, rows, cols):
        if rr - e < 0 or cc - e < 0 or rr + e > H - 1 or cc + e > W - 1:
            continue
        rois.append(data[f, rr - e:rr + e + 1, cc - e:cc + e + 1].astype(float))
        meta.append((f, rr, cc))
    meta = np.array(meta, dtype=int) if meta else np.empty((0, 3), int)
    return rois, meta


def _worker(payload):
    """多进程入口: (fit_func, fit_args) -> 估计数组或 None。"""
    fit_func, fit_args = payload
    try:
        return fit_func(*fit_args)
    except Exception:
        return None


def run_pool(fit_func, rois, fit_args, n_jobs=1, desc="fitting", verbose=True):
    """对每个 ROI 调用 fit_func(roi, *fit_args), 支持多进程。

    返回与 rois 等长的列表, 失败项为 None。
    """
    payload = [(fit_func, (r, *fit_args)) for r in rois]
    t0 = time.time()
    if n_jobs > 1 and len(rois) > 0:
        with ProcessPoolExecutor(max_workers=int(n_jobs)) as ex:
            ests = list(tqdm(ex.map(_worker, payload, chunksize=64),
                             total=len(rois), disable=not verbose, desc=desc))
    else:
        ests = [_worker(p) for p in tqdm(payload, disable=not verbose, desc=desc)]
    dt = time.time() - t0
    if verbose:
        per = dt / max(len(rois), 1) * 1000
        print(f"{desc} done in {dt:.2f}s ({per:.2f} ms/loc)")
    return ests
