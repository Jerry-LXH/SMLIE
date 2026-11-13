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
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from tqdm import tqdm


def extract_rois(stack, cand, edge):
    """从 (T, H, W) 栈中抠出 (2*edge+1) x (2*edge+1) 的 ROI, 过滤越界候选。

    返回 (rois: list[np.ndarray], meta: (n,3) int [frame, row, col])。
    ROI 以 float32 保存 (定位精度足够), 相比 float64 可把这块内存减半,
    对动辄几十万候选的多进程场景很重要。
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
        rois.append(data[f, rr - e:rr + e + 1, cc - e:cc + e + 1].astype(np.float32))
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


def _payload_iter(fit_func, rois, fit_args, start, stop):
    """惰性生成 (fit_func, args) 任务, 避免一次性构造几十万条 tuple。"""
    for i in range(start, stop):
        yield (fit_func, (rois[i], *fit_args))


def run_pool(fit_func, rois, fit_args, n_jobs=1, desc="fitting", verbose=True):
    """对每个 ROI 调用 fit_func(roi, *fit_args), 支持多进程。

    返回与 rois 等长的列表, 失败项为 None。

    健壮性: Windows 下若系统提交内存不足(例如其它进程泄漏), 子进程可能
    在导入 numpy/scipy 时就直接被杀, 表现为 BrokenProcessPool。这里按块
    提交任务, 一旦检测到进程池不可用, 就自动改用线程池把剩余 ROI 跑完,
    不再让整个定位流程失败。
    """
    n = len(rois)
    t0 = time.time()
    ests = [None] * n

    if n_jobs > 1 and n > 0:
        n_jobs = int(n_jobs)
        block = max(8192, 512 * n_jobs)   # 分块提交, 便于中途降级
        done = 0
        ex = None
        try:
            ex = ProcessPoolExecutor(max_workers=n_jobs)
            for start in range(0, n, block):
                stop = min(start + block, n)
                ests[start:stop] = list(tqdm(
                    ex.map(_worker,
                           _payload_iter(fit_func, rois, fit_args, start, stop),
                           chunksize=64),
                    total=stop - start, disable=not verbose, desc=desc))
                done = stop
            ex.shutdown(wait=True)
            ex = None
        except (BrokenProcessPool, OSError, RuntimeError) as exc:
            if ex is not None:
                try:
                    ex.shutdown(wait=False, cancel_futures=True)
                except Exception:
                    pass
            if verbose:
                print(f"[run_pool] 进程池不可用 ({type(exc).__name__}: {exc}); "
                      f"改用线程池继续剩余 {n - done} 个 ROI. "
                      f"若提示'页面文件太小', 请检查是否有进程占满系统提交内存.")
            with ThreadPoolExecutor(max_workers=n_jobs) as tex:
                for start in range(done, n, block):
                    stop = min(start + block, n)
                    ests[start:stop] = list(tqdm(
                        tex.map(_worker,
                                _payload_iter(fit_func, rois, fit_args, start, stop),
                                chunksize=64),
                        total=stop - start, disable=not verbose, desc=desc))
    else:
        ests = [_worker(p) for p in tqdm(
            _payload_iter(fit_func, rois, fit_args, 0, n),
            total=n, disable=not verbose, desc=desc)]

    dt = time.time() - t0
    if verbose:
        per = dt / max(n, 1) * 1000
        print(f"{desc} done in {dt:.2f}s ({per:.2f} ms/loc)")
    return ests
