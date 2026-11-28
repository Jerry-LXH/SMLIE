"""trace 数据的存取 —— 目前只有池化 npz 的读取。

池化文件的契约
--------------
``traces_pooled.npz`` 要把帧数不同的视野装进一个 ``(N, Tmax)`` 数组, 所以短的
trace 是**在尾部补 NaN** 对齐的, 真实长度存在 ``n_frames`` 里。

读回来必须先按 ``n_frames`` 截断, 否则::

    np.cumsum([1, 2, nan, ...])  -> 全是 nan
    changepoint.segment(...)    -> 一个切点都没有, 整条判成"一个暗段"

而且**不抛任何异常** —— 实测 139 条 SPARK trace 全部得到 ``n_states=1``。
:func:`load_pooled` 就是来堵这个坑的: 它只返回截断后的 trace,
``len(trace['frames'])`` 就是真实帧数。
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..datatypes import Trace

POOLED_KEYS = ('intensity', 'bg_flux', 'emitter_id', 'file_id', 'files',
               'n_frames', 't_s', 'exposure_s')


def _lengths(intensity, n_frames):
    """每条 trace 的真实帧数: 优先 ``n_frames``, 缺了就找最后一个有限值。"""
    n, T = intensity.shape
    if n_frames is not None:
        return np.clip(np.asarray(n_frames, int).reshape(-1), 0, T)
    finite = np.isfinite(intensity)
    has = finite.any(axis=1)
    last = T - 1 - np.argmax(finite[:, ::-1], axis=1)
    return np.where(has, last + 1, 0)


def load_pooled(path, *, channel=None, verbose=False):
    """读 ``traces_pooled.npz`` -> ``(traces, meta)``。

    Parameters
    ----------
    path : str | Path
        池化文件, 例如 ``TRACE_pooled/SPARK_6mW/traces_pooled.npz``。
    channel : str, optional
        写进每条 trace 的 ``channel`` 字段; 默认取父目录名 (组名),
        便于 :func:`smlm.trace.plot_trace_fields` 用 ``groups`` 分组。
    verbose : bool
        打印条数与帧数分布。

    Returns
    -------
    traces : list[Trace]
        每条含 ``emitter_id / file_id / file / channel / frames / t_s /
        intensity / bg_flux / exposure_s``, 长度已按 ``n_frames`` 截断。
    meta : dict
        ``path`` / ``files`` (来源目录名) / ``t_s`` (公共时间轴) /
        ``exposure_s`` / ``n_frames`` (每条的真实帧数) / ``n_columns`` (数组宽度)。

    Examples
    --------
    >>> traces, meta = load_pooled('TRACE_pooled/SPARK_6mW/traces_pooled.npz')
    >>> len(traces), meta['exposure_s']
    (139, 0.1)
    """
    path = Path(path)
    with np.load(path) as z:
        I = np.asarray(z['intensity'], np.float32)
        if I.ndim != 2:
            raise ValueError(f'{path.name}: intensity 必须是 (N, T), 得到 {I.shape}')
        B = (np.asarray(z['bg_flux'], np.float32)
             if 'bg_flux' in z.files else np.full_like(I, np.nan))
        ids = (np.asarray(z['emitter_id'], int) if 'emitter_id' in z.files
               else np.arange(len(I)))
        fid = (np.asarray(z['file_id'], int) if 'file_id' in z.files
               else np.zeros(len(I), int))
        files = [str(s) for s in z['files']] if 'files' in z.files else []
        nf = np.asarray(z['n_frames'], int) if 'n_frames' in z.files else None
        t_s = (np.asarray(z['t_s'], float) if 't_s' in z.files
               else np.arange(I.shape[1], dtype=float))
        ex = float(z['exposure_s']) if 'exposure_s' in z.files else (
            float(t_s[1] - t_s[0]) if t_s.size > 1 else np.nan)

    lens = _lengths(I, nf)
    name = str(channel) if channel is not None else path.parent.name

    traces = []
    bad = 0
    for i in range(len(I)):
        n = int(lens[i])
        y = I[i, :n]
        bad += int(np.isnan(y).any())
        traces.append(Trace(
            emitter_id=int(ids[i]),
            file_id=int(fid[i]),
            file=files[fid[i]] if 0 <= fid[i] < len(files) else '',
            channel=name,
            frames=np.arange(n),
            t_s=t_s[:n],
            intensity=y,
            bg_flux=B[i, :n],
            exposure_s=ex,
        ))

    meta = dict(path=path, files=files, t_s=t_s, exposure_s=ex,
                n_frames=lens, n_columns=int(I.shape[1]))
    if verbose:
        lo, hi = (int(lens.min()), int(lens.max())) if lens.size else (0, 0)
        print(f'  {path.parent.name}: {len(traces)} 条, 帧数 {lo}~{hi} '
              f'(数组 {I.shape[1]} 列, 尾部 NaN 已按 n_frames 截掉)')
    if bad:
        print(f'  [!] {bad} 条 trace 在 n_frames 以内仍有 NaN, '
              '状态判定会直接报错; 检查池化时的 intensity')
    return traces, meta
