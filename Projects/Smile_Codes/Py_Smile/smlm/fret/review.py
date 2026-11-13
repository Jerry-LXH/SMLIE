# smlm/fret/review.py
"""人工复核的**状态**部分: 候选队列、保留/剔除、缓存、落盘。

本模块**不 import ipywidgets** —— 复核逻辑(headless 也能跑、能写进测试)与
交互层分开, 见 ``widgets.py``。

为什么要落盘
------------
``KEEP_IDS`` 原来只是 kernel 里的一个 list, 一旦重启就没了, 而人工看几十条 trace
是很贵的操作。所以每次决定都写一份 ``<src>_review.json``(带参数指纹), 下次自动
加载 —— 参数变了指纹就不匹配, 会提示而不是静默沿用。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import os
import numpy as np

from .datatypes import StepEvent
from . import viz


def _fingerprint(objs) -> str:
    s = json.dumps(objs, sort_keys=True, default=str)
    return hashlib.sha1(s.encode()).hexdigest()[:12]


@dataclass
class ReviewSession:
    """一次复核的完整状态。

    Parameters
    ----------
    candidates : list[StepEvent]
        待复核的候选(已按 ``frame`` 排序)。
    trace_a, trace_d : list[Trace]
        受体 / 供体通道的 trace, 下标 = emitter_id。
    centres : dict
        ``{'a': (N,2), 'd': (N,2)}``, 各通道的孔径中心 ``[row, col]``。
    images : dict
        ``{'a': (T,H,W), 'd': (T,H,W)}`` 漂移校正后的图像栈。
    """
    candidates: list
    trace_a: list
    trace_d: list
    centres: dict = field(default_factory=dict)
    images: dict = field(default_factory=dict)
    ex_time: float = 1.0
    movie_kwargs: dict = field(default_factory=dict)
    trace_kwargs: dict = field(default_factory=dict)
    mle: dict = field(default_factory=dict)      # {通道: {emitter_id: (T,)}}
    index: int = 0
    decisions: dict = field(default_factory=dict)   # emitter_id -> bool
    path: str = None

    _movies: dict = field(default_factory=dict, repr=False)
    _traces: dict = field(default_factory=dict, repr=False)

    # ---------- 队列 ----------
    def __len__(self):
        return len(self.candidates)

    @property
    def finished(self):
        return self.index >= len(self.candidates)

    def current(self):
        return None if self.finished else self.candidates[self.index]

    @property
    def keep_ids(self):
        return sorted(i for i, v in self.decisions.items() if v)

    @property
    def reject_ids(self):
        return sorted(i for i, v in self.decisions.items() if not v)

    # ---------- 决策 ----------
    def decide(self, keep: bool):
        ev = self.current()
        if ev is None:
            return
        self.decisions[int(ev.emitter_id)] = bool(keep)
        self.index += 1
        if self.path:
            self.save()

    def back(self):
        if self.index == 0:
            return
        self.index -= 1
        ev = self.candidates[self.index]
        self.decisions.pop(int(ev.emitter_id), None)
        if self.path:
            self.save()

    # ---------- 渲染(带缓存) ----------
    def movie(self, ev=None):
        ev = ev or self.current()
        if ev is None:
            return [], np.zeros(0, int)
        i = int(ev.emitter_id)
        if i not in self._movies:
            kw = dict(half=20, step=10,
                      radius_d=self.trace_d[i].meta.get('radius', 5),
                      radius_a=self.trace_a[i].meta.get('radius', 4))
            kw.update(self.movie_kwargs)
            self._movies[i] = viz.render_movie_pngs(
                self.images['d'], self.images['a'],
                self.centres['d'][i], self.centres['a'][i], **kw)
        return self._movies[i]

    def trace_png(self, ev=None):
        ev = ev or self.current()
        if ev is None:
            return b''
        i = int(ev.emitter_id)
        if i not in self._traces:
            kw = dict(ex_time=self.ex_time)
            kw.update(self.trace_kwargs)
            self._traces[i] = viz.render_trace_png(
                self.trace_a[i], self.trace_d[i], ev,
                mle_a=(self.mle.get('a') or {}).get(i),
                mle_d=(self.mle.get('d') or {}).get(i), **kw)
        return self._traces[i]

    # ---------- 落盘 ----------
    def params(self):
        return {'n_candidates': len(self.candidates),
                'ex_time': self.ex_time,
                'movie': self.movie_kwargs,
                'trace': self.trace_kwargs}

    def save(self, path=None):
        path = path or self.path
        if not path:
            return None
        d = {'version': 1, 'fingerprint': _fingerprint(self.params()),
             'params': self.params(), 'index': int(self.index),
             'keep': self.keep_ids, 'reject': self.reject_ids}
        d0 = os.path.dirname(os.path.abspath(path))
        if d0:
            os.makedirs(d0, exist_ok=True)
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(d, f, ensure_ascii=False, indent=1)
        self.path = path
        return path

    def load(self, path=None, strict=True):
        """读回已存的决策。参数指纹不匹配时返回 False(除非 ``strict=False``)。"""
        path = path or self.path
        if not path or not os.path.exists(path):
            return False
        with open(path, encoding='utf-8') as f:
            d = json.load(f)
        if strict and d.get('fingerprint') != _fingerprint(self.params()):
            print(f'[review] {os.path.basename(path)} 的参数指纹不匹配, 未加载')
            return False
        self.path = path
        self.decisions = {int(k): True for k in d.get('keep', [])}
        self.decisions.update({int(k): False for k in d.get('reject', [])})
        self.index = int(d.get('index', 0))
        return True

    # ---------- 汇总 ----------
    def summary(self):
        return (f'已审 {len(self.decisions)} / {len(self.candidates)} 条   '
                f'保留 {len(self.keep_ids)} / 剔除 {len(self.reject_ids)}')
