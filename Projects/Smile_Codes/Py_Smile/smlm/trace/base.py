from abc import ABC, abstractmethod
from ..datatypes import Trace, StateInfo, StateResult
import numpy as np

class TraceExtractor(ABC):
    """从 emitter 位置 + 原始图像栈提取稠密强度轨迹。"""

    @abstractmethod
    def extract(self, emitters, movie: np.ndarray,
                n_frames: int | None = None) -> list[Trace]:
        ...

# ======================================================================
# 分析器基类
# ======================================================================
class StateAnalyzer(ABC):
    method_name = "base"

    def __init__(self, **opts):
        defaults = self.default_opts()
        unknown = set(opts) - set(defaults)          # 用户传了、默认里没有的键
        if unknown:
            raise TypeError(
                f"{self.method_name} 不认识参数 {sorted(unknown)}. "
                f"可用: {sorted(defaults)}"
            )
        self.opts = {**defaults, **opts}

    def param_specs(self) -> dict:
        """参数规格: name -> (default, help)。子类用 super() 累加。"""
        return dict(
            use            = ("intensity", "分析哪个字段: intensity / raw"),
            bleach_tail    = (10,          "判定漂白时检查的尾部帧数"),
            bg_threshold   = (None,        "背景阈值: None=自动, 0=禁用, >0=固定值"),
            bg_auto_sigma  = (3.0,         "自动背景阈值的 σ 倍数"),
            min_state_sep  = (2.0,         "相邻态最小间隔 (× noise)"),
            bg_quantile    = (0.10,        "估计背景时取的低分位数"),
            verify_background = (True,     "是否二次校验背景态"),
            verbose        = (False,       "逐条打印进度"),
        )

    def default_opts(self) -> dict:
        # 从 specs 自动抽出默认值, 不再手写第二份
        return {k: v[0] for k, v in self.param_specs().items()}

    # ---- 子类唯一必须实现的方法 --------------------------------------
    @abstractmethod
    def _fit_sequence(self, y: np.ndarray) -> np.ndarray:
        """返回逐帧状态标签 (>0, 编号任意)。基类负责其余全部后处理。"""
        ...

    # ---- 通用后处理流水线 -------------------------------------------
    def analyze(self, trace) -> StateResult:
        o = self.opts
        y = np.asarray(getattr(trace, o["use"]), dtype=float)
        F = y.size

        seq = self._fit_sequence(y)
        seq = _sort_by_mean(seq, y)
        seq = merge_close_states(seq, y, o)
        sframe, bframe, seq = detect_bleach(seq, y, o["bleach_tail"])

        if o["verify_background"]:
            seq, sframe, bframe = self._verify_background(seq, y, sframe, bframe)

        seq = _renumber(seq)
        n_act = int(seq.max()) if np.any(seq > 0) else 0
        active_fr = F if bframe is None else bframe

        # stateInfo + 阶梯拟合
        fit = np.zeros(F)
        infos = []
        for j in range(1, n_act + 1):
            m = seq == j
            mi, si = float(y[m].mean()), float(y[m].std())
            infos.append(StateInfo(j, mi, si, float(m.sum()) / max(active_fr, 1)))
            fit[m] = mi
        if np.any(seq == 0):
            fit[seq == 0] = float(y[seq == 0].mean())

        # transitions
        chg = np.where(np.diff(seq) != 0)[0] + 1
        tr = np.array([[int(f), int(seq[f])] for f in chg], dtype=int).reshape(-1, 2)

        # background
        if np.any(seq == 0):
            bg_val = float(y[seq == 0].mean())
        elif infos:
            bg_val = min(i.mean_intensity for i in infos)
        else:
            bg_val = 0.0

        lifetime = None if bframe is None else int(bframe - sframe)

        return StateResult(
            emitter_id=getattr(trace, "emitter_id", -1),
            trace=y, sequence=seq, transitions=tr, n_states=n_act,
            state_info=infos, fit_trace=fit, background=bg_val,
            trace_bgsub=y - bg_val, fit_trace_bgsub=fit - bg_val,
            bleach_frame=bframe, start_frame=sframe, lifetime=lifetime,
            method=self.method_name,
            meta={"frames": np.asarray(trace.frames), "used_field": o["use"]},
        )

    def analyze_many(self, traces) -> list[StateResult]:
        out = []
        n = len(traces)
        for i, t in enumerate(traces):
            if self.opts["verbose"]:
                print(f"  trace {i+1}/{n} ({self.method_name})")
            out.append(self.analyze(t))
            if (i + 1) % 10 == 0 or i + 1 == n:
                print(f"Progress: {100*(i+1)/n:.1f}%")
        return out

    def _verify_background(self, seq, y, sframe, bframe):
        """若所有 active 态均值都低于背景阈值, 判定整条为背景。"""
        if not np.any(seq > 0):
            return seq, sframe, bframe
        thr = _bg_threshold(y, self.opts, robust_noise(y))
        if thr == -np.inf:
            return seq, sframe, bframe
        u = np.unique(seq[seq > 0])
        if all(y[seq == s].mean() <= thr for s in u):
            return np.zeros_like(seq), 0, None
        return seq, sframe, bframe

# ======================================================================
# 通用工具
# ======================================================================
def robust_noise(y: np.ndarray) -> float:
    """基于逐帧差分的 MAD 噪声估计 (对阶跃鲁棒)。"""
    if y.size < 2:
        return max(float(np.std(y)), 1e-6)
    dy = np.diff(y)
    n = np.median(np.abs(dy - np.median(dy))) / 0.6745 / np.sqrt(2)
    if n < 1e-6:
        n = np.std(y) * 0.1
    return max(float(n), 1e-9)


def _bg_threshold(y, opts, noise):
    bt = opts["bg_threshold"]
    if bt is None:
        sy = np.sort(y)
        n_low = max(round(y.size * opts["bg_quantile"]), 1)
        bg_level = sy[:n_low].mean()
        return bg_level + opts["bg_auto_sigma"] * noise
    if bt == 0:
        return -np.inf
    return float(bt)


def _sort_by_mean(seq, y):
    """把状态标签按均值升序重排 (1=最暗)。"""
    u = np.unique(seq[seq > 0])
    if u.size == 0:
        return seq
    means = np.array([y[seq == s].mean() for s in u])
    order = np.argsort(means)
    out = seq.copy()
    for rank, s in enumerate(u[order], start=1):
        out[seq == s] = -rank            # 用负号避免覆盖冲突
    out[out < 0] *= -1
    return out


def _renumber(seq):
    u = np.unique(seq[seq > 0])
    out = seq.copy()
    for j, lab in enumerate(u, start=1):
        out[seq == lab] = j
    return out


def merge_close_states(seq, y, opts):
    """合并背景态 + 均值差 < min_state_sep*noise 的相邻态。"""
    seq = seq.copy()
    if not np.any(seq > 0):
        return seq
    noise = robust_noise(y)
    bg_thr = _bg_threshold(y, opts, noise)

    # 1) 把所有 <=背景阈值 的态并成一个
    u = np.unique(seq[seq > 0])
    smean = np.array([y[seq == s].mean() for s in u])
    bg_labels = u[smean <= bg_thr]
    if bg_labels.size > 1:
        for lab in bg_labels[1:]:
            seq[seq == lab] = bg_labels[0]

    # 2) 迭代合并最接近的两个态
    while True:
        u = np.unique(seq[seq > 0])
        if u.size < 2:
            break
        smean = np.array([y[seq == s].mean() for s in u])
        slen = np.array([np.sum(seq == s) for s in u])
        order = np.argsort(smean)
        u, smean, slen = u[order], smean[order], slen[order]
        diffs = np.diff(smean)
        idx = int(np.argmin(diffs))
        if diffs[idx] < opts["min_state_sep"] * noise:
            a, b = u[idx], u[idx + 1]
            if slen[idx] >= slen[idx + 1]:
                seq[seq == b] = a
            else:
                seq[seq == a] = b
        else:
            break
    return seq


def detect_bleach(seq, y, min_tail):
    """
    识别起始帧与漂白帧。最暗的态视为背景标签:
    前导背景 -> 0; 若尾部背景连续 >= min_tail, 则记漂白帧并置 0。
    返回 (start_frame, bleach_frame|None, seq)  索引 0 基。
    """
    seq = seq.copy()
    F = seq.size
    u = np.unique(seq[seq > 0])
    if u.size == 0:
        seq[:] = 0
        return F, F, seq

    smean = np.array([y[seq == s].mean() for s in u])
    bg_label = u[int(np.argmin(smean))]

    if u.size == 1:                       # 单态: 无漂白
        return 0, None, seq

    sframe = 0
    if seq[0] == bg_label:
        first = np.where((seq != bg_label) & (seq > 0))[0]
        if first.size == 0:
            seq[:] = 0
            return F, F, seq
        sframe = int(first[0])
        seq[:sframe] = 0

    last = np.where((seq != bg_label) & (seq > 0))[0]
    if last.size == 0:
        seq[seq > 0] = 0
        return F, F, seq
    last = int(last[-1])

    bframe = None
    if (F - 1 - last) >= min_tail:
        bframe = last + 1
        seq[bframe:] = 0
    return sframe, bframe, seq


# ======================================================================
# 注册表 / 工厂
# ======================================================================
_ANALYZERS: dict[str, type[StateAnalyzer]] = {}   # 空表, 由装饰器填充

def register_analyzer(cls):
    _ANALYZERS[cls.method_name] = cls
    return cls

def get_analyzer(method="changepoint", **opts) -> StateAnalyzer:
    key = method.lower()
    if key not in _ANALYZERS:
        raise ValueError(f"未知方法 '{method}'. 可用: {list(_ANALYZERS)}")
    return _ANALYZERS[key](**opts)

def analyze_traces(traces, method="changepoint", **opts) -> list[StateResult]:
    return get_analyzer(method, **opts).analyze_many(list(traces))