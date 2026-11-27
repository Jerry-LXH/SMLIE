import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
#  逐字段曲线 (raw / bg / snr / ...)
# ---------------------------------------------------------------------------
# 只给几个已知字段标单位, 其余留空。
# 注意 raw/intensity/bg_flux 在 mode='matched' 下是 a.u. (≈ 峰值×2πσ²), 不是光子数;
# 只有 mode='hard' 时它们才等于孔径内的光子数。
_FIELD_UNITS = {
    "raw": "a.u.",
    "intensity": "a.u.",
    "fit": "a.u.",
    "bg_flux": "a.u.",
    "bg": "ADU/px",
    "bg_std": "ADU/px",
    "snr": "SNR",
}


def _one_trace(obj):
    """把传入的单个 trace 包成 list; 序列原样返回。"""
    if isinstance(obj, dict):
        return [obj]
    return list(obj)


_PALETTE = ("tab:blue", "tab:red", "tab:green", "tab:orange", "tab:purple",
            "tab:brown", "tab:pink", "0.4")
_LINESTYLES = ("-", "--", ":", "-.")

# 状态色带: state 是**亮度等级** (0 = 最暗), 等级越高绿色越深;
# 漂白段用 bleach_frame 单独叠加斜纹, 不再占用状态码。
_BAND_LEVEL_CMAP = "Greens"
_BAND_ALPHA = 0.45


def _level_colors(n_levels):
    """等级 -> 颜色: 0 = 灰, 1..k-1 = 逐渐加深的绿。"""
    colors = {0: "0.80"}
    if n_levels <= 1:
        return colors
    cmap = plt.get_cmap(_BAND_LEVEL_CMAP)
    for j in range(1, n_levels):
        frac = 0.35 + 0.55 * (j - 1) / max(n_levels - 2, 1)
        colors[j] = cmap(frac)
    return colors


def _state_runs(mask):
    """连续 ``True`` 段: ``[(start, stop_exclusive), ...]``。"""
    mask = np.asarray(mask, bool).reshape(-1)
    if not mask.any():
        return []
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def _field_label(name):
    unit = _FIELD_UNITS.get(name)
    return f"{name} ({unit})" if unit else name


def _clamp_bottom(ax, bottom):
    """把纵轴下限压到 ``bottom``, 上限保持自适应。

    多条 trace 并排时, 若各自自适应, 零点位置会在行与行之间跳动, 看不出
    "有没有信号"; 统一从 0 起画才比得出来。数据整体低于 ``bottom`` 时给一个
    最小高度, 避免上下限相反。
    """
    _, hi = ax.get_ylim()
    if hi <= bottom:
        hi = float(bottom) + 1.0
    ax.set_ylim(float(bottom), hi)


def _normalize_group(spec):
    """分组说明 -> ``dict(name, fields, color, linestyle)``。

    ``spec`` 可以是 ``'raw'``、``('bg', 'bg_std')``, 或
    ``dict(fields=..., color=..., label=..., linestyle=...)``。
    """
    if isinstance(spec, str):
        spec = {"fields": (spec,)}
    elif isinstance(spec, (list, tuple)):
        spec = {"fields": tuple(spec)}
    if not isinstance(spec, dict):
        raise TypeError("分组必须是 str / tuple / dict, 得到 "
                        f"{type(spec).__name__}")
    fields = spec.get("fields")
    if isinstance(fields, str):
        fields = (fields,)
    fields = tuple(fields or ())
    if not fields:
        raise ValueError("每个分组至少要有一个字段")
    return {"name": spec.get("label") or spec.get("name") or "+".join(fields),
            "fields": fields,
            "color": spec.get("color"),
            "linestyle": spec.get("linestyle")}


def _normalize_groups(groups):
    if isinstance(groups, dict):                 # {组名: 字段说明}
        prepared = []
        for name, spec in groups.items():
            item = dict(spec) if isinstance(spec, dict) else {'fields': spec}
            item.setdefault('label', name)        # 字典的键默认当组名
            prepared.append(item)
        groups = prepared
    if not groups:
        raise ValueError("groups / fields 不能为空")
    return [_normalize_group(g) for g in groups]


def _line_styles(group, group_index, default_color=None):
    """把一组的字段展开成逐字段的 (color, linestyle)。"""
    n = len(group["fields"])
    color = group["color"] if group["color"] is not None else default_color
    if color is None:
        color = _PALETTE[group_index % len(_PALETTE)]
    if isinstance(color, str):
        colors = [color] * n
    else:
        colors = list(color)
        if len(colors) == 1:
            colors = colors * n
        elif len(colors) != n:
            raise ValueError(f"分组 {group['name']!r} 有 {n} 个字段, "
                             f"却给了 {len(colors)} 个颜色")
    if group["linestyle"] is not None:
        styles = [group["linestyle"]] * n
    elif n > 1 and len(set(colors)) == 1:
        styles = [_LINESTYLES[i % len(_LINESTYLES)] for i in range(n)]
    else:
        styles = ["-"] * n
    return list(zip(colors, styles))


def _draw_bands(ax, xv, codes, bleach_frame=-1):
    """按亮度等级画背景色带: 0 = 灰, ≥1 = 绿 (越多级越深); 漂白段叠斜纹。"""
    codes = np.asarray(codes).reshape(-1)
    n_levels = int(codes.max()) + 1 if codes.size else 1
    for code, color in _level_colors(n_levels).items():
        for start, stop in _state_runs(codes == code):
            ax.axvspan(xv[start], xv[stop - 1], color=color, alpha=_BAND_ALPHA,
                       lw=0, zorder=0)
    if codes.size and bleach_frame is not None and 0 <= bleach_frame < codes.size:
        ax.axvspan(xv[bleach_frame], xv[-1], facecolor="none",
                   edgecolor="0.25", hatch="///", alpha=_BAND_ALPHA, lw=0,
                   zorder=1)


def _format_scalar(name, value):
    """``name = value``, 整数去掉小数点。"""
    raw = np.asarray(value)
    if raw.ndim != 0:
        raise ValueError(f"字段 {name!r} 的形状是 {raw.shape}, 不是单个数字; "
                         "逐帧数组请放进 groups / fields")
    val = raw.item()
    if isinstance(val, (bool, np.bool_)):
        text = str(bool(val))
    elif isinstance(val, (int, np.integer)):
        text = str(int(val))
    elif isinstance(val, (float, np.floating)):
        text = f"{float(val):g}"
    else:
        text = str(val)
    return f"{name} = {text}"


def show_trace_movies(stack, traces, emitters, *, indices=None, select=None,
                      n_max=12, half=12, frames=None, verbose=True,
                      **movie_kw):
    """把筛出来的 trace 直接放成 emitter 小视频。

    ``traces[i]`` 与 ``emitters[i]`` 一一对应（都由同一份 centers 生成），所以
    trace 的索引**直接就是** emitter 的索引，不需要再建一张映射表。

    Parameters
    ----------
    stack : ImageStack
        漂移校正后的单通道图像。
    traces : sequence of Trace
    emitters : sequence
        与 ``traces`` 一一对应的 emitter 列表。
    indices : sequence of int, optional
        要放的 trace 索引。给了就用它，忽略 ``select``。
    select : dict, optional
        传给 :func:`smlm.trace.states.select_traces` 的筛选条件, 例如
        ``dict(n_on=(2, 6), bleach_frame=(0, None))``。
    n_max : int
        最多放几条 (GIF 越大越卡)。
    half, frames
        透传 ``smlm.viz.show_emitter_movies`` (窗口半径 / 帧范围)。
    verbose : bool
        输出时逐条打印该 trace 的状态统计。

    Returns
    -------
    IPython.display.Image
        notebook 里直接显示。

    Examples
    --------
    >>> show_trace_movies(corr_image, fill_traces, emitters,
    ...                   select=dict(n_on=(2, 6)), n_max=8,
    ...                   half=12, frames=slice(0, 400))
    """
    from ..viz import show_emitter_movies      # 延迟导入: viz 之间不互相依赖
    from .states import select_traces

    if not isinstance(traces, (list, tuple)):
        traces = list(traces)
    if indices is None:
        indices = select_traces(traces, **(select or {}))
    indices = [int(i) for i in indices][:max(int(n_max), 1)]
    if not indices:
        raise ValueError('没有 trace 可展示: 放宽 select 条件或先算 state')
    for i in indices:
        if not 0 <= i < len(emitters):
            raise IndexError(f'index {i} 超出 0..{len(emitters) - 1} '
                             '(traces 与 emitters 必须一一对应)')

    if verbose:
        for i in indices:
            print('  ' + _trace_summary(i, traces[i]))

    return show_emitter_movies(stack, emitters, idx=indices, half=half,
                               frames=frames, **movie_kw)


def _trace_summary(index, trace):
    """``#3  n_on=2  n_sw=3  bleach=95  on=0.38``。"""
    parts = [f'#{index}']
    if 'n_on' in trace:
        parts.append(f'n_on={int(trace["n_on"])}')
    if 'n_switches' in trace:
        parts.append(f'n_sw={int(trace["n_switches"])}')
    if 'bleach_frame' in trace:
        bf = int(trace['bleach_frame'])
        parts.append('bleach=-' if bf < 0 else f'bleach={bf}')
    if 'on_fraction' in trace:
        parts.append(f'on={float(trace["on_fraction"]):.2f}')
    return '  '.join(parts)


def plot_trace_fields(traces, fields=("raw", "bg"), *, groups=None, indices=None,
                      n_max=6, x=None, xlabel="Frame", title=None, figsize=None,
                      linewidth=0.8, color=None, legend=True, yscale=None,
                      scalars=None, bands=None, ylim_bottom=0.0):
    """画 trace 里若干个**逐帧字段**的曲线, 支持把多个字段叠在同一张图上。

    每行一条 trace, 每列一个**分组**; 同列共用时间轴。单个数 (``neighbours`` /
    ``n_bg_pixels`` / ``radius`` ...) 用 ``scalars`` 放进图例。

    Parameters
    ----------
    traces : Trace | dict | sequence
        单条 trace 或 trace 列表 (``list[Trace]``)。
    fields : sequence of str
        每个字段各占一列 (``groups`` 为 ``None`` 时使用)。
        字段必须是一次一条 trace 的逐帧数组 (长度 = 帧数); 标量请用 ``scalars``。
    groups : dict | sequence, optional
        分组说明。一组的字段叠在同一张图上::

            groups={'signal': dict(fields=('raw',), color='tab:blue'),
                    'background': dict(fields=('bg', 'bg_std'), color='tab:red')}

        每项可以是 ``'raw'`` (单字段)、``('bg', 'bg_std')`` (多字段),
        或 ``dict(fields=..., color=..., label=..., linestyle=...)``。
        一个组只给一个颜色时, 组内多条线用线型 (实线/虚线/点线) 区分;
        给一串颜色则逐字段上色。
    scalars : str | sequence of str, optional
        要显示成 ``名字 = 数值`` 的**单个数字**字段, 追加到该行每个子图的图例里。
    bands : str, optional
        状态字段名 (如 ``'state'``), 在每个子图里按**亮度等级**画背景色带
        (0 = 灰, ≥1 = 绿, 越多级越深); 若该 trace 有 ``bleach_frame``,
        漂白段再叠一层斜纹。
    indices : sequence of int, optional
        画哪几条; 默认前 ``n_max`` 条。
    n_max : int
        ``indices`` 为 ``None`` 时最多画几条。
    x : array-like, optional
        横轴; 默认用该 trace 的 ``frames``。
    color : str | sequence, optional
        组没有指定颜色时用的默认颜色。
    legend : bool
        是否列出字段名。``scalars`` 给了就一定画图例 (只列标量)。
    yscale : str, optional
        透传给 ``ax.set_yscale`` (如 ``'log'``)。给了它就不再压 ``ylim_bottom``
        (对数轴没有 0)。
    ylim_bottom : float | None
        纵轴下限 (默认 ``0``): 上限自适应, 下限统一压到 0 —— 多行并排时零点位置
        不跳, 才看得出"有没有信号"。``None`` = 完全自适应 (旧行为)。

    Returns
    -------
    matplotlib.figure.Figure

    Examples
    --------
    >>> plot_trace_fields(fill_traces, indices=range(4),
    ...                   groups={'raw': ('raw',),
    ...                           'bg': dict(fields=('bg', 'bg_std'),
    ...                                      color='tab:red')})
    """
    import matplotlib.pyplot as plt

    traces = _one_trace(traces)
    if not traces:
        raise ValueError("traces 是空的")

    specs = _normalize_groups(fields if groups is None else groups)
    # 颜色/线型先算好: 分组写错时在画图之前就报错
    group_styles = [_line_styles(g, i, color) for i, g in enumerate(specs)]

    indices = (list(range(min(len(traces), int(n_max)))) if indices is None
               else list(indices))
    if not indices:
        raise ValueError("indices 是空的")

    wanted = [name for g in specs for name in g["fields"]]
    scalar_names = ([scalars] if isinstance(scalars, str)
                    else list(scalars or []))

    for idx in indices:
        if not 0 <= idx < len(traces):
            raise IndexError(f"index {idx} 超出 0..{len(traces) - 1}")
        for name in wanted:
            if name not in traces[idx]:
                raise KeyError(f"trace #{idx} 没有字段 {name!r}; "
                               f"可用: {sorted(traces[idx])}")
            arr = np.asarray(traces[idx][name])
            if arr.ndim != 1:
                raise ValueError(
                    f"字段 {name!r} 的形状是 {arr.shape}, 不是逐帧序列 "
                    f"(长度 = 帧数); 单个数字请放进 scalars=(...)")
        for name in scalar_names:
            if name not in traces[idx]:
                raise KeyError(f"trace #{idx} 没有字段 {name!r}; "
                               f"可用: {sorted(traces[idx])}")
            if np.ndim(traces[idx][name]) != 0:
                raise ValueError(
                    f"scalars 里的 {name!r} 不是单个数字 "
                    f"(shape={np.shape(traces[idx][name])}); "
                    f"逐帧数组请放进 groups / fields")
    if bands is not None:
        for idx in indices:
            if bands not in traces[idx]:
                raise KeyError(f"trace #{idx} 没有字段 {bands!r}; "
                               f"可用: {sorted(traces[idx])}")
            if np.ndim(traces[idx][bands]) != 1:
                raise ValueError(f"bands 里的 {bands!r} 必须是逐帧数组, "
                                 f"得到 shape={np.shape(traces[idx][bands])}")

    nrows, ncols = len(indices), len(specs)
    if figsize is None:
        figsize = (4.2 * ncols + 0.8, 1.7 * nrows + 0.9)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False,
                             sharex="col")
    from matplotlib.lines import Line2D
    for r, idx in enumerate(indices):
        trace = traces[idx]
        xv = np.asarray(trace["frames"]) if x is None else np.asarray(x)
        # 标量: 该行每个子图的图例都带一份, 单看任一张图都能自解释
        scalar_text = [_format_scalar(n, trace[n]) for n in scalar_names]
        scalar_handles = [Line2D([], [], linestyle="none", marker="",
                                 label=text) for text in scalar_text]
        for c, (group, styles) in enumerate(zip(specs, group_styles)):
            ax = axes[r][c]
            if bands is not None:
                _draw_bands(ax, xv, np.asarray(trace[bands]),
                            bleach_frame=int(trace.get("bleach_frame", -1)))
            for name, (col, ls) in zip(group["fields"], styles):
                ax.plot(xv, np.asarray(trace[name], float), lw=linewidth,
                        color=col, linestyle=ls, label=_field_label(name))
            ax.set_title(f"#{trace.get('emitter_id', idx)}  {group['name']}",
                         fontsize=9)
            ax.tick_params(labelsize=7)
            if len(group["fields"]) == 1:
                ax.set_ylabel(_field_label(group["fields"][0]), fontsize=8)
            else:
                units = {_FIELD_UNITS.get(n) for n in group["fields"]}
                ax.set_ylabel(group["name"] if len(units) > 1 or None in units
                              else units.pop(), fontsize=8)
            handles, labels = ax.get_legend_handles_labels()
            if not (legend and len(group["fields"]) > 1):
                # 单字段时字段名已经在 y 轴标签里, 不重复; 标量仍要显示
                handles, labels = [], []
            handles = handles + scalar_handles
            labels = labels + scalar_text
            if handles:
                ax.legend(handles, labels, fontsize=7, frameon=False,
                          loc="best")
            if yscale:
                ax.set_yscale(yscale)
            if ylim_bottom is not None and not yscale:
                _clamp_bottom(ax, ylim_bottom)
            if r == nrows - 1:
                ax.set_xlabel(xlabel)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
#  群体统计
# ---------------------------------------------------------------------------
def _hist_int(ax, values, color, xlabel):
    """整数分布用整数 bin, 避免 0.5 处的空档。"""
    lo, hi = int(values.min()), int(values.max())
    ax.hist(values, bins=np.arange(lo - 0.5, hi + 1.5, 1.0),
            color=color, edgecolor="white", alpha=0.85)
    if hi - lo <= 25:
        ax.set_xticks(range(lo, hi + 1))
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("traces", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, axis="y", alpha=0.25)


def plot_state_summary(traces, *, bins="auto", figsize=(9.0, 6.2), title=None):
    """群体统计: 闪烁次数 / 亮度等级数 / ON 占比 / 漂白比例。

    只读 ``ThresholdStates`` 与 ``ChangepointStates`` **共用**的那几个字段, 所以两条
    路线的结果可以直接对比着看::

        plot_state_summary(threshold_traces)       # 阈值法
        plot_state_summary(changepoint_traces)     # 变点法

    Parameters
    ----------
    traces : Trace | dict | sequence
        trace 列表 (需要先跑过状态判定)。
    bins : int | str
        ``on_fraction`` 直方图的分箱; 离散量 (``n_on`` / ``n_states``) 自动用整数分箱。
    figsize, title
        图形尺寸与总标题。

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    traces = _one_trace(traces)
    if not traces:
        raise ValueError("traces 是空的")
    need = ("n_on", "n_switches", "n_states", "on_fraction", "bleach_frame")
    missing = [k for k in need if k not in traces[0]]
    if missing:
        raise KeyError(f"trace 缺少字段 {missing}; "
                       "先跑 ThresholdStates 或 ChangepointStates")

    n_on = np.array([t["n_on"] for t in traces], float)
    n_states = np.array([t["n_states"] for t in traces], int)
    on_frac = np.array([t["on_fraction"] for t in traces], float)
    bleached = np.array([t["bleach_frame"] for t in traces]) >= 0

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    ax_n, ax_s, ax_f, ax_b = axes.ravel()

    _hist_int(ax_n, n_on, "tab:blue", "ON episodes (n_on)")
    ax_n.set_title(f"Blinks  (median {np.median(n_on):.0f})", fontsize=10)

    _hist_int(ax_s, n_states, "tab:purple", "brightness levels (n_states)")
    ax_s.set_title(f"Levels  (max {int(n_states.max())})", fontsize=10)

    ax_f.hist(on_frac, bins=bins, color="tab:green", edgecolor="white",
              alpha=0.85)
    ax_f.axvline(np.median(on_frac), color="crimson", lw=1.4, ls="--",
                 label=f"median = {np.median(on_frac):.2f}")
    ax_f.legend(fontsize=8, frameon=False)
    ax_f.set_xlabel("ON fraction (alive frames)", fontsize=8)
    ax_f.set_ylabel("traces", fontsize=8)
    ax_f.set_title("Duty fraction", fontsize=10)
    ax_f.tick_params(labelsize=7)

    n_bl = int(bleached.sum())
    n_no = int(bleached.size - n_bl)
    ax_b.pie([n_bl, n_no],
             labels=[f"bleached ({n_bl})", f"not ({n_no})"],
             colors=["tab:red", "0.85"], autopct="%1.0f%%", startangle=90,
             wedgeprops=dict(edgecolor="white"))
    ax_b.set_title(f"Bleaching  (n={bleached.size})", fontsize=10)

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig
