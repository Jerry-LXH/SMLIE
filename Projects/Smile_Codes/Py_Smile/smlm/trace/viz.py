import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def plot_state_result(res, show_bgsub=False, ax=None, title=None):
    """
    单条 trace 的状态分析可视化。数据全部取自 StateResult 自身。
    上: 信号 + 阶梯拟合 + 状态色带 + 漂白线
    下: 状态序列 (阶梯)
    """
    frames = np.asarray(res.meta.get("frames", np.arange(res.trace.size)))
    y = res.trace_bgsub if show_bgsub else res.trace
    fit = res.fit_trace_bgsub if show_bgsub else res.fit_trace

    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(11, 5), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
    ax_sig, ax_seq = ax

    # 状态色带
    n = res.n_states
    colors = cm.viridis(np.linspace(0.15, 0.9, max(n, 1)))
    for j in range(1, n + 1):
        m = res.sequence == j
        if m.any():
            ax_sig.fill_between(frames, y.min(), y.max(), where=m,
                                color=colors[j - 1], alpha=0.12, step="mid")

    ax_sig.plot(frames, y, color="0.35", lw=0.8, label="signal")
    ax_sig.plot(frames, fit, color="crimson", lw=1.6, label="fit")

    # 状态均值虚线
    for info in res.state_info:
        lvl = info.mean_intensity - (res.background if show_bgsub else 0)
        ax_sig.axhline(lvl, color=colors[info.label - 1], ls=":", lw=0.8, alpha=0.7)

    # 转变点 & 漂白线
    for f, _ in res.transitions:
        ax_sig.axvline(frames[f], color="0.6", lw=0.5, ls="--", alpha=0.5)
    if res.bleach_frame is not None:
        ax_sig.axvline(frames[res.bleach_frame], color="black", lw=1.4,
                       label=f"bleach @ {frames[res.bleach_frame]}")

    ax_sig.set_ylabel("Intensity" + (" (bg-sub)" if show_bgsub else ""))
    ax_sig.legend(loc="upper right", fontsize=8, ncol=3)
    ttl = title or f"Emitter #{res.emitter_id}  |  {res.method}  |  " \
                   f"{res.n_states} states  |  lifetime={res.lifetime}"
    ax_sig.set_title(ttl, fontsize=10)

    # 状态序列
    ax_seq.step(frames, res.sequence, where="mid", color="tab:blue", lw=1.2)
    ax_seq.fill_between(frames, 0, res.sequence, step="mid", alpha=0.25)
    ax_seq.set_ylabel("state")
    ax_seq.set_xlabel("Frame")
    ax_seq.set_yticks(range(0, n + 1))
    ax_seq.set_ylim(-0.3, n + 0.3)

    plt.tight_layout()
    return ax


def plot_state_gallery(results, n_max=12, ncols=2, show_bgsub=False):
    """一次画多条 trace 的缩略图。数据全部取自 StateResult 自身。"""
    k = min(len(results), n_max)
    nrows = int(np.ceil(k / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6.5 * ncols, 2.2 * nrows),
                             squeeze=False)
    for i in range(k):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        res = results[i]                      # 不再需要 traces[i]
        frames = np.asarray(res.meta.get("frames", np.arange(res.trace.size)))
        y   = res.trace_bgsub    if show_bgsub else res.trace
        fit = res.fit_trace_bgsub if show_bgsub else res.fit_trace
        ax.plot(frames, y,   color="0.4",     lw=0.6)
        ax.plot(frames, fit, color="crimson", lw=1.2)
        if res.bleach_frame is not None:
            ax.axvline(frames[res.bleach_frame], color="k", lw=1.0)
        ax.set_title(f"#{res.emitter_id}  S={res.n_states}  τ={res.lifetime}",
                     fontsize=8)
    # 关掉多余的空子图
    for i in range(k, nrows * ncols):
        r, c = divmod(i, ncols)
        axes[r][c].axis("off")
    plt.tight_layout()
    return fig


def summarize_results(results):
    """返回汇总 dict, 便于统计漂白/态数分布。"""
    return dict(
        n_traces=len(results),
        n_states=np.array([r.n_states for r in results]),
        lifetimes=np.array([r.lifetime if r.lifetime is not None else np.nan
                            for r in results]),
        bleached=np.array([r.bleach_frame is not None for r in results]),
    )

def plot_summary(results, bins=20):
    """把 summarize_results 的统计信息画成 2x2 面板。"""
    s = summarize_results(results)
    n_states  = s["n_states"]
    lifetimes = s["lifetimes"]
    bleached  = s["bleached"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    ax_ns, ax_lt, ax_bl, ax_sc = axes.ravel()

    # (1) 态数分布 —— 离散, 用整数 bin
    if n_states.size:
        lo, hi = int(n_states.min()), int(n_states.max())
        edges = np.arange(lo - 0.5, hi + 1.5, 1.0)
        ax_ns.hist(n_states, bins=edges, color="tab:blue",
                   edgecolor="white", alpha=0.85)
        ax_ns.set_xticks(range(lo, hi + 1))
    ax_ns.set_title(f"State count  (n={s['n_traces']})", fontsize=10)
    ax_ns.set_xlabel("number of states")
    ax_ns.set_ylabel("count")

    # (2) 寿命分布 —— 忽略 NaN (未漂白/未定义)
    lt_valid = lifetimes[~np.isnan(lifetimes)]
    if lt_valid.size:
        ax_lt.hist(lt_valid, bins=bins, color="tab:green",
                   edgecolor="white", alpha=0.85)
        med = np.median(lt_valid)
        ax_lt.axvline(med, color="crimson", lw=1.4, ls="--",
                      label=f"median = {med:.0f}")
        ax_lt.legend(fontsize=8)
    ax_lt.set_title(f"Lifetime  (valid n={lt_valid.size})", fontsize=10)
    ax_lt.set_xlabel("lifetime (frames)")
    ax_lt.set_ylabel("count")

    # (3) 漂白比例 —— 饼图
    n_bl = int(bleached.sum())
    n_no = int(bleached.size - n_bl)
    if bleached.size:
        ax_bl.pie([n_bl, n_no],
                  labels=[f"bleached ({n_bl})", f"not ({n_no})"],
                  colors=["tab:red", "0.8"],
                  autopct="%1.0f%%", startangle=90,
                  wedgeprops=dict(edgecolor="white"))
    ax_bl.set_title("Bleaching fraction", fontsize=10)

    # (4) 态数 vs 寿命 —— 散点, 看两者是否相关
    mask = ~np.isnan(lifetimes)
    if mask.any():
        # 加一点水平抖动, 避免整数态数重叠成竖线
        jitter = (np.random.rand(mask.sum()) - 0.5) * 0.15
        ax_sc.scatter(n_states[mask] + jitter, lifetimes[mask],
                      s=18, alpha=0.5, color="tab:purple")
    ax_sc.set_title("States vs lifetime", fontsize=10)
    ax_sc.set_xlabel("number of states")
    ax_sc.set_ylabel("lifetime (frames)")

    fig.suptitle("Population summary", fontsize=12)
    plt.tight_layout()
    return fig