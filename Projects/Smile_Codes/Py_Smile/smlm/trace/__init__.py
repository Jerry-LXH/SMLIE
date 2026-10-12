# --- 数据类型 ---
from ..datatypes import Trace, StateInfo, StateResult

# --- 框架核心 (工厂 / 顶层函数) ---
from .base import (
    StateAnalyzer,
    get_analyzer,
    analyze_traces,
    register_analyzer,
)

from .extract import ApertureTraceExtractor

# --- 触发内置算法注册 (关键!) ---
# import 这一句会执行 methods/__init__.py, 里面的 @register_analyzer 随之运行,
# 各算法才会进入注册表, get_analyzer 才能找到它们。
from . import methods   # noqa: F401

# --- 可视化 ---
from .viz import plot_state_result, plot_state_gallery, summarize_results, plot_summary

# 可选: 明确对外暴露的接口
__all__ = [
    "Trace", "StateInfo", "StateResult",
    "StateAnalyzer", "get_analyzer", "analyze_traces", "register_analyzer",
    "plot_state_result", "plot_state_gallery", "summarize_results",
]