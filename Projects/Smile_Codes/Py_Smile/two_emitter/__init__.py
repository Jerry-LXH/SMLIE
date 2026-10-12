"""two_emitter —— 探究"亚衍射极限双发射体"问题的临时工具包。

问题
----
两个间距 ``d`` **低于衍射极限**的 **静态**发射体:

- 单帧里两个 PSF 重叠, 无法（或很难）直接分开;
- 但它们各自有**独立的 ON/OFF 闪烁**, 不同帧会"轮流单亮"。

那么: 利用闪烁的独立性恢复位置 / 间距, 能否比时间积分、短窗口等
常规方法更准确? 准确到什么程度?

先做 **level0**: 假设每帧"谁在发光"完全已知, 在给定单帧光子数与 PSF 下,
两个分子最多能靠多近还能被分开? (答案: 远低于衍射极限)

模块
----
- :mod:`two_emitter.dataset` : 数据集统一 schema (保存 / 加载) + ROI 等几何辅助
- :mod:`two_emitter.level0`  : CRLB 理论极限 + oracle(已知态) 估计器
- :mod:`two_emitter.level1`  : 状态未知时的 4 态联合 HMM (前向-后向 + EM), 含帧内占空建模
- :mod:`two_emitter.metrics` : 位置 / 间距误差指标
- :mod:`two_emitter.viz`     : 显示辅助

数据集默认落在 ``E:\\CodingProj\\Data\\Simulated_Datasets\\two_emitter_level0``;
由 ``Simulation_Smile/two_emitter_gen.ipynb`` 生成。
"""
from .dataset import (SCHEMA, DATASET_ROOT, ON_STATE, DEFAULT_HALF_PX,
                      TwoEmitterDataset, save_two_emitter_dataset,
                      list_datasets, load_all, write_manifest, load_manifest)
from . import dataset, level0, level1, metrics, viz

__all__ = [
    "SCHEMA", "DATASET_ROOT", "ON_STATE", "DEFAULT_HALF_PX",
    "TwoEmitterDataset", "save_two_emitter_dataset",
    "list_datasets", "load_all", "write_manifest", "load_manifest",
    "dataset", "level0", "level1", "metrics", "viz",
]

__version__ = "0.1.0"
