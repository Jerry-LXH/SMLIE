# Py_Smile (smlm) 关键数据类型与数据流说明

> 本文档总结 `Py_Smile/smlm` 包中定义在 `datatypes.py` 的关键数据类型，以及它们在整个单分子定位 / 时间轨迹分析流程中的数据流动方式。
> 所有数据类型的定义见 `smlm/datatypes.py`，处理管线各阶段见 `smlm/{io,camera,imgprocess,localize,drift,emitter,trace}`，端到端调用示例见 `test_proess.ipynb`。
> **注意**: 图像处理类 Transform (预处理/插值/去噪/反卷积/SOFI/RIED/漂移校正) 统一收敛到
> `smlm/imgprocess`；旧的 `smlm/preprocess` 与 `smlm/sofi` 目录已删除。
> `smlm/drift` 只负责漂移**估计**，校正侧 (`DriftCorrector`) 在 `smlm/imgprocess/temporal.py`。

---

## 1. 关键数据类型一览（定义于 `smlm/datatypes.py`）

| 类型 | 产生阶段 | 一句话说明 |
|---|---|---|
| `RawDataset` | IO 读取 | 所有 loader 的统一输出：原始 ADU 图像 + 头信息 |
| `ImageStack` | 相机标定 | 光子单位图像栈，携带 `pixel_size` 等下标定信息 |
| `Candidates` | 检测 | 整数像素级候选点 |
| `Localizations` | 定位 | 亚像素定位结果，**列式存储**（每个数组长度 = 分子总数） |
| `DriftTrace` | 漂移估计 | 逐帧累积漂移曲线（估计与校正之间的契约） |
| `Emitter` | 轨迹连接 | 单个发光体的连接轨迹，通过 `loc_idx` 回查源定位 |
| `EmitterStats` | 统计 | 每个 emitter 的汇总统计量 |
| `Trace` | 轨迹提取 | 稠密强度时间轨迹（含 OFF 帧） |
| `StateInfo` | 轨迹分析 | 单个离散态的统计（均值 / 方差 / 占空比） |
| `StateResult` | 轨迹分析 | 整条 trace 的态序列分析结果（变点 / 漂白 / 寿命） |

### 1.1 各类型详细说明

**`RawDataset`** — 所有 loader（SIF / TIF）的统一输出。
- `data`: `(T, H, W)` float32 原始 ADU 计数（与 `ImageStack` 一致，帧在**第一维**）
- `exposure_time`, `em_gain`：相机头信息，供 EMCCD 转光子用
- `metadata`: dict（含 `pixel_size` 等）
- `n_frames`: 返回 `data.shape[0]`（帧数）

**`ImageStack`** — camera 层的输出约定（`Camera.process()` 统一包装）。
- `data`: `(T, H, W)` float32 **光子数**（帧在**第一维**）
- `pixel_size`: nm/pixel，定位转 nm 时使用
- 属性：`n_frames`, `shape`
- 坐标约定：`data[t, row(y), col(x)]`；`y` 对应行、`x` 对应列

**`Candidates`** — 检测输出，整数像素坐标。
- `frame`(int), `y`(int, row), `x`(int, col)

**`Localizations`** — 定位输出，SMLM 通用列式存储，**核心数据结构**。
- `frame`(int), `x`/`y`(亚像素, pixel), `photons`, `bg`(光子/pixel), `sigma`(PSF 宽度), `lp`(横向定位精度 nm)
- `z` / `z_lp`：轴向坐标与轴向定位精度（nm，来自 z 标定）；**2D 定位时为 None**
- `extra`: dict，存放拟合残差 / CRLB / 各向异性 σx、σy 等
- 支持 `locs[mask]` / `locs[idx]` / `locs[slice]` 索引（同步切所有列与 `extra` 中等长数组）
- `to_dataframe()` 转 pandas

**3D 定位方案（`smlm/localize/loc/` 包，按方案分文件）**
- `common.py`：共享的 ROI 提取 + 多进程 worker 池
- `fisher.py`：共享的 Fisher 矩阵机制 `fisher_from_jacobian`（各模型的 Jacobian 留在方案文件）
- `gaussian2d.py`：`MLELocalizer`（2D 各向同性高斯，z=None）
- `astigmatic.py`：`AstigmaticLocalizer`（3D 像散/柱面镜，椭圆高斯；σx/σy 存 `extra`，z 经 `z_calibration` 恢复，否则 None）
- `doublehelix.py`：`DoubleHelixLocalizer`（3D 双螺旋，预留骨架）

**`DriftTrace`** — 漂移估计与校正之间的契约。
- `dx`/`dy`: `(n_frames,)` 累积漂移，单位 pixel，索引 0-based 帧号；`dx`=列方向、`dy`=行方向
- 校正约定：帧 t 坐标**减去** `(dx[t], dy[t])`；图像按 `(-dy, -dx)` 反向平移
- `total_drift`, `to_dataframe()`, `DriftTrace.zeros(n_frames)`（零漂移占位）

**`Emitter`** — 轨迹连接输出（`loc_idx` 是关键设计）。
- `frame`, `x`, `y`, `loc_idx`（指向源 `Localizations` 的行索引，**避免冗余存储**）
- `bleached`: 是否片内漂白
- 派生属性：`on_frame` / `last_frame` / `n_locs` / `span` / `n_dark` / `duty_cycle`（闪烁、漂白、常亮统一用这些表达）
- 方法：`xy`, `centroid`, `covariance()`, `eccentricity()`

**`EmitterStats`** — `collect_statistics()` 输出（向量化，从源 locs 回查 photons/bg/sigma）。
- 每项长度 = emitter 数：`centroid`, `survival_sec`, `n_locs`, `duty_cycle`, `bleached`, `eccentricity`, `photons_per_sec`, `photons_sum`, `sigma_mean`, `bg_mean`

**`Trace`** — 强度轨迹提取输出。
- `emitter_id`, `frames`（稠密帧轴，含 OFF 帧）, `intensity`（**分析主对象**：扣背景后的信号）, `raw`, `bg`, `snr`, `x`/`y`（提取中心）, `meta`

**`StateInfo`** — `label`, `mean_intensity`, `std_intensity`, `occupancy`（相对 active 段的占空比）。

**`StateResult`** — 轨迹态分析输出。
- `sequence`: 逐帧状态标签（0 = OFF/背景）
- `transitions`: `(n,2)` = `[frame_index, new_state]`
- `fit_trace`, `trace_bgsub`, `fit_trace_bgsub`, `background`
- `bleach_frame` / `start_frame` / `lifetime`（漂白帧、起始帧、寿命，索引 0-based）
- `state_info`, `n_states`, `method`, `meta`

---

## 2. 端到端数据流

```mermaid
flowchart TD
    A[原始 SIF/TIF] -->|"load_raw()"| B[RawDataset<br/>ADU (T,H,W)]
    B -->|"Camera.process()"| C[ImageStack<br/>光子 (T,H,W) + pixel_size]
    C -->|"Pipeline(Crop, SelectFrames)"| D[ImageStack<br/>裁剪/选帧后]
    D -->|"LocalMaxDetector.detect()"| E[Candidates<br/>整数像素候选]
    E -->|"MLELocalizer.localize()"| F[Localizations<br/>亚像素+光子+sigma+lp]
    F -->|"filter_locs()"| G[Localizations<br/>物理量范围过滤]
    G -->|"CloudDriftEstimator.estimate()"| H[DriftTrace<br/>逐帧累积漂移]
    H -->|"imgprocess.DriftCorrector.apply()"| I[(corr_image: ImageStack<br/>corr_locs: Localizations)]
    I -->|"analyze_emitters()<br/>GapClosingLinker+Filter"| J[list[Emitter]<br/>loc_idx 回查源]
    J -->|"collect_statistics()"| K[EmitterStats]
    J -->|"ApertureTraceExtractor.extract()"| L[list[Trace]<br/>稠密强度序列]
    L -->|"analyze_traces()<br/>get_analyzer('changepoint')"| M[list[StateResult]<br/>态序列+漂白+寿命]
    I -.->|"也可以直接用 corr_locs"| J
```

### 2.1 各阶段与数据类型对应

| 阶段 | 模块 | 输入 → 输出 |
|---|---|---|
| 1. 读取 | `smlm/io` | 文件 → `RawDataset` |
| 2. 转光子 | `smlm/camera` | `RawDataset` → `ImageStack` |
| 3. 预处理/图像处理 | `smlm/imgprocess` | `ImageStack` → `ImageStack`（`dataclasses.replace` 拷贝；含 Crop/SelectFrames/傅里叶上采样/高斯去噪/RL3D/稀疏反卷积/OTF重加权/SOFI/RIED） |
| 4. 检测 | `smlm/localize` (detect) | `ImageStack` → `Candidates` |
| 5. 定位 | `smlm/localize` (loc) | `ImageStack` + `Candidates` → `Localizations` |
| 5'. 一体化+缓存 | `smlm/localize` (`ClassicalPipeline`) + `smlm/io/cache` (`run_cached`) | `ImageStack` → `Localizations`（带 `.locs.npz` 磁盘缓存） |
| 6. 过滤 | `smlm/localize` (filter) | `Localizations` → `Localizations` |
| 7. 漂移 | `smlm/drift` (估计) + `smlm/imgprocess` (校正) | `Localizations` → `DriftTrace`；再用 `imgprocess.DriftCorrector` 对 `ImageStack`+`Localizations` 校正 |
| 8. Emitter | `smlm/emitter` | `Localizations` → `list[Emitter]`（→ `EmitterStats`） |
| 9. 提取轨迹 | `smlm/trace` (extract) | `list[Emitter]` + 图像栈 → `list[Trace]` |
| 10. 分析轨迹 | `smlm/trace` (methods) | `list[Trace]` → `list[StateResult]` |

---

## 3. 关键约定（容易踩坑）

- **轴顺序已统一**：`RawDataset.data` 与 `ImageStack.data` 均为 `(T, H, W)`，帧在第一维。SIF 底层 `imageData` 是 `(H, W, T)`，由 `SIFLoader` 转成 `(T, H, W)`。
- **坐标语义**：全程 `x` = 列(col)、`y` = 行(row)；`Candidates` 为整数像素，`Localizations` 为 0-based 亚像素。
- **单位**：检测 / 定位 / 漂移在 **pixel** 下进行；`ImageStack.pixel_size`（nm/pixel）只在可视化或换算 nm 时使用。
- **`Emitter.loc_idx`** 是指回源 `Localizations` 的行索引，因此 `photons/sigma/bg` 等不回存到 Emitter，需用 `locs.photons[emitter.loc_idx]` 取回。
- **`DriftTrace`** 的 `dx`/`dy` 是"坐标要减去的量"，图像平移方向与之相反。
- **轨迹 vs 定位粒度**：`Localizations` 每个点一条记录；`Trace` 每条 emitter 一条稠密序列（含暗帧），两者通过 `emitter_id`/`loc_idx` 关联。
- **注册机制**：`smlm/trace` 的 `methods/__init__.py` 必须在导入时执行，`@register_analyzer` 装饰器才会把算法填入注册表，`get_analyzer()` 才能找到（当前内置 `changepoint`，HMM 等可扩展）。

---

## 4. 快速上手（对应 `test_proess.ipynb`）

```python
from smlm.io import load_raw
from smlm.camera import EMCCD
from smlm.imgprocess import Pipeline, Crop, SelectFrames
from smlm.localize import ClassicalPipeline, LocalMaxDetector, MLELocalizer, filter_locs
from smlm.io import run_cached
from smlm.drift import CloudDriftEstimator
from smlm.imgprocess import DriftCorrector      # 校正侧在 imgprocess
from smlm.emitter import analyze_emitters, collect_statistics
from smlm.trace import ApertureTraceExtractor, analyze_traces

raw     = load_raw("xxx.sif")                              # RawDataset
image   = EMCCD(conv_factor=5.75, em_gain=raw.em_gain).process(raw, pixel_size=106.6)  # ImageStack
stack   = Pipeline(Crop(y=(100,400), x=(100,400)), SelectFrames(stop=500))(image)      # ImageStack
method  = ClassicalPipeline(LocalMaxDetector(...), MLELocalizer(...))
locs    = run_cached(method, stack, src)                   # Localizations（带缓存）
locs    = filter_locs(locs)                                # Localizations
drift   = CloudDriftEstimator(...).estimate(locs)          # DriftTrace
corr_image, corr_locs, keep = DriftCorrector(crop=(257,257)).apply(drift, image=stack, locs=locs)
emitters, _ = analyze_emitters(corr_locs, n_frames=None, ...)   # list[Emitter]
stats   = collect_statistics(emitters, locs, raw.exposure_time, 0)  # EmitterStats
traces  = ApertureTraceExtractor(edge=3).extract(emitters, corr_image.data, corr_image.n_frames)  # list[Trace]
results = analyze_traces(traces, method="changepoint", penalty=3, bleach_tail=50)       # list[StateResult]
```
