# Py_Smile (smlm) 关键数据类型与数据流说明

> 本文档总结 `Py_Smile/smlm` 包中定义在 `datatypes.py` 的关键数据类型，以及它们在整个单分子定位 / 时间轨迹分析流程中的数据流动方式。
> 所有数据类型的定义见 `smlm/datatypes.py`；**各子包/模块的分工见第 3 节**；端到端调用示例见
> `smTest.ipynb`（单通道）与 `smTest_dual.ipynb`（双通道 FRET）。
> **注意**: 图像处理类 Transform (预处理/插值/去噪/反卷积/SOFI/RIED/漂移校正) 统一收敛到
> `smlm/imgprocess`；旧的 `smlm/preprocess` 与 `smlm/sofi` 目录已删除。
> `smlm/drift` 只负责漂移**估计**，校正侧 (`DriftCorrector`) 在 `smlm/imgprocess/temporal.py`。

---

## 1. 关键数据类型一览（定义于 `smlm/datatypes.py`）

| 类型 | 产生阶段 | 一句话说明 |
|---|---|---|
| `RawDataset` | IO 读取 | 所有 loader 的统一输出：原始 ADU 图像 + 头信息 |
| `ImageStack` | 相机标定 | 光子单位图像栈，携带 `pixel_size` 等下标定信息 |
| `MultiImageStack` | 多通道取框 | 同形的多通道图像栈 + 通道名/偏移/坐标配准 |
| `Candidates` | 检测 | 整数像素级候选点 |
| `Localizations` | 定位 | 亚像素定位结果，**列式存储**（每个数组长度 = 分子总数） |
| `DriftTrace` | 漂移估计 | 逐帧累积漂移曲线（估计与校正之间的契约） |
| `Emitter` | 轨迹连接 | 单个发光体的连接轨迹，通过 `loc_idx` 回查源定位 |
| `EmitterStats` | 统计 | 每个 emitter 的汇总统计量 |
| `Trace` | 轨迹提取 | 稠密强度时间轨迹（含 OFF 帧）+ 状态字段（见 1.1） |

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

**`MultiImageStack`** — 多通道载体（双视角取框后的统一容器）。各通道必须**同形**，
因为下游是逐像素/逐坐标对齐处理，不做形变重采样。
- `data`: `(C, T, H, W)`；`channels`: 通道名列表；`reference`: 参考通道下标
- `offsets`: `(C, 2)` 各通道在探测器上的原点 `(y0, x0)`（整数裁切时吸收掉粗偏移）
- `reg`: `{通道: DualViewCalibration}`，参考通道不在其中；`geometry()` 据此构造
  `multiplex.DualViewGeometry`，做"通道 local ↔ 参考通道坐标"的换算
- `metadata`: 取框来源（`analysis_crops` / `roi` / `margin`）等
- 方法/属性：`n_channels`, `n_frames`, `shape`, `ref_name`, `channel(name)`,
  `channel_offset(name)`, `crop(y=..., x=...)`, `subset(channels)`,
  `to_reference(xy, channel)`, `with_data(data)`

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

**`Trace`** — 强度轨迹（单通道 / 双通道共用）。字段很多，分两类：

- **取亮度**：`frames`（稠密帧轴，含 OFF 帧）、`raw`、`bg`、`bg_std`、`bg_flux`、
  `n_bg_pixels`、`neighbours`、`intensity`（**分析主对象**：扣背景后的信号）、`snr`；
  外加几何/口径元信息 `x`/`y`、`channel`、`radius`、`mode`、`sw`/`sw2`、`area`、`sigma`
- **状态**：由 `ThresholdStates` / `ChangepointStates` 这类 filler 直接写回同一条 trace
  （**没有**单独的 `StateResult` 类型了）：`state`（逐帧等级）、`fit`、`n_states`、
  `bleach_frame`、`alive`、`n_on`、`n_switches`、`n_cp`、`on_frames`、`on_fraction`
  
  完整字段表见 `smlm.trace.TRACE_FIELDS`（状态部分是 `STATE_FIELDS`）；
  `to_table(traces)` 一次性导出成 pandas，`select_traces(traces, n_on=..., ...)` 按这些
  标量筛分子。

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
    I -->|"analyze_emitters_sm()<br/>GapClosingLinker + 合并 + 过滤"| J[list[Emitter]<br/>loc_idx 回查源]
    J -->|"collect_statistics()"| K[EmitterStats]
    J -->|"TraceFillPipeline().run(fillers)"| L[list[Trace]<br/>raw/bg/intensity/snr]
    L -->|"ThresholdStates / ChangepointStates"| M[Trace 状态字段<br/>state/n_on/bleach_frame]
    M -->|"select_traces(n_on=..., ...)"| N[筛过的 list[Trace]]
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
| 9. 提取轨迹 | `smlm/trace` (`fill` + `signal`/`annulus`/`background`/`calculate`) | `list[Emitter]` + 图像栈 → `list[Trace]` |
| 10. 分析轨迹 | `smlm/trace` (`states`) | `list[Trace]` → 同一条 trace 上补 `state`/`n_on`/`bleach_frame` |

双通道（双视角）流程在第 8~10 步之外多一层：`smlm/multiplex` 负责取框、通道间配准、
多通道 emitter 聚类与多通道 trace 提取；`smlm/fret` 在其上做"红/绿 trace 对 → 事件 →
γ/α/β → E"。见第 3 节两个子包的分工。

---

## 3. 子模块一览

```
smlm/
├── datatypes.py      全局数据类型（第 1 节）
├── viz.py            单通道阶段的通用绘图
├── io/               读文件 / 定位缓存 / 导出
├── camera/           ADU → photon
├── imgprocess/       ImageStack → ImageStack 的 Transform（预处理/去噪/反卷积/SOFI/漂移校正）
├── localize/         检测 + 定位（含 3D 方案 loc/）
├── drift/            漂移**估计**
├── emitter/          定位 → 轨迹（连接/合并/过滤/统计）
├── trace/            取亮度 + 状态判定（filler 流水线）
├── multiplex/        多通道（双视角）配准与多通道 trace
└── fret/             双通道 FRET：事件 / 筛选 / γαβ / E
```

### 3.1 `smlm/datatypes.py`、`smlm/viz.py`（顶层）

| 文件 | 内容 |
|---|---|
| `datatypes.py` | 第 1 节那 9 个类型；全包共用，不放任何算法 |
| `viz.py` | 单通道阶段的通用图：`show_frames`、`average_frames`、`overlay_points`、`plot_drift`、`plot_distributions`、`plot_trajectories`、`plot_emitter_stats`、`plot_bleaching_summary`、`emitter_roi`、`show_emitter_movies`、`preview_emitters`、`save_emitters`、`check_trace` |

### 3.2 `smlm/io` — 读取 / 缓存 / 导出

| 文件 | 内容 |
|---|---|
| `__init__.py` | `load_raw(path, **kw)`：按扩展名分发到 SIF/TIF，返回 `RawDataset` |
| `base.py` | `Loader` 抽象基类 |
| `sif.py` / `_readsif.py` | `SIFLoader`：Andor SIF 解析；`_readsif` 只管容器格式，不解释头信息语义 |
| `tif.py` | `TIFLoader`、`save_tif` |
| `cache.py` | 定位结果的**磁盘缓存**：`run_cached`（跑+存+读）、`save_cached_locs`、`load_locs`、`save_locs`、`result_path`；指纹 `data_fingerprint`/`image_signature`/`config_hash` 决定缓存是否命中 |
| `video.py` | `save_movie`、`preview_movie`、`auto_clim`：栈（可叠加 trace）导出 MP4，帧率默认 `1/exposure_time`；灰度量化走实测 LUT，逐字节可复现 |

### 3.3 `smlm/camera` — ADU → photon

| 文件 | 内容 |
|---|---|
| `base.py` | `Camera.process(raw)`：统一出口，返回 `ImageStack`（光子单位 + `pixel_size`） |
| `emccd.py` | `EMCCD(conv_factor, em_gain)`：含 EM 增益的 ADU → photon |
| `scmos.py` | `sCMOS(gain_map, offset_map)`：逐像素 gain/offset 校正 |

### 3.4 `smlm/imgprocess` — 图像 → 图像

所有 `Transform` 都是 `ImageStack` 进、`ImageStack` 出（`dataclasses.replace` 语义），可以用 `Pipeline` 串起来。

| 文件 | 内容 |
|---|---|
| `base.py` | `Transform` / `Pipeline` 抽象 |
| `preprocess.py` | `Crop`（空间裁切）、`SelectFrames`（时间选帧） |
| `temporal.py` | 时间维处理：`BleachCorrect`（漂白校正）、**`DriftCorrector`（漂移校正侧）**、`TemporalMean` |
| `interp.py` | 傅里叶上/下采样（`FourierUpsample`、`fourier_*`） |
| `denoise.py` | `GaussianFilter` |
| `dehotpixel.py` | `Dehotpixel`：sCMOS 热点像素 |
| `deconv.py` | `RL3D`（加速 Richardson-Lucy）、`SparseDeconv`、`OtfReweight` |
| `sofi.py` | `Sofit` / `CumulantSofit` / `RiedSofit`：SOFI / RIED 超分辨 |
| `core.py` | SOFI/RIED 的纯数学（与 Transform 解耦，便于单测） |
| `viz.py` | SOFI/RIED 结果对比显示 |

### 3.5 `smlm/localize` — 检测 + 定位

| 文件 | 内容 |
|---|---|
| `base.py` | `Detector` / `Localizer` / `LocalizationMethod` / `ClassicalPipeline`（检测+定位合成一个可缓存的方法） |
| `detect.py` | `LocalMaxDetector`：k·σ 阈值 + 局部极大值 |
| `filter.py` | `filter_locs`：光子 / σ / 背景等物理量范围过滤 |
| `pipeline.py` | `localize_stack(...)`：一体化 + 磁盘缓存（内部转调 `io.cache`） |
| `loc/common.py` | ROI 提取 + 多进程 worker 池 |
| `loc/fisher.py` | `fisher_from_jacobian`：各模型共享的 Fisher 信息机制 |
| `loc/gaussian2d.py` | `MLELocalizer`：2D 各向同性高斯（`z=None`） |
| `loc/astigmatic.py` | `AstigmaticLocalizer`：3D 像散（σx/σy 存 `extra`，z 靠标定恢复） |
| `loc/doublehelix.py` | `DoubleHelixLocalizer`：3D 双螺旋（骨架） |

### 3.6 `smlm/drift` — 只负责**估计**

| 文件 | 内容 |
|---|---|
| `base.py` | `DriftEstimator` 抽象 |
| `cloud.py` | `CloudDriftEstimator`：把定位点云分段求平均图、互相关，估出逐帧累积漂移 → `DriftTrace` |

校正侧在 `smlm/imgprocess/temporal.py` 的 `DriftCorrector`：帧 t 的坐标**减去** `(dx[t], dy[t])`，图像按 `(-dy, -dx)` 反向平移。

### 3.7 `smlm/emitter` — 定位 → 轨迹

| 文件 | 内容 |
|---|---|
| `linking.py` | `Linker` / `GapClosingLinker`：按帧距把定位连成轨迹（允许暗隙 `max_gap`） |
| `merging.py` | `EmitterMerger`：把定位器拆成两条的同一个分子合回去 |
| `filtering.py` | `MinDurationFilter` / `FirstFrameFilter` / `EccentricityFilter` / `FilterPipeline` / `ecc_pvalue` |
| `statistics.py` | `collect_statistics(emitters, locs, ...)` → `EmitterStats` |
| `pipeline.py` | `analyze_emitters_sm(...)`：连接 → 合并 → 过滤 一条龙，返回 `(emitters, emitters_raw)` |
| `base.py` | `_TrackBuilder`：连接器共用的骨架 |

### 3.8 `smlm/trace` — 取亮度 + 状态判定

核心抽象是 **filler 流水线**：每个 `TraceFiller` 只往 trace 上填一个/一组字段，声明 `requires` / `provides`，由 `TraceFillPipeline` 排序后执行；所以"换信号算法"和"换背景算法"互不影响。

| 文件 | 内容 |
|---|---|
| `fill.py` | `TraceBox` / `FillContext` / `TraceFiller` / `TraceFillPipeline` / `register_field` / `to_table` / `as_centers_yx` |
| `signal.py` | `ApertureRaw`：信号阶段（孔径内原始积分，`raw` 口径） |
| `annulus.py` | `AnnulusRing`：背景阶段（分子周围环形背景逐帧估计） |
| `aperture.py` | 孔径积分工具：`ApertureSpec`、`Aperture`、`aperture_kernel`、`radius_from_sigma`、`sigma_per_molecule`、`crop_centered`、`clean_sigma` |
| `background.py` | 背景纯函数：`AnnulusBackground`、`ring_mask`、`clip_median(_std)`、`valid_mask`、`ring_budget`、`moving_average`、`model_weights`、`emitter_frame_photons`、`emitter_sigma` |
| `calculate.py` | 只做算术的两个阶段：`IntensityFromRawBg`（`intensity = raw - bg·sw/sw2`）、`PoissonSNR`；`standard_fillers` 是默认组合 |
| `states.py` | 状态判定：`ThresholdStates`（阈值二值化）、`ChangepointStates`（变点分段）、`select_traces`、`summarize_states`、`runs_of` |
| `changepoint.py` | 变点检测的纯算法（只依赖 numpy）：`segment`、`binary_seg`、`refine_cp`、`agglomerative`、`forward`、`robust_noise` |
| `viz.py` | `plot_trace_fields`（多字段叠加）、`plot_state_summary`、`show_trace_movies` |

### 3.9 `smlm/multiplex` — 多通道（双视角）

分工原则：**图像只做漂移平移**（不重采样，避免引入逐分子 κ 抖动）；**坐标层**才做完整配准。

| 文件 | 内容 |
|---|---|
| `roi.py` | `MultichannelROI`：通道 0 的框 + 各通道 `(dy, dx)` 偏移，定义"同一物理窗口"；`show_roi(stack)` 出确认图，`windows` 给出各通道窗口，`shape` 给统一尺寸 |
| `bundle.py` | `MultiChannelCalibrationBundle`：多通道标定的落盘/载入（npz + manifest + `BUNDLE_FORMAT` 版本） |
| `transform.py` | `DualViewCalibration`：单通道对的 2D 配准模型（poly/affine），`transform()` / `apply_locs()` |
| `coords.py` | `DualViewGeometry`：多通道 local ↔ 参考通道坐标、公共可用区；`map_channel_xy`、`calib_windows_from_calibs` |
| `fit.py` | 从配对点稳健拟合：`fit_transform`（4·MAD）、`cv_select_order`、`poly_terms`、`design`、`n_params` |
| `pair.py` | `PairSet` / `pair_channels` / `pair_multichannels` / `merge_pairs`：一对一配对（三条约束）与多文件合并 |
| `pipeline.py` | `pair_multichannel_files`（逐文件配对）、`fit_multichannel_pairs`（合并拟合）、`multichannel_stack_from_calibration`（按 bundle 取框）、`localize_channel`、`mean_frame` |
| `drift.py` | 多通道漂移：`CropWindow`、`estimate_drift_multi`（一次估计、按几何换算到各通道）、`correct_multi`、`drift_crop_window`、`drift_requirements`、`merge_registered` |
| `emitters.py` | `analyze_multichannel_emitters`：**逐通道各自聚类**（用本通道自己的定位），`channel_locs` / `channel_names` / `emitter_centers_xy` / `psf_sigma_from_locs`（逐通道拟合宽度中位 → 取 trace 用的 σ） |
| `traces.py` | `extract_multichannel_traces`：每通道可指定"用哪个通道的 emitter 当中心"和"用哪个通道的背景屏蔽"，输出 `{通道: list[Trace]}`。`sigma` 收 `{通道: px}`，**用 `psf_sigma_from_locs` 的结果**：`0.2165·λ/NA` 的理论值在有轻微离焦时明显偏小（实测 0.93 vs 1.43 px），孔径半径被夹到 2 px、匹配滤波幅度衰减到 ~0.6，所有绝对门限都跟着失真 |
| `viz.py` | `plot_multichannel_crop`、`plot_multichannel_roi`、`plot_multichannel_localizations`、`plot_multi_drift`、`plot_drift_fov`、`plot_channel_traces`、`plot_displacement_field`、`plot_residual_hist`、`plot_residual_map`、`plot_overview` |

### 3.10 `smlm/fret` — 双通道 FRET

约定：`trace_a` 永远是**受体（红）**，`trace_d` 是**供体（绿）**。
一次完整分析（见 `smTest_dual.ipynb`）走：`select_pairs` → `calibrate` →
`fret_efficiency` → `save_bundle`；落盘目录是 `<sif 主名>_fret/`，五个文件名固定在
`fret.FILES` 里，下游按名字找。

| 文件 | 内容 |
|---|---|
| `events.py` | `EventSpec` / `find_step` / `scan_steps` / `post_level` / `red_state`：受体侧把**每个亮段的结束**当候选台阶（迟滞二值化，闪烁的暗隙不再冒充漂白），供体侧用亮帧中位；带 7 条物理门限。事件上附 `red_state`（平稳/闪烁/闪光结束/没真灭）与 `red_on`（受体逐帧 ON 掩码） |
| `pairs.py` | `PairCriteria` / `select_pairs`：漏斗式筛选（event → red → duty → flags → spread → attr → neigh → after），逐级记录否决原因。`spread` **只在受体亮着的帧上**量；`blinking` 时 `after` 按 `blink_relax` 放宽 |
| `calibrate.py` | `calibrate` / `FretCalibration`（可落盘）/ `gamma_points` / `local_step` / `estimate_alpha` / `estimate_beta` / `estimate_gamma` / `robust_line` |
| `efficiency.py` | `EfficiencySpec` / `fret_efficiency` / `window_mean`：先算 `S_A`，再给 `E`（用 γ）与 `E_don`（不用 γ 的交叉检验）。`I_A^pre` / `I_D^pre` 只平均**受体亮着的帧**（受体闪时，熄灭间隙里供体已去淬灭，混进来会把 E 压低） |
| `datatypes.py` | `StepEvent` / `FretResult` / `CRITERIA` / `CRITERIA_STRICT` / `passed`：FRET 专用类型与判据名 |
| `io.py` | `save_csv`（逐分子） / `save_traces_csv`（全帧长表） / `load_traces_csv`（读回成按分子分块的数组） / `TRACE_COLS` |
| `export.py` | `save_bundle(out_dir, ...)`：一次把 `fret.csv` / `E_distribution.png` / `good_traces_all.png` / `traces_all.csv` / `good_traces.csv` 全部写好（`include_all=False` 就不写未筛选的那张表）；只做编排，写表在 `io`、画图在 `viz`。`FILES` 是对外承诺的文件名 |
| `viz.py` | `plot_pair` / `plot_pair_gallery`（含 `which='kept'/'rejected'`）/ `plot_gamma` / `plot_efficiency` / `plot_efficiency_vs_brightness` |

---

## 4. 关键约定（容易踩坑）

- **轴顺序已统一**：`RawDataset.data` 与 `ImageStack.data` 均为 `(T, H, W)`，帧在第一维。SIF 底层 `imageData` 是 `(H, W, T)`，由 `SIFLoader` 转成 `(T, H, W)`。
- **坐标语义**：全程 `x` = 列(col)、`y` = 行(row)；`Candidates` 为整数像素，`Localizations` 为 0-based 亚像素。
- **单位**：检测 / 定位 / 漂移在 **pixel** 下进行；`ImageStack.pixel_size`（nm/pixel）只在可视化或换算 nm 时使用。
- **`Emitter.loc_idx`** 是指回源 `Localizations` 的行索引，因此 `photons/sigma/bg` 等不回存到 Emitter，需用 `locs.photons[emitter.loc_idx]` 取回。
- **`DriftTrace`** 的 `dx`/`dy` 是"坐标要减去的量"，图像平移方向与之相反。
- **轨迹 vs 定位粒度**：`Localizations` 每个点一条记录；`Trace` 每条 emitter 一条稠密序列（含暗帧），两者通过 `emitter_id`/`loc_idx` 关联。
- **trace 状态是策略类**：`smlm/trace/states.py` 的 `ThresholdStates` / `ChangepointStates` 都是"吃 `list[Trace]`、就地补 `state`/`n_on`/`bleach_frame` 字段"，没有注册表也不用装饰器；`select_traces(traces, n_on=..., n_switches=...)` 按统计量筛。

---

## 5. 快速上手

对应 notebook：
- `smTest.ipynb` —— 单通道：读 sif → 定位 → 漂移 → emitter → trace 状态
- `smTest_dual.ipynb` —— 双通道 FRET：ROI → 漂移 → emitter → trace 对 → 判据 → γ/α/E
- `smTest_dual_calibration.ipynb` —— 双通道标定（`MultichannelROI` + `MultiChannelCalibrationBundle`）
- `smTest_dual_batch.ipynb` —— 批量处理（每个视野一键落盘到 `<sif 主名>_fret/`）

```python
from smlm.io import load_raw, run_cached
from smlm.camera import EMCCD
from smlm.imgprocess import Pipeline, Crop, SelectFrames, DriftCorrector
from smlm.localize import LocalMaxDetector, MLELocalizer, filter_locs, localize_stack
from smlm.drift import CloudDriftEstimator
from smlm.emitter import analyze_emitters_sm, collect_statistics
import smlm.trace as tr

raw   = load_raw("xxx.sif")                                          # RawDataset
image = EMCCD(conv_factor=5.75, em_gain=300).process(raw, pixel_size=106.6)   # ImageStack
stack = Pipeline(Crop(y=(0, 512), x=(0, 512)), SelectFrames())(image)         # ImageStack

locs = localize_stack(stack, k_sigma=2.6, det_kw=DET_KW, loc_kw=LOC_KW,
                      filt_kw=FILT_KW, src_path="xxx.sif", read_cache=True)   # Localizations
drift = CloudDriftEstimator(...).estimate(locs[locs.frame < stack.n_frames])  # DriftTrace
corr_image, corr_locs, keep = DriftCorrector(crop=(257, 257)).apply(drift, image=stack, locs=locs)
emitters, _ = analyze_emitters_sm(corr_locs, n_frames=None, search_radius=0.6,
                                  max_gap=1000, min_frames=10, do_merge=True)  # list[Emitter]
stats  = collect_statistics(emitters, locs, raw.exposure_time, 0)             # EmitterStats

# trace: fill pipeline 里叠一组"取亮度的方法"，状态类也是 filler
centers_yx = np.array([[e.centroid[1], e.centroid[0]] for e in emitters])   # [row, col]
fillers = [tr.ApertureRaw(mode='hard', radius=3), tr.AnnulusRing(gap=2, width=4, excl_radius=3),
           tr.IntensityFromRawBg(), tr.PoissonSNR(),
           tr.ThresholdStates()]              # 补 state / n_on / bleach_frame
traces  = tr.TraceFillPipeline(corr_image.data, centers_yx, sigma=ph, tag='r').run(fillers)
tr.select_traces(traces, n_on=(2, None))
```

`sigma=ph` 用**拟合宽度**（`np.median(corr_locs.sigma)`，多通道下用
`mp.psf_sigma_from_locs(multiLocs)`），不要用 `0.2165·λ/NA` 的理论值：实测有轻微
离焦/色差时理论值偏小 40% 以上，孔径半径被 `[2,7]` px 夹住（2 px vs 3 px），
匹配滤波幅度衰减到 ~0.6，`intensity` 的绝对标度随之失真。
