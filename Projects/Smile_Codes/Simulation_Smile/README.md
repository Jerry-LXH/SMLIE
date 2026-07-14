# simulation_smile

可组合的显微镜成像模拟框架。覆盖 **镜头/PSF、照明光斑、荧光分子光动力学（闪烁/漂白）、相机、噪声** 全链路，
设计为在 Jupyter notebook 中逐环节调试、按需拼装。

> 当前状态：**架构骨架**。核心接口（`core/`）已实现可用；各成像环节为带完整签名与文档的占位（`NotImplementedError`），
> 后续在 notebook 中逐个填充实现。

## 设计：配置 + 组件 + 管线

三层解耦，让分类好的函数/类天然可组合：

| 层 | 角色 | 例子 |
|----|------|------|
| **配置层** `@dataclass` | 纯数据、可序列化的参数 | `PSFConfig`, `CameraConfig`, `NoiseConfig` |
| **组件层** `Stage` | 由配置构造的实现类，统一 `apply(ctx)` 接口 | `PupilPSF`, `BlinkingModel`, `SCMOSCamera`, `Poisson` |
| **管线层** `Pipeline` | 按顺序串联组件，渲染单帧/影像序列 | `SimulationPipeline` |

### 可组合的关键：统一的 `Stage` 协议 + 共享 `SimContext`

每个组件都实现一个方法：

```python
def apply(self, ctx: SimContext) -> SimContext: ...
```

所有组件读写同一个 `SimContext`（含 `grid / emitters / field / image / rng / meta`）。
因为大家共享同一个对象、各管各的字段，所以 stage 可以**随意增删、替换、重排**。
中间产物写入 `ctx.meta`（如 `ctx.meta['psf']`），方便在 notebook 里取出可视化调试。

固定的物理执行顺序：

```
photophysics → illumination → optics → camera → noise
```

任一环节都可省略（如去掉 photophysics 即得静态单帧）。

## 目录

```
simulation_smile/
  core/         # 地基：Stage 协议、SimContext、Grid、Emitters、backend(xp)、registry、coords、units
  optics/       # 镜头 + PSF：pupil（光瞳）、zernike（像差）、psf（PupilPSF/GaussianPSF）、lens
  illumination/ # 照明光斑：UniformField / GaussianSpot / TIRF
  fluorophore/  # 光动力学：photophysics（闪烁/漂白状态机）、emitters（几何）、spectra
  camera/       # 相机：pixelate（下采样）、sensors（IdealCamera/EMCCD/SCMOSCamera）
  noise/        # 噪声：shot（Poisson）、readout（读出）、background（背景/暗电流）
  pipeline/     # scene（场景）、pipeline（builder + render）、timeline（时间序列引擎）
  io/           # 多帧 TIFF、配置序列化
  viz/          # PSF/光瞳/帧/movie/光强曲线绘图
  configs.py    # 聚合所有配置：import configs as C
  presets.py    # 预设场景：smlm_storm() / widefield() / confocal()
notebooks/      # 分环节调试 notebook（00 快速上手 ~ 05 整链）
tests/          # PSF / 光动力学 / 端到端 测试
```

## 安装

```bash
cd Projects/Smile_Codes/Simulation_Smile
pip install -e .          # 核心仅依赖 numpy/scipy/matplotlib/tifffile
pip install -e ".[dev]"   # 含 pytest / pyyaml
# 可选 GPU：pip install -e ".[gpu]"  然后 simulation_smile.set_backend("cupy")
```

## 三种组装方式

```python
from simulation_smile import presets, SimulationPipeline, configs as C
from simulation_smile.optics import PupilPSF, GaussianPSF
from simulation_smile.illumination import GaussianSpot
from simulation_smile.fluorophore import BlinkingModel
from simulation_smile.camera import SCMOSCamera
from simulation_smile.noise import Poisson, ReadNoise
import numpy as np

# 1) 预设一键起步
movie = presets.smlm_storm().render_movie(n_frames=200)

# 2) 手动拼装，替换任意环节
pipe = (SimulationPipeline(grid=C.Grid(pixel_size=108e-9, size=64))
        .with_optics(PupilPSF(C.PSFConfig(na=1.49, wavelength=668e-9, zernike={5: 0.3})))
        .with_illumination(GaussianSpot(C.IlluminationConfig(waist=20e-6)))
        .with_photophysics(BlinkingModel(C.FluorophoreConfig(k_on=0.5, k_off=20, k_bleach=0.01)))
        .with_camera(SCMOSCamera(C.CameraConfig(qe=0.95, gain=2.0)))
        .with_noise([Poisson(), ReadNoise(C.NoiseConfig(read_sigma=1.5))]))

# 3) 单独调试某一环节（不进管线）
psf = PupilPSF(C.PSFConfig(na=1.49)).psf_stack(z=np.linspace(-1e-6, 1e-6, 21))
```

## GPU 后端

所有数值代码通过 `core.backend.xp` 取数组模块。一处切换全局 CPU/GPU：

```python
import simulation_smile as sm
sm.set_backend("cupy")   # 之后所有 stage 在 GPU 上运行
```

## 参考实现

- 光瞳→PSF / 噪声 / 非均匀背景：`ML_for_Smile/Paper_Codes/DeepSTORM3D-master/DeepSTORM3D/physics_utils.py`
- PSF / Zernike / 随机点 / 下采样：`Codes_revised_JL/+gen/`（`gen_PSF3D.m` 等）
- 照明高斯拟合：`Codes_revised_JL/+optics/+illumi/`
