# Blinking SMLM：统计推断与信息学习实验室

目标是理解统计模型、推断和信息极限，不是实现完整 SMLM 重建系统。默认中文讲解、英文术语，CPU、小数据，无真实图像下载、无 GPU 配置、无大型训练。

## 开始

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m jupyterlab
```

依赖安装由你在需要时执行；仓库本身只含文本与未执行 notebook。PyTorch 是可选学习阶段中最大的依赖，前三个 notebook 只需 NumPy/Matplotlib/Jupyter。M3 Air 默认 CPU；小模型使用 float64，避免为速度牺牲数值诊断。

按顺序打开 `notebooks/01_hmm.ipynb` 至 `06_latent_crlb.ipynb`。讲义见 `notes/01_hmm.md`、`notes/02_learning_path.md`、`notes/03_information.md`、`notes/04_reading.md`。从根目录运行 `python -m unittest discover -s tests` 验证 HMM。

## 建议节奏：8 周，每周约 5–7 小时，可按掌握程度调整

1. 概率预备：条件概率、边缘化、Bayes、期望、协方差、对数似然；推导 Poisson score。
2. HMM：手算三帧，再用前向后向验证穷举。区分 filtering/smoothing/Viterbi。
3. 参数学习：EM 的 expected counts；FHMM 两个 blinking emitter，解释后验依赖。
4. VI：推导 ELBO；理解 mean-field 丢失哪些依赖，数值比较精确后验。
5. 自动微分：固定观测求 score，有限差分核对；再理解 reparameterization。
6. VAE：小型 Poisson VAE；分清生成模型与识别网络，不把 decoder 输出当后验。
7. Fisher/CRLB：解析 Poisson FIM、MC score、nuisance 参数、奇异性。
8. 未知 blinking 的 observed FIM；阅读 SVAE/DMM/SBI，提出可解释模型设计。

每一阶段都应能：不用代码写出概率模型；说明算法实际求什么；核对一个可穷举实例；解释失败条件。不要以“跑通代码”代替掌握。

当前实现范围：HMM/FHMM/VI 为教学级实验；VAE 是连续隐变量练习，不是 SMLM 定位器；CRLB 的空间模型是 1D 像素积分 Gaussian PSF 的统计代理。真实项目的相机噪声、漂白、漂移、未知 emitter 数目和 2D/3D PSF 留给远程工作站。
