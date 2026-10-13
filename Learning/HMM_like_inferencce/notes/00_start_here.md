# 讲义阅读入口

这些讲义以统计推断为中心，不要求先实现 SMLM 定位项目。它们从直觉、符号、推导、完整数值例子走到 notebook。

公式统一使用独立的纯文本块：不含 LaTeX 命令，不依赖数学渲染。正文采用普通术语和 Unicode 字母；方括号是索引，例如 y[t] 为第 t 帧。理论帧号从 1 开始，Python 数组索引从 0 开始。

## 先学 HMM，不要一次读完所有主题

打开 [HMM 第一课](01_hmm.md)，第一遍只到第 8 节。准备纸笔，用 (0,3,1) 三帧观测跟着算每一步。读完后应能解释前向量、后向量、posterior 为什么是不同东西。

第二遍读 HMM 的后半部分：相邻状态概率、Viterbi、EM、FHMM。再运行 notebook 01–02。

## 接下来读什么？

1. [VI、ELBO 与 VAE](02_learning_path.md)：第 1–7 节先解决近似分布与 ELBO，再看网络。对应 notebook 03–04。
2. [自动微分](06_autograd.md)：先能自己核对 Poisson score 和有限差分，再读 VAE 或 Fisher 的梯度代码。
3. [Fisher 与 CRLB](03_information.md)：先学单参数 Poisson，再学像素 FIM、nuisance 和未知 blinking。对应 notebook 05–06。
4. [原始阅读来源](04_reading.md)：在已理解教学例子后进入教材与论文，不把论文方法名当基础解释。

## 每课的学习标准

- 能说出已观测量、隐变量、未知参数分别是什么。
- 能解释公式中的乘法、求和、归一化分别来自哪里。
- 能在可手算的模型上核对结果。
- 能说出模型在哪种情况下会失败，而不只是跑通 notebook。

原来的三份简略讲义保存在 `archive_brief`，用于对照；第一遍以本目录的新版本为准。本次只修改讲义，不改变 notebook 的算法或训练规模。代码此前的运行验证范围仍见 `05_validation.md`。
