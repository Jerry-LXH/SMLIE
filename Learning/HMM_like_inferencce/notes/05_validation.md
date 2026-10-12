# 验证记录（2026-10-10）

- HMM 三项单元测试通过：穷举 likelihood/posterior 与 forward-backward 一致；xi 的边缘一致；2000 帧 log-domain 稳定；单帧边界正常。
- notebooks 01–03 的数值代码通过：EM likelihood 单调；FHMM posterior covariance 非零；ELBO + KL = log evidence。
- notebooks 01–03 执行时移除了绘图语句，因为可用运行时无 Matplotlib。绘图与完整 Jupyter kernel 执行尚未验证。
- 全部六个 notebook 的 Python cell 语法检查通过。
- notebooks 04–06 未做运行验证：可用 Python 无 PyTorch。所含 analytic/finite-difference/enumeration 检查需安装依赖后运行，不应把语法通过视作数值验证。
- 没有安装依赖、下载数据或训练大型模型；notebook 不带输出。
