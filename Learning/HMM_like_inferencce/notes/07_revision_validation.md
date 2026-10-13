# 讲义修订检查（2026-10-12）

目标目录：`/Users/jerryling/Documents/Lab_on_chip/Learning/HMM_like_inferencce/notes`。

本次重写 HMM、VI/VAE、Fisher/CRLB 三份主讲义，新增阅读入口与自动微分讲义。原简略版本归档至 `archive_brief`。

公式使用纯文本算式块，避免 inline LaTeX 与数学渲染兼容问题。新讲义检查无 LaTeX 命令或数学分隔符，代码围栏闭合。

使用 Python 标准库核对三帧 HMM 的 forward/backward、路径穷举 evidence、所有帧的 smoothing 和 xi 边缘一致性，以及两状态 ELBO 数值例子。无需 NumPy、PyTorch 或新依赖。

本次未重新运行原 notebook，也未改变算法源文件。此前代码运行验证范围以 `05_validation.md` 为准。Markdown 结构已检查；未进行当前应用渲染的截图验收。
