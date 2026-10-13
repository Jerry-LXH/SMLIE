# 自动微分入门：从链式法则到 VAE 与 Fisher

本课单独解释 notebook 04–06 的梯度计算。默认 CPU、float64、小张量。不是训练大型网络的教程。

## 1. 自动微分在做什么？

设输入 x 先平方得到 u，再计算 f=3u+2。链式法则：

```text
u = x^2
f = 3u + 2
df/dx = (df/du) × (du/dx) = 3 × 2x
x=2 时，df/dx=12
```

自动微分将运算记录成图，为每个局部操作使用已知导数，再按链式法则组合。它不是在每个参数上重新运行有限差分，也不是自动理解公式含义。

## 2. PyTorch 的最小例子

```python
import torch

torch.set_default_dtype(torch.float64)
x = torch.tensor(2.0, requires_grad=True)
f = 3 * x.square() + 2
f.backward()
print(x.grad)  # 12
```

`requires_grad=True` 表示需要对该量求导。`backward()` 将梯度累积到叶张量的 `.grad`。不同调用会累积，训练循环需要 `optimizer.zero_grad()`；通常每次优化迭代重新建立计算图。

如果只需要返回导数，而不写入 `.grad`：

```python
x = torch.tensor(2.0, requires_grad=True)
f = 3 * x.square() + 2
(df_dx,) = torch.autograd.grad(f, x)
```

## 3. gradient、Jacobian、Hessian 分别是什么？

标量输出对多个参数求导得到 gradient；向量输出对多个参数求导得到 Jacobian；标量输出的二阶导数矩阵叫 Hessian。

在 CRLB notebook 中，rate 是“每像素一个数”的向量，对位置、强度、背景求导得到 Jacobian；log likelihood 是一个标量，对这些参数求导得到 score gradient。

shape 应先画清楚：

```text
P 个像素，D 个参数
rate shape：P
rate Jacobian shape：P × D
score shape：D
score outer product shape：D × D
```

## 4. detach 切断哪条路径？

`detach()` 返回不沿原有计算图继续传梯度的张量视图。它用于明确哪些量在当前求导任务中应视为固定。

Fisher score 中观测必须固定：

```python
lam = torch.tensor(4.0, requires_grad=True)
y = torch.tensor(3.0)  # 这份观测不变
lp = torch.distributions.Poisson(lam).log_prob(y)
(score,) = torch.autograd.grad(lp, lam)
# score = 3/4 - 1 = -0.25
```

模拟时可在 `torch.no_grad()` 中生成 photons，或对结果 detach。不要把模拟器的路径梯度当作固定观测的 likelihood score。

反过来，VAE 的连续 z 通常需要保留 encoder 到 z 的梯度。如果在那里随意 detach，就会阻断 encoder 从 reconstruction 项获得梯度。

## 5. 有限差分是验证，不是同一算法

```text
导数近似 = (f(θ+h) - f(θ-h)) / (2h)
```

h 太大有截断误差，太小有浮点消减误差。float64 下试几个 h，寻找结果稳定区间；不要只用一个 h 判断。

```python
def loglik(lam, y):
    return torch.distributions.Poisson(lam).log_prob(y)

h = 1e-5
fixed_y = torch.tensor(3.0)
center = torch.tensor(4.0, requires_grad=True)
(s,) = torch.autograd.grad(loglik(center, fixed_y), center)
fd = (loglik(center.detach()+h, fixed_y)
      - loglik(center.detach()-h, fixed_y)) / (2*h)
assert torch.allclose(s, fd, atol=1e-8, rtol=1e-6)
```

## 6. 验证随机目标的梯度要固定随机源

VAE 每次重新采 ε，loss 会变化。用两份不同 ε 的 loss 做有限差分，差异可能主要来自 Monte Carlo noise。

验证 reparameterized gradient 时，固定 ε，再分别计算正负扰动的 loss。训练时仍可以每步换 ε；验证导数与训练的采样策略不是同一件事。

## 7. 对连续 z 与离散 z 的不同处理

Gaussian latent 可写 z=μ+σ×ε，梯度沿连续运算传播。二值 blinking 的采样本身不是普通连续可微操作。

小规模 blinking 模型可精确求和：前向算法由乘、加、logsumexp 等可微操作构成，参数梯度穿过这些运算，不需要穿过 sampled discrete z。

## 8. 阅读代码前必须回答

- 输出是 scalar loss、vector rate，还是 distribution 参数？
- 需要对哪个变量求导？哪些数据应固定？
- 张量 shape 是什么？参数单位是什么？
- likelihood 是否对应真实噪声模型？
- 是否有有限差分、解析解或穷举可以核对？

自动微分会准确执行错误目标的求导。先写清概率模型，再让 PyTorch 帮你算导数。
