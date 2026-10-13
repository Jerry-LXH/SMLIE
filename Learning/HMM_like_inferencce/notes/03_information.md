# 第三课：Fisher 信息与 CRLB，从一个 Poisson 数开始

本课的关键是区分“对什么求导”“对什么取平均”“界约束什么”。先做可以手算的标量例子，再扩展到位置、多参数与 hidden blinking。

## 1. Fisher 在衡量什么？

参数决定数据分布。如果微小改变参数会明显改变观测分布，相近参数较容易区分；如果观测分布几乎不变，就较难区分。Fisher information 衡量这种局部敏感性，是实验与模型的性质。

它不是网络 loss 越低就越高的指标，也不是一个测量已经包含“真位置”的证明。

定义 score 时，数据固定，参数变化：

```text
score(y; θ) = derivative with respect to θ of log p(y; θ)
```

Fisher 是在固定真参数下，反复采集数据时 score 平方的平均：

```text
F(θ) = E over Y drawn from p(·;θ): score(Y;θ)^2
```

先对 likelihood 求导，再对重复实验的数据分布平均，是两步不同操作。

## 2. 手算 Poisson：估计平均 photons

设一次观测为 Poisson(λ)，目标参数就是 λ。

```text
log p(y;λ) = y log λ - λ - log(y!)
score(y;λ) = y/λ - 1
```

y 固定，因此阶乘项对 λ 的导数是零。若 y 大于 λ，score 为正，这份数据倾向提高 λ；若小于 λ，倾向降低 λ。

```text
E[score] = E[Y]/λ - 1 = 0

F(λ) = E[(Y/λ - 1)^2]
     = E[(Y-λ)^2] / λ^2
     = Var(Y) / λ^2
     = 1/λ
```

这里用到了 Poisson 均值与方差均为 λ。独立重复 n 次，Fisher 变成 n/λ。

不要由 1/λ 下降推断“光子越多定位越差”。这里参数是 photon 均值，描述其绝对误差；相对误差的界仍会下降。估计位置时参数和导数不同。

## 3. CRLB 约束的是哪种误差？

经典标量 CRLB 在适当正则条件及无偏条件下给出：

```text
Var(θ_hat) >= 1/F(θ)
```

无偏指在同一个真参数下重复实验，估计结果的均值等于真参数。对 Poisson 的 n 次观测：

```text
Var(λ_hat) >= λ/n
SD(λ_hat)  >= sqrt(λ/n)
```

样本均值正好达到这个界。举例 λ=4、n=100，方差界 0.04，标准差界 0.2，相对标准差 5%。

一个简短推导说明无偏条件的作用：

```text
E[θ_hat(Y)] = θ
对 θ 求导，允许交换导数与积分后：E[θ_hat × score] = 1
因为 E[score]=0：E[(θ_hat-θ) × score] = 1
Cauchy–Schwarz：1 <= Var(θ_hat) × E[score^2]
整理：Var(θ_hat) >= 1/F
```

若估计有偏，左边一般不再是 1。边界、参数相关支持集、不可辨识性也会影响正则条件。

CRLB 不是每次实验误差的上限，不是 posterior 区间，也不是所有神经网络 RMSE 的通用下界。有偏估计的 MSE 等于 variance 加 bias²；不能只比较其 variance 与无偏 CRLB。

## 4. 多参数 Fisher 矩阵

真实定位常有位置 x、signal N、背景 b。score 是三维梯度，Fisher 是 outer product 的平均：

```text
θ = (x,N,b)
s[i] = derivative of log likelihood with respect to θ[i]
F[i,j] = E[s[i] × s[j]]
Cov(θ_hat) >= inverse(F)  （矩阵半正定意义）
```

非对角项表示不同参数的 score 相关。背景增加与信号增加可能部分解释同一种图像变化，因此它们会与位置估计互相影响。

位置精度界要取：

```text
位置 variance bound = inverse(F)[x,x]
位置 SD bound       = sqrt(inverse(F)[x,x])
```

这一般不同于 1/F[x,x]。后者对应其余参数已知的标量问题。

## 5. 从像素 Poisson likelihood 推出解析 FIM

给定参数，各像素独立，rate 为 λ[p]：

```text
log likelihood = sum over p:
                   (y[p] log λ[p] - λ[p] - log(y[p]!))

s[i] = sum over p: (y[p]/λ[p] - 1) × derivative of λ[p] wrt θ[i]
```

计算两个 score 的乘积期望时，不同像素交叉项因独立且均值为零消失。同一像素项的平均平方是 1/λ[p]。因此：

```text
F[i,j] = sum over p:
           derivative(λ[p],θ[i]) × derivative(λ[p],θ[j]) / λ[p]

J[p,i] = derivative(λ[p],θ[i])
F = transpose(J) × diag(1/λ) × J
```

J 是 rate 对参数的 Jacobian，不是任意训练 loss 的 Hessian。

位置模型举例：

```text
λ[p] = b + N × h[p](x)
derivative(λ[p],x) = N × derivative(h[p](x),x)
```

h 是 PSF 像素积分。你的 optics 知识在 h 上，统计推导从 rate 开始。这里要求 λ 为正；零均值的参数边界不应直接代入一般公式。

## 6. Monte Carlo + autodiff 的完整流程

选定真参数 θ0，模拟 M 份独立数据。每份数据固定不变，求该观测的 likelihood gradient：

```text
y[m] ~ p(·;θ0)
s[m] = gradient wrt θ of log p(y[m];θ), evaluated at θ0
F_MC = sum_m outer(s[m],s[m]) / M
```

对应流程：

```python
samples = simulate(theta0).detach()
scores = [grad_log_likelihood(theta0, y) for y in samples]
Fmc = sum(outer(s, s) for s in scores) / len(scores)
```

模拟只是生成固定观测，不对采样器求导。VAE 训练希望连续 latent 的路径传梯度；Fisher score 要对固定数据求导，这两种任务不能混用。

先核对单参数 Poisson 的解析 1/λ，再核对像素解析 FIM。没有解析结果时，有限差分和 tiny exact enumeration 很重要。

## 7. MC 结果怎样知道是否可信？

检查 score 均值是否接近零，F 是否对称、是否半正定。score outer product 的平均天然半正定，但接近奇异时仍难稳定求逆。

把样本分成多个独立 blocks，每 block 单独估计 F，用各 block estimates 的样本标准差除以 block 数的平方根，得到矩阵各元素的 MC 标准误差。增加样本看结果是否在不确定性范围内收敛。

标准误差应基于独立实验单位。如果原实验是一整段相关 blinking 序列，一份 MC 样本应是整段序列；不能将其中每一帧当成独立重复。

正则条件下 Fisher 也等于负 log likelihood Hessian 的期望。但一份数据的负 Hessian 不等于 expected Fisher，也未必半正定，不能随便互换。

## 8. nuisance 参数为何会降低有效信息？

nuisance parameter 不是“不重要的参数”，而是目标不是它、但它未知。例如我们关心位置，却不能假装背景已知。

```text
ψ = 目标参数，η = nuisance 参数
F_eff = F_ψψ - F_ψη × inverse(F_ηη) × F_ηψ
```

在分块可逆时，目标参数 covariance 的界等于 inverse(F_eff)。第二项扣掉 nuisance 可解释的变化。

纯矩阵例子：

```text
F = [[4,2],
     [2,2]]

第二参数已知：第一参数 variance bound = 1/4 = 0.25
第二参数未知：F_eff = 4 - 2 × (1/2) × 2 = 2
              第一参数 variance bound = 1/2 = 0.50
```

这个例子只说明矩阵关系，不代表具体 SMLM 实验。

## 9. 单位、参数化与奇异性

x 使用 nm 时，F[x,x] 单位为 nm 的负二次方，最终位置 SD 单位才是 nm。不同参数尺度会影响 condition number，应检查有意义的尺度化参数，而不是仅报告一个原始条件数。

若用 u=log N 保证强度为正，自动微分算到的是 u 的 Fisher。还原物理参数 covariance 的界需乘参数变换的 Jacobian：

```text
J = derivative of physical θ with respect to transformed u
物理参数 covariance bound = J × inverse(F_u) × transpose(J)
```

F 奇异可能说明参数不可辨识。不要自动加 ridge 或用 pseudoinverse 宣称不可辨识参数有有限 CRLB。先解释哪些参数方向不能由数据区分。

## 10. 未知 blinking：必须对 marginal likelihood 求导

模拟时知道 z，真实观测没有 z。这对应两个不同实验：

```text
complete-data likelihood：p(y,z;θ)
observed-data likelihood：p(y;θ) = sum over z: p(y,z;θ)
```

若直接用模拟 z 算 score，等于额外观测到 blinking 真状态，通常得到过于乐观的信息。真实实验应对后者求导。

小 HMM 可以让 autodiff 穿过前向算法：先精确计算 log marginal likelihood，再对参数求 gradient。不能先挑一条最可能路径再当真值求导。

适当条件下：

```text
observed score(y) = E[complete score(y,Z) | y]
F_complete - F_observed
  = E[Cov(complete score | Y)]
```

右边半正定，表示 hidden z 带来的 missing information。未知状态将一些参数敏感性平均掉了。

这比较的是同一参数化与实验下的期望信息矩阵，不保证两次有限 MC 估计严格排序。参数依赖的转移与 prior 也必须包含在 likelihood 中；若当它们已知，得到的是不同的 bound。

## 11. 不要混淆三个量

ELBO 是 log evidence 的下界；Fisher 是 likelihood score outer product 的期望；posterior variance 是给定某次观测后，参数 posterior 的方差。

用 ELBO gradient 替代 marginal likelihood score，一般不是原模型 Fisher。posterior interval 不自动等于 CRLB，只有某些正则大样本情形下才可能近似关联。

Fisher 还是局部量：标签交换、多峰歧义、完全重合、未知 emitter 数量都可能超出普通固定维度 CRLB 的解释范围。

## 12. 学习顺序和自测

先做 notebook 05：解析像素 FIM、MC score、finite difference。再做 notebook 06：marginal likelihood 与路径穷举比较，区分 observed 和 complete 信息。

1. 求 score 时固定什么？固定观测。
2. Fisher 的期望对什么取？固定真参数下重复实验的数据。
3. 位置 SD 是哪个矩阵元素？逆 Fisher 的位置对角项再开平方。
4. 模拟 z 能直接用于未知 blinking 的界吗？不能，那是 complete-data 实验。
5. F 奇异怎么办？先研究可辨识性，不自动通过正则化创造信息。
6. 为何先算简单 Poisson？它能解析核对，能发现错误。
