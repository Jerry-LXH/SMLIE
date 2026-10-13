# 第二课：从精确 HMM 到 VI、ELBO 和 VAE

先完成 HMM 第一课的三帧例子。本课重点不是记住几个方法名，而是理解：为什么要近似一个分布、怎样衡量近似、网络实际负责什么。算式均为独立纯文本块，`E_q` 表示按分布 q 取平均。

## 1. 为什么 HMM 之后需要 VI？

单 emitter HMM 可以精确边缘化。多个 emitter 的联合状态很多，而共享观测导致 posterior 耦合。即使 joint distribution 写得很简单，其 posterior 也可能很难计算。

真正的困难不是“z 未知”，因为推断本来就处理未知变量；困难是求 posterior 需要对大量状态求和，或做难算的连续积分。

**Variational inference（变分推断，VI）**选择一个容易表示的分布 q，让它尽量接近真实 posterior。原问题变成在一个分布族内优化参数。例如每个 Bernoulli 分布只需一个 on 概率；每个 Gaussian 分布需要均值和协方差。

q 是近似 posterior，不是物理 prior。生成模型规定数据如何产生，q 描述看到数据后如何反推 latent variables。

## 2. “接近”怎么定义：KL divergence

```text
KL(q || p) = sum over z: q(z) × log(q(z)/p(z))
```

用 q 的概率给 log 比值加权。连续变量把和换积分。KL 非负，两个分布相同才为零；但它不对称，不是通常的距离。

如果 q 在 p 为零的区域有正概率，KL 会无穷，因此近似分布的支持集也很重要。KL 非负可以用 Jensen 不等式证明；第一遍可以先接受这个性质，再看 textbook 的证明。

VI 常最小化 q 到 posterior 的 KL。受限 family 会造成 approximation gap；优化算法本身也可能尚未找到 family 内最优解。

## 3. ELBO 为什么出现？逐行推导

Bayes 公式给出：

```text
log posterior = log joint - log evidence
log p(z|y)    = log p(y,z) - log p(y)
```

代入 KL 的定义：

```text
KL(q || posterior)
  = E_q[log q(z) - log p(z|y)]
  = E_q[log q(z) - log p(y,z) + log p(y)]
  = E_q[log q(z)] - E_q[log p(y,z)] + log p(y)
```

最后一步里，log evidence 对 z 是常数，取平均不会改变它。整理：

```text
ELBO = E_q[log p(y,z)] - E_q[log q(z)]
log evidence = ELBO + KL(q || posterior)
```

ELBO 是 Evidence Lower Bound（证据下界）。KL 非负，所以 ELBO 不大于 log evidence。这里是 log evidence 的下界，不是定位误差或 CRLB 的下界。

固定生成参数时，log evidence 不随 q 变化。最大化 ELBO 等价于最小化这个 KL。好处是 ELBO 不直接要求算 evidence；joint 与 q 的期望仍可能需要 MC，并不是所有积分都自动消失了。

## 4. 一个两状态数值例子

假设已观测某个 y，模型给出的两项 joint probability 为 0.03 和 0.07。

```text
joint = (0.03, 0.07)
evidence = 0.10
posterior = (0.30, 0.70)
log evidence = log(0.10) ≈ -2.302585
```

先选择 q=(0.5,0.5)：

```text
ELBO = 0.5 × log(0.03/0.5) + 0.5 × log(0.07/0.5)
     ≈ -2.389762

KL = 0.5 × log(0.5/0.3) + 0.5 × log(0.5/0.7)
   ≈ 0.087177

ELBO + KL ≈ -2.302585
```

换成 q=(0.3,0.7)，两个 log 比值都等于 log(0.1)，ELBO 等于 log evidence，gap 为零。

这个小例子说明：ELBO 低于 evidence 的差额，正是 q 与 posterior 的 KL。复杂模型里 q 未必能表示精确 posterior，因此即使优化完成，gap 也可能仍大于零。

## 5. ELBO 的另一种拆法

```text
joint = prior × likelihood
ELBO = E_q[log likelihood] - KL(q || prior)
```

第一项衡量按 q 认为可能的隐状态，平均能否解释数据；第二项惩罚 q 偏离生成模型 prior。第二项中的 p 是 prior，而上一节 gap 中的 p 是 posterior，不要混淆两个 KL。

Poisson likelihood 的第一项是 Poisson log probability，不是通用平方误差。Gaussian 且方差固定时，才可把相关项写成负平方误差加常数。

## 6. mean-field 丢掉什么？

假设两个 emitter 的精确 posterior 只允许 (1,0) 或 (0,1)，各 50%：恰好一个 on，但不知道是哪一个。

每个 emitter 的 marginal on 概率均为 0.5。若将两个 marginal 相乘，四种状态各 25%，错误地产生 (0,0) 与 (1,1)。这说明 marginals 不能描述 joint dependence。

这个严格零概率例子会使该乘积与真实 posterior 的 KL 无穷；真实 Poisson 的概率通常是软的，问题仍存在但不一定无穷。

mean-field 是限制 q 的分解方式。最优 mean-field q 不一定就是精确 marginals 的乘积，因为优化的是 reverse KL，而不是逐个 marginal matching。

时序 FHMM 中：

```text
逐帧 mean-field：q = product over m,t: q_mt(z[t,m])
  丢掉时间依赖，也丢掉 emitter 间依赖。

链结构近似：q = product over m: q_m(z[1:T,m])
  每个 emitter 保留整条轨迹的依赖，但不同 emitter 间仍独立。
```

后一种是 structured variational family。每条链内部可以调用前向后向。这解释了“先学 HMM”与“后学结构化 VI”之间的联系。

## 7. 坐标更新为什么需要其他链的期望？

mean-field 可以轮流优化每一个因子：固定其他因子，更新当前因子。其通用更新结构是：

```text
log q_i(z_i) = E under other factors[log p(y,z)] + 常数
```

常数负责归一化。可以把 ELBO 中含当前因子的项与熵项组合成负 KL，最优分布即是右边指数化再归一化后的分布。

在 FHMM 中，观测 rate 是多个 emitter 的和。log likelihood 包含 log(rate)，因此“先求其他信号均值，再代入 log”通常不等于“对 log likelihood 求期望”。不能仅从观测 count 中减平均信号，就称为精确坐标 VI。

## 8. 网络为什么出现：amortized inference

普通 VI 每个新样本都单独优化 q 的参数。识别网络则把 y 映射到 q 的参数，训练后重复使用这个映射，称为 amortized inference。

约定 θ 是生成模型参数，φ 是识别网络参数。固定 θ 优化 q 是 inference；同时优化 θ 与 φ 是 model learning 加 approximate inference。

网络节约新样本的推断成本，但也可能有 amortization gap：网络输出的 q 比该 family 对这个样本单独优化的最优 q 更差。网络再大也不能突破一个过于受限的 q family。

## 9. VAE 的 encoder 和 decoder 各输出什么？

VAE 常用连续 latent、标准正态 prior，以及 Gaussian approximate posterior。

```text
prior：z ~ Normal(0,I)
encoder：y → μ(y), σ(y)
q：Normal(μ(y), diagonal covariance σ(y)^2)
decoder：z → likelihood 参数
```

encoder 输出近似 posterior 的参数，不一定是一个确定隐状态；decoder 输出观测分布的参数，例如每个像素的 Poisson rate。它们并不自动等于“真实状态”和“无噪声真图像”。

KL 到标准正态可以解析计算：

```text
KL = 0.5 × sum over latent dimensions j:
             (μ[j]^2 + σ[j]^2 - 1 - log(σ[j]^2))
```

来源是将两个 Gaussian 的 log density 相减，再利用 q 下的均值与方差求期望。

## 10. reparameterization 怎样让梯度穿过去？

```text
ε ~ Normal(0,I)
z = μ + σ × ε
```

随机源 ε 的分布不依赖网络参数。固定一次 ε 后，z 是 μ、σ 的可微函数，可以沿 encoder、z、decoder、log likelihood 反传。下一次训练再采新的 ε，平均目标近似原来的期望。

Gaussian reparameterization 不能直接用于二值 blinking。小模型可精确对离散状态求和，大模型需其他梯度方法；连续松弛会改变问题或产生近似偏差，必须说明。

## 11. 从普通 VAE 到物理约束网络

一个可解释方向是保留 Markov prior、物理 PSF 像素积分和 photon likelihood，让网络预测 variational parameters 或识别势，而不是随意替代全部生成机制。

“势”是不一定归一化的局部权重函数。网络可以给每个状态一项局部权重，再由结构化算法组合局部势和转移结构，算成一个一致的联合分布。Johnson 的 SVAE 利用这一思路与共轭图模型结构结合，不是给任意 VAE 添加时序就叫 SVAE。

DMM/DKF 使用网络表示非线性状态空间生成与识别模型。名称包含 Kalman，不意味着可做精确线性 Gaussian Kalman inference。SBI 则重点面对能模拟、但 likelihood 难求的模型，通过模拟学习 posterior、likelihood 或 ratio。

有物理 decoder 仍可能存在不可辨识性：重合 emitter、标签交换、背景与弱信号混淆。可解释性需要明确参数语义，还需要 calibration 和错设检验。

## 12. 自测与 notebook 对应

1. q 是 prior 吗？不是，它近似给定数据的 posterior。
2. ELBO 为什么是下界？log evidence 等于 ELBO 加非负 KL。
3. ELBO 上升一定说明参数更真实吗？不一定，近似 gap 和模型 evidence 都可能改变。
4. mean-field 为什么不能表示“恰好一个 on”？该约束需要 emitter 间依赖。
5. encoder 输出什么？q 的参数；decoder 输出 likelihood 的参数。
6. 为什么离散 blinking 不能直接用 Gaussian reparameterization？离散采样不是上述连续可微变换。

先做 notebook 03 的精确枚举比较，再做 notebook 04 的连续 Poisson VAE。后者用于学习梯度和 ELBO，不是 blinking SMLM 定位器。自动微分细节见 `06_autograd.md`。
