# 从 HMM 到物理约束网络

## VI 和 ELBO 的推导

从 log p(y)=log Σ_z q(z)p(y,z)/q(z) 用 Jensen 得 ELBO≤log p(y)，要求 q 支持集内 p(y,z)>0，且 q 能覆盖需要的部分。更直接的恒等式：

L(q,θ)=E_q[log pθ(y,z)]−E_q[log q(z)]；log pθ(y)−L=KL(q||pθ(z|y))≥0。

连续 latent 把和换积分。模型证据、likelihood、posterior、ELBO 各是什么必须说清。优化 q 是 inference；同时优化 θ 是 variational learning。提高 ELBO 既可能提高 evidence，也可能缩小 approximation gap，不能仅凭 ELBO 比较真参数误差。

mean-field 坐标更新：log q_i*(z_i)=E_(q_−i)[log p(y,z)]+const。链结构的 q_i 对整条 emitter 轨迹取分布，保留时间结构。HMM 的动态规划可以成为 structured VI 内的精确子程序。

## VAE

pθ(y,z)=p(z)pθ(y|z)，encoder qφ(z|y) 做 amortized inference，目标 E_q log pθ(y|z)−KL(qφ||p)。Gaussian q 使用 z=μφ+σφε, ε~N(0,I)；固定 ε 后沿确定计算图反传。

Monte Carlo score 计算中 y 必须固定；VAE 训练中则让连续 z 的 reparameterized 路径传梯度。这是两种不同求导目的。离散 z 无法直接做 Gaussian reparameterization，可精确求和、score-function estimator 或连续松弛；松弛会改变目标/引入偏差，不能默默替代 blinking 状态。

## 三条方法学路线的区别

SVAE：网络输出 conjugate potentials，把识别网络与 graphical-model inference 组合；不要把任意时序 VAE 都叫 Johnson 的 SVAE。

DMM/DKF：非线性 transition/emission 与 structured recognition 网络组成状态空间模型。这里的 Kalman 名称不保证模型仍能做精确线性 Gaussian Kalman inference。

SBI：当 likelihood 难求但 simulator 可用时，学习 posterior、likelihood 或 likelihood ratio。已有可微且可精确边缘化的小 HMM 时，先用精确 likelihood 作基准；SBI 不等于自动得到 Fisher 或校准 posterior。

## 可解释性设计任务（暂不实现全项目）

画出 z_tm → λ_tp → y_tp；位置 r_m、brightness a_m、背景 b_p → λ_tp；A → z_tm。保留物理 emission 和明确的状态语义，让网络学习 q 的参数/势函数。网络若自由修改 PSF 或 transition，需标明增加了哪些可辨识性问题。

每个候选模型必须回答：哪些参数具有物理单位？哪些是 nuisance？哪些由数据可辨识？网络近似的是何种分布？后验区间是否覆盖真值？低 photon、重叠、背景错设时会怎样？先用两个 emitter 的可穷举模型检验，再考虑远程机器。
