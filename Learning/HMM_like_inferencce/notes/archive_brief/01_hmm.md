# HMM：第一课与推导

## 1. 从概率问题出发

观测 y 是 photons，z 是未知 off/on。参数 θ={π,A,λ}：π 初始分布，A 转移概率，λ 各状态的期望 count。θ 暂时已知时求 p(z|y) 是 inference；θ 未知时估计 θ 是 learning。二者不要混同。

HMM 有两个条件独立假设：下一状态只依赖当前状态；当前观测给定当前状态后，与其他状态及观测独立。它们是模型假设，不是由“时间连续”自动成立的事实。

\[
p(z_{1:T},y_{1:T}|θ)=π_{z_1}\prod_{t=2}^T A_{z_{t-1},z_t}\prod_{t=1}^T p(y_t|z_t,θ).
\]

转移矩阵采用行表示起点，列表示终点。A_ij=p(z_t=j|z_(t−1)=i)。off 也可产生背景 photons。最小模型 y_t|z_t=k ~ Poisson(λ_k)，其中 λ_off=b，λ_on=b+a。

未知 z 的 likelihood 必须是所有路径联合概率之和。K 个状态、T 帧共有 K^T 条路径；动态规划利用共享子问题避免穷举。

## 2. 前向递推：为什么成立？

定义 α_t(k)=p(y_1:t,z_t=k)。在 z_(t−1) 上边缘化，再使用条件独立：

\[
α_1(k)=π_k B_1(k),\quad α_t(k)=B_t(k)\sum_j α_{t-1}(j)A_{jk},\quad B_t(k)=p(y_t|z_t=k).
\]

对最后状态求和得 p(y)=Σ_k α_T(k)。归一化 α_t 得 filtering：p(z_t|y_1:t)。前向递推 O(TK²)，不是从 noisy photons 先阈值判 on/off。

连续乘很小的概率会下溢。实现用 log α 与 logsumexp；另一合法办法是每帧 scaling 并记录 normalization。普通 softmax 不能替代 likelihood normalization 的 bookkeeping。

## 3. 后向递推与 smoothing

β_t(k)=p(y_(t+1):T|z_t=k)，末端 β_T(k)=1。

\[
β_t(i)=\sum_j A_{ij}B_{t+1}(j)β_{t+1}(j).
\]

由于给定 z_t 后，过去和未来条件独立：

\[
γ_t(k)=p(z_t=k|y)=α_t(k)β_t(k)/p(y).
\]

相邻状态联合后验：

\[
ξ_t(i,j)=α_t(i)A_{ij}B_{t+1}(j)β_{t+1}(j)/p(y).
\]

必须检查 Σ_k γ_t(k)=1、Σ_j ξ_t(i,j)=γ_t(i)、Σ_i ξ_t(i,j)=γ_(t+1)(j)。后一组比只检查总和更能发现 off-by-one bug。

filtering 使用当前与过去；smoothing 使用整段序列；Viterbi 求整条最可能路径 argmax_z p(z|y)。逐帧 argmax γ 并非 Viterbi，甚至可能拼出模型不允许的跳转。Viterbi 把求和换成 max，并保存 backpointer；它不提供 posterior uncertainty。

## 4. 三帧手算题（先不要运行 notebook）

π=(0.8,0.2)，A=[[0.9,0.1],[0.2,0.8]]，λ=(1,4)，y=(0,3,1)。

写出 8 条 z 路径。每条路径概率=初始概率×两个转移概率×三个 emission 概率。求和得到 evidence；加总 z_2=on 的路径再除 evidence 得 γ_2(on)。独立用 α/β 算一次。最后解释为何 y_2=3 仍不能确定 on。

## 5. Baum–Welch / EM

E step：用当前 θ 的前向后向求 γ,ξ。M step 最大化 Q(θ)=E_old[log pθ(y,z)]，不是把 argmax 状态当真值。

\[
π_k^{new}=γ_1(k),\quad A_{ij}^{new}=\frac{\sum_{t=1}^{T-1}ξ_t(i,j)}{\sum_{t=1}^{T-1}γ_t(i)},\quad λ_k^{new}=\frac{\sum_t γ_t(k)y_t}{\sum_t γ_t(k)}.
\]

多条独立序列的 counts 相加，但不跨序列计转移。这个 λ 更新只适用于自由的 state-specific Poisson rate；若 λ 由 PSF、位置、共享背景约束，则 M step 通常需要数值优化。零 occupancy、标签交换、局部最优都要检查。

## 6. 对 blinking SMLM 的适用边界

二状态 HMM 的离散停留时间为 geometric。漂白可加 absorbing state；复杂 dark states 可用多状态或 semi-Markov。曝光内发生多次切换时，“每帧一个二值状态”的 emission 可能错设，应考虑帧内 occupation。高密度时多个 emitter 的独立先验链被共同像素观测耦合，这正是 FHMM 的动机。

毕业检查：能从 joint 分布推导 α、β、γ、ξ；能说明精确算法成本；能解释为何 photon budget 降低会使 likelihood 变弱而 prior 影响增大，但 prior 不能凭空制造可辨识性。
