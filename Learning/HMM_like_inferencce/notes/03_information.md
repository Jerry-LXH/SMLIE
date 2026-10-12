# Fisher 信息和 CRLB：计算流程

## 1. 明确实验与参数

固定 θ，重复实验 Y~pθ。score s(Y)=∇θ log pθ(Y)，F(θ)=Eθ[ssᵀ]。正则条件成立时 E[s]=0，F=−E[Hessian log pθ(Y)]。单次 observed Hessian 不是 expected Fisher；MC score outer product 是 PSD，单次负 Hessian 未必 PSD。

无偏 estimator 满足 Cov(θhat) ⪰ F^−1。位置单位为 nm 时 F_xx 单位为 nm^−2，sqrt([F^−1]_xx) 才是 nm 的 SD bound。固定参数、随机 prior 的 Bayesian bound 和 biased estimator 的界不能混用。

## 2. 独立 Poisson 情形

log p(y|θ)=Σ_p[y_p log λ_p−λ_p−log(y_p!)]。
score_i=Σ_p(y_p/λ_p−1)∂_i λ_p。
利用 E[y_p]=λ_p、Var(y_p)=λ_p、不同像素独立，得 F_ij=Σ_p(∂_i λ_p)(∂_j λ_p)/λ_p。这提供解析参考以校验自动微分和 MC，而不是依赖一次数值逆。

## 3. MC + autodiff 实操

- 选 θ0、单位、已知参数和 nuisance 参数；写 normalized likelihood。
- 在 θ0 模拟独立实验，保存且 detach 样本。一个完整 correlated sequence 算一个样本。
- 对同一观测求 ∇θ log pθ(y)|θ0。不要梯度穿过采样器。
- 平均 score outer products。检查 score 均值，用独立 blocks 估计 MC SE；提高 M 看收敛。
- 有解析解时核对；无解析解时做有限差分和 tiny exact enumeration。
- 检查 F 对称性、特征值、单位缩放、病态方向，再用 solve 解线性系统。不要把 pseudoinverse 当作不可辨识参数的有限 CRLB。

## 4. nuisance parameters

θ=(ψ,η)。若 η 未知，有效信息 F_eff=F_ψψ−F_ψη F_ηη^−1 F_ηψ；若 η 已知，只有 F_ψψ。矩阵分块可逆且正则时 [F^−1]_ψψ=F_eff^−1。位置、brightness、背景通常相关，因此只取 F_xx 倒数不够。

若 positivity 用 log N 参数化，F 是新坐标下的信息。还原物理参数 covariance 必须用 Jacobian：Cov(θhat) bound =J F_u^−1 Jᵀ，不能直接把 log 参数 SD 当 photon 或 nm。

## 5. hidden blinking 的信息

pθ(y)=Σ_z pθ(y,z)。观测 score=Eθ[complete score|y]（在可交换导数与和的条件下）。全方差分解给 missing-information identity：F_complete−F_observed=E[Cov(complete score|Y)] ⪰0。

因此把模拟 z 当真值会算 complete-data bound；物理实验没有观测 z 时应对 forward likelihood 求导。这不是某一条模拟轨迹 conditional FIM 与 observed FIM 的逐样本排序，而是期望信息的关系。

以 ELBO 梯度替代 marginal likelihood score，一般不再是原模型 Fisher。对多 emitter 可精确 marginalize 的规模很小；扩展时必须量化近似误差。

## 6. CRLB 的使用边界

低 photon 时局部无偏条件可能不适合实际 estimator；posterior variance 也不等于 CRLB。多个等价 emitter 标签、位置重合、未知数量、参数边界都会影响正则性。Fisher 是局部量，不能揭示全部多峰 ambiguity。增加 blinking dynamics 并不自动超过独立模型界，应先统一实验设计和 nuisance 假设。
