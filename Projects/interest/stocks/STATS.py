import numpy as np
import pandas as pd
from scipy.stats import skew, kurtosis, jarque_bera
from scipy.optimize import minimize
from statsmodels.tsa.stattools import adfuller
from statsmodels.stats.diagnostic import acorr_ljungbox, het_arch
from statsmodels.tsa.api import VAR
from arch import arch_model
import matplotlib.pyplot as plt
import warnings

def turn_series(series):
    '''
    Transform series into pd.series
    '''
    if isinstance(series, np.ndarray):
        series = pd.Series(series)
    elif not isinstance(series, pd.Series):
        raise TypeError("输入必须为 pandas Series 或 NumPy array")
    return series

def analyze_basic_stats(ts_array):
    """
    输入：ts_array 为一维 NumPy 数组
    输出：字典形式的统计分析结果
    """
    ts_array = turn_series(ts_array)

    # 1-4 阶矩
    mean_val = np.mean(ts_array)
    var_val = np.var(ts_array, ddof=1)  # 无偏方差
    skew_val = skew(ts_array)
    kurt_val = kurtosis(ts_array)  # 默认 Fisher，正态分布为 0

    # ADF 检验
    adf_result = adfuller(ts_array, autolag='AIC')
    adf_stat = adf_result[0]
    adf_p = adf_result[1]

    # JB 检验
    jb_stat, jb_p = jarque_bera(ts_array)

    # 打包结果
    results = {
        'Mean': mean_val,
        'Variance': var_val,
        'Skewness': skew_val,
        'Kurtosis': kurt_val,
        'ADF Statistic': adf_stat,
        'ADF p-value': adf_p,
        'JB Statistic': jb_stat,
        'JB p-value': jb_p
    }

    return results

def test_serial_dependence(series, lags=10):
    """
    对时间序列执行 Ljung–Box Q 检验 和 ARCH-LM 检验。
    参数：
        series: 时间序列，np.ndarray 或 pd.Series
        lags: 滞后阶数，默认为10
    返回：
        dict：包括 Ljung–Box 和 ARCH 检验的统计量与 p 值
    """
    series = turn_series(series)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Ljung–Box Q 检验
        lb_test = acorr_ljungbox(series, lags=[lags], return_df=True)
        lb_stat = lb_test['lb_stat'].values[0]
        lb_pvalue = lb_test['lb_pvalue'].values[0]

        # ARCH-LM 检验
        arch_test = het_arch(series, nlags=lags)
        arch_stat = arch_test[0]
        arch_pvalue = arch_test[1]

    return {
        'Ljung–Box Q Statistic': lb_stat,
        'Ljung–Box Q p-value': lb_pvalue,
        'ARCH LM Statistic': arch_stat,
        'ARCH LM p-value': arch_pvalue
    }

def bivariate_var_analysis(series1, series2, maxlags=10):
    """
    对两个时间序列进行双变量 VAR 分析，使用 AIC 选择最优滞后阶数

    参数:
    - series1, series2: Pandas Series，两个等长时间序列
    - maxlags: 最大滞后阶数（用于模型比较）

    返回:
    - var_model: 拟合后的 VAR 模型
    - lag_order: 最优滞后阶数
    """
    series1 = turn_series(series1)
    series2 = turn_series(series2)
    # 构建 DataFrame
    data = pd.concat([series1, series2], axis=1)
    data.columns = ['Series1', 'Series2']
    data = data.dropna()

    # 建模：选择滞后阶数
    model = VAR(data)
    lag_selection = model.select_order(maxlags)
    lag_order = lag_selection.aic

    print(f"选出的最优滞后阶数 (by AIC): {lag_order}")

    # 拟合 VAR 模型
    var_model = model.fit(lag_order)

    print("\n模型摘要：")
    print(var_model.summary())

    return var_model, lag_order

def vcc_mgarch_fit_exp(series1, series2, p, q, lambda_=0.94):
    """
    对两个时间序列进行VCC-MGARCH建模，返回条件协方差矩阵序列 H_t

    参数：
        series1, series2: numpy arrays 或 pd.Series, 要求已对齐
        lambda_: float, EWMA 权重参数，默认0.94

    返回：
        H_list: 条件协方差矩阵序列
        D_list: 条件标准差对角阵序列
        R_list: 条件相关矩阵序列
    """
    # 1. 统一数据格式
    series1 = np.asarray(series1) * 100
    series2 = np.asarray(series2) * 100
    T = len(series1)

    # 2. 分别拟合 GARCH(1,1)
    am1 = arch_model(series1, vol='GARCH', p=p, q=p).fit(disp="off")
    am2 = arch_model(series2, vol='GARCH', p=p, q=q).fit(disp="off")
    
    sigma1 = am1.conditional_volatility
    sigma2 = am2.conditional_volatility

    # 3. 标准化残差
    z1 = (am1.resid / sigma1)
    z2 = (am2.resid / sigma2)
    print(test_serial_dependence(z1))
    print(test_serial_dependence(z2))

    # 4. 初始化相关矩阵 R_t（EWMA）
    z_stack = np.vstack([z1, z2])
    q_t = np.cov(z_stack[:, :30])  # 用前30期初始化
    R_list, D_list, H_list = [], [], []

    for t in range(T):
        # EWMA 更新协方差矩阵
        z_t = z_stack[:, t].reshape(-1, 1)
        q_t = lambda_ * q_t + (1 - lambda_) * (z_t @ z_t.T)
        
        # 转为相关矩阵
        diag_q = np.diag(np.sqrt(np.diag(q_t)))
        inv_diag_q = np.linalg.inv(diag_q)
        r_t = inv_diag_q @ q_t @ inv_diag_q

        # 构造 D_t 对角阵
        d_t = np.diag([sigma1[t], sigma2[t]])

        # 构造 H_t = D_t R_t D_t
        h_t = d_t @ r_t @ d_t

        R_list.append(r_t)
        D_list.append(d_t)
        H_list.append(h_t)

    return np.array(H_list), np.array(D_list), np.array(R_list), am1, am2
    
def dcc_loglikelihood(params, z):
    a, b = params
    T = z.shape[1]
    Q_bar = np.cov(z)
    Q_t = Q_bar.copy()
    ll = 0.0

    for t in range(1, T):
        z_prev = z[:, t - 1].reshape(-1, 1)
        Q_t = (1 - a - b) * Q_bar + a * (z_prev @ z_prev.T) + b * Q_t

        # 标准化 Q_t 成 R_t
        D_inv = np.diag(1 / np.sqrt(np.diag(Q_t)))
        R_t = D_inv @ Q_t @ D_inv

        # 计算 log-likelihood（只与相关矩阵相关的部分）
        z_t = z[:, t].reshape(-1, 1)
        ll += np.log(np.linalg.det(R_t)) + (z_t.T @ np.linalg.inv(R_t) @ z_t)

    return ll.item()  # 返回标量


def fit_dcc_parameters(series1, series2, p=1, q=1):
    # 1. GARCH 拟合得到标准化残差
    series1 = np.asarray(series1) * 100
    series2 = np.asarray(series2) * 100

    am1 = arch_model(series1, vol='GARCH', p=p, q=q).fit(disp="off")
    am2 = arch_model(series2, vol='GARCH', p=p, q=q).fit(disp="off")

    sigma1 = am1.conditional_volatility
    sigma2 = am2.conditional_volatility

    z1 = am1.resid / sigma1
    z2 = am2.resid / sigma2
    print(test_serial_dependence(z1))
    print(test_serial_dependence(z2))
    z = np.vstack([z1, z2])

    # 2. 最小化负的 log-likelihood
    bounds = [(1e-6, 1 - 1e-6), (1e-6, 1 - 1e-6)]  # 防止 a+b > 1
    constraints = ({'type': 'ineq', 'fun': lambda x: 1 - x[0] - x[1]})

    result = minimize(dcc_loglikelihood, x0=[0.01, 0.98], args=(z,),
                      method='SLSQP', bounds=bounds, constraints=constraints)

    if result.success:
        a_opt, b_opt = result.x
        print(f"Optimal a = {a_opt:.4f}, b = {b_opt:.4f}, a + b = {a_opt + b_opt:.4f}")
    else:
        raise RuntimeError("DCC-GARCH 参数估计失败")

    return a_opt, b_opt, am1, am2

def vcc_mgarch_fit_dcc(series1, series2, p, q, a=0.01, b=0.98):
    """
    使用 GARCH(1,1) 和 DCC-GARCH 方法对两个时间序列进行建模

    参数：
        series1, series2: numpy arrays 或 pd.Series, 已对齐且平稳
        p, q: GARCH(p,q) 的阶数
        a, b: DCC-GARCH 中对相关矩阵演化的权重参数（要求 a + b < 1）

    返回：
        H_list: 条件协方差矩阵序列
        D_list: 条件标准差对角阵序列
        R_list: 条件相关矩阵序列
        am1, am2: GARCH 模型拟合对象
    """
    # 转换格式并缩放
    series1 = np.asarray(series1) * 100
    series2 = np.asarray(series2) * 100
    T = len(series1)

    # GARCH 拟合
    am1 = arch_model(series1, vol='GARCH', p=p, q=q).fit(disp="off")
    am2 = arch_model(series2, vol='GARCH', p=p, q=q).fit(disp="off")

    sigma1 = am1.conditional_volatility
    sigma2 = am2.conditional_volatility

    # 标准化残差
    z1 = am1.resid / sigma1
    z2 = am2.resid / sigma2
    z_stack = np.vstack([z1, z2])  # shape: (2, T)

    # 初始化 Q_0 为标准化残差协方差矩阵
    Q_bar = np.cov(z_stack)
    Q_t = Q_bar.copy()

    R_list, D_list, H_list = [], [], []

    for t in range(T):
        # DCC-GARCH 更新 Q_t
        if t > 0:
            z_prev = z_stack[:, t - 1].reshape(-1, 1)
            Q_t = (1 - a - b) * Q_bar + a * (z_prev @ z_prev.T) + b * Q_t

        # 标准化 Q_t 得到 R_t
        diag_q = np.diag(np.sqrt(np.diag(Q_t)))
        inv_diag_q = np.linalg.inv(diag_q)
        R_t = inv_diag_q @ Q_t @ inv_diag_q

        # 构造 D_t, H_t
        d_t = np.diag([sigma1[t], sigma2[t]])
        H_t = d_t @ R_t @ d_t

        R_list.append(R_t)
        D_list.append(d_t)
        H_list.append(H_t)

    return np.array(H_list), np.array(D_list), np.array(R_list), am1, am2