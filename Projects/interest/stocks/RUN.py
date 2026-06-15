import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from STATS import *

def run_analysis(name):
    '''
    run the 2-list-analysis
    '''
    if name == '总体':
        data_us = pd.read_excel('data/美股'+name+'日数据.xls').dropna()
        data_china = pd.read_excel('data/A股'+name+'日数据.xls').dropna()
    else:
        data_us = pd.read_excel('data/美股'+name+'日数据.xlsx').dropna()
        data_china = pd.read_excel('data/A股'+name+'日数据.xlsx').dropna()
    date_name_us = '日期'
    date_name_china = '日期'
    value_name_us = '收盘价(总股本加权平均)\n[单位] 元'
    value_name_china = '收盘价(总股本加权平均)\n[单位] 元'
    data_us[date_name_us] = pd.to_datetime(data_us[date_name_us])
    data_china[date_name_china] = pd.to_datetime(data_china[date_name_china])
    common_values = set(data_us[date_name_us]) & set(data_china[date_name_china])
    data_us = data_us[data_us[date_name_us].isin(common_values)]
    data_china = data_china[data_china[date_name_china].isin(common_values)]
    data_us.set_index(date_name_us, inplace=True) 
    data_china.set_index(date_name_china, inplace=True)
    list_us = np.array(data_us[value_name_us].to_list())
    list_china =np.array( data_china[value_name_china].to_list())
    r_us = np.log(list_us[1:] / list_us[:-1])
    r_china = np.log(list_china[1:] / list_china[:-1])
    
    ### BASIC STAT Part
    index = ['US_'+ name, 'China_'+ name]
    basic_us = analyze_basic_stats(r_us)
    basic_china = analyze_basic_stats(r_china)
    self_test_us = test_serial_dependence(r_us)
    self_test_china = test_serial_dependence(r_china)
    self_test_us = {f'{k}_1': v for k, v in self_test_us.items()}
    self_test_china = {f'{k}_1': v for k, v in self_test_china.items()}
    stats_basic = pd.DataFrame([basic_us, basic_china], index=index)
    stats_self_1 = pd.DataFrame([self_test_us, self_test_china], index=index)

    ### VAR Part
    maxlags = 15
    var_model, lag_order = bivariate_var_analysis(r_us,r_china, maxlags=maxlags)
    coefs = var_model.coefs 
    p = coefs.shape[0]
    k = coefs.shape[1]
    rows = []
    for eq in range(k):  # eq 是方程的索引
        row = []
        # 遍历每个滞后阶
        for lag in range(maxlags):
            if lag >= p:  # 如果当前滞后阶大于最优滞后阶数，将系数设为 0
                coeffs = [0] * k  # 设置为零
            else:
                coeffs = coefs[lag, eq]  # 获取当前滞后阶的系数，形状是 (k,)
            
            row.extend(coeffs)  # 拼接当前滞后阶的系数
        rows.append(row)  # 将行加入列表

    print(rows)
    col_names = [f"L{lag+1}_{var}" for lag in range(maxlags) for var in ['1','2']]
    stats_var = pd.DataFrame(rows, index=index, columns=col_names)


    residuals = var_model.resid  # DataFrame: shape (n_obs, n_vars)
    resid_us = np.array(residuals.iloc[:, 0].tolist())  # 第一个变量的残差序列
    resid_china = np.array(residuals.iloc[:, 1].tolist())  # 第二个变量的残差序列
    self_test_us = test_serial_dependence(resid_us)
    self_test_china = test_serial_dependence(resid_china)
    self_test_us = {f'{k}_2': v for k, v in self_test_us.items()}
    self_test_china = {f'{k}_2': v for k, v in self_test_china.items()}
    stats_self_2 = pd.DataFrame([self_test_us, self_test_china], index=index)

    ### DCC GARCH Part
    a_fit, b_fit, am1, am2 = fit_dcc_parameters(resid_us, resid_china)
    H, D, R, _, _ = vcc_mgarch_fit_dcc(resid_us, resid_china, p=1, q=1, a=a_fit, b=b_fit)

    basic_us = am1.params.to_dict()
    basic_china = am2.params.to_dict()
    fit_self = pd.DataFrame([basic_us, basic_china], index=index)

    basic_us = {'a_dcc':a_fit, 'b_dcc':b_fit}
    basic_china = {'a_dcc':a_fit, 'b_dcc':b_fit}
    fit_dcc = pd.DataFrame([basic_us, basic_china], index=index)

    resid_us = am1.resid  # 第一个变量的残差序列
    resid_china = am2.resid  # 第二个变量的残差序列
    sigma1 = am1.conditional_volatility
    sigma2 = am2.conditional_volatility
    z1 = am1.resid / sigma1
    z2 = am2.resid / sigma2
    self_test_us = test_serial_dependence(z1)
    self_test_china = test_serial_dependence(z2)
    self_test_us = {f'{k}_3': v for k, v in self_test_us.items()}
    self_test_china = {f'{k}_3': v for k, v in self_test_china.items()}
    stats_self_3 = pd.DataFrame([self_test_us, self_test_china], index=index)



    stats = pd.concat([stats_basic, stats_self_1, stats_var, stats_self_2, fit_self, fit_dcc, stats_self_3], axis=1)
    return stats