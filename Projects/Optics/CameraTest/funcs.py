import os
import numpy as np
import tifffile as tiff
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def time_average_all_tif(folder_path):
    """
    从指定文件夹读取所有 TIFF 图像，计算逐像素平均值与标准差。
    若遇到坏文件或非 TIFF 文件会自动跳过。
    返回:
        avg_I  —— 平均图像 (float32)
        sigma_I —— 标准差图像 (float32)
    """

    # 支持 .tif 和 .tiff
    tif_files = [f for f in os.listdir(folder_path)
                 if f.lower().endswith(('.tif', '.tiff'))]

    if not tif_files:
        raise ValueError(f"❌ 目录 {folder_path} 下未找到任何 .tif / .tiff 文件")

    I_collection = []

    for fname in sorted(tif_files):
        img_path = os.path.join(folder_path, fname)
        try:
            I = tiff.imread(img_path).astype(np.float32)
            I_collection.append(I)
        except Exception as e:
            print(f"⚠️ 跳过文件 {fname}: {e}")
            continue

    if len(I_collection) == 0:
        raise RuntimeError(f"❌ 所有文件均无法读取，请检查 {folder_path}")

    # 转成 3D 数组： (N, H, W)
    I_array = np.stack(I_collection, axis=0)

    # 平均值与标准差（注意不是方差）
    avg_I = np.mean(I_array, axis=0).astype(np.float32)
    sigma_I = np.std(I_array, axis=0, ddof=1).astype(np.float32)  # 用样本标准差

    print(f"✅ 成功处理 {len(I_collection)} 张图像，共形状 {I_array.shape[1:]}")

    return avg_I, sigma_I

def linear_fit_plot(x, y, xlabel='X', ylabel='Y', title='线性拟合'):
    """
    使用numpy进行线性拟合和可视化
    """
    # 线性拟合：y = a*x + b
    x = np.array(x)
    y = np.array(y)
    a, b = np.polyfit(x, y, 1)
    
    # 计算拟合值
    y_fit = a * x + b
    
    # 计算R²
    residuals = y - y_fit
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # 创建图形
    plt.figure(figsize=(6, 4))
    
    # 绘制散点图和拟合线
    plt.scatter(x, y, alpha=0.7, label='Data Points', color='blue', s=50)
    plt.plot(x, y_fit, 'r-', linewidth=2, label=f'fit line: y = {a:.4f}x + {b:.4f}')
    
    # 添加图例和标签
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title(title, fontsize=14)
    
    # 添加拟合信息文本框
    textstr = f'k: {a:.4f}\nb: {b:.4f}\nR² = {r_squared:.4f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='top', bbox=props)
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    return a, b, r_squared


def fit_asymptote_fixed_slope(x, y, fixed_slope=0.5, method='tail_fit'):
    """
    拟合渐近线（已知斜率）
    
    Parameters:
    -----------
    x, y : 数据点
    fixed_slope : 固定的渐近线斜率
    method : 拟合方法
        - 'tail_fit' : 使用尾部数据拟合
        - 'weighted' : 加权拟合，尾部权重更高
        - 'nonlinear' : 非线性模型拟合
    """
    
    x = np.array(x)
    y = np.array(y)
    
    if method == 'tail_fit':
        # 方法1A: 使用尾部20%的数据进行拟合
        n_tail = int(0.2 * len(x))
        if n_tail < 3:  # 确保有足够的数据点
            n_tail = min(3, len(x))
        
        x_tail = x[-n_tail:]
        y_tail = y[-n_tail:]
        
        # 固定斜率拟合截距: b = mean(y) - k * mean(x)
        b = np.mean(y_tail) - fixed_slope * np.mean(x_tail)
        
        print(f"尾部拟合方法:")
        print(f"  使用最后 {n_tail} 个数据点")
        print(f"  固定斜率 k = {fixed_slope}")
        print(f"  拟合截距 b = {b:.6f}")
        
        return b, x_tail, y_tail
        
    elif method == 'weighted':
        # 方法1B: 加权拟合，给尾部数据更高权重
        weights = np.linspace(0.1, 1, len(x))  # 线性权重，尾部权重为1
        
        # 加权平均值计算截距
        weighted_mean_x = np.average(x, weights=weights)
        weighted_mean_y = np.average(y, weights=weights)
        b = weighted_mean_y - fixed_slope * weighted_mean_x
        
        print(f"加权拟合方法:")
        print(f"  权重范围: {weights[0]:.2f} - {weights[-1]:.2f}")
        print(f"  固定斜率 k = {fixed_slope}")
        print(f"  拟合截距 b = {b:.6f}")

        plt.figure(figsize=(8, 6))
        plt.scatter(x, y, alpha=0.6, color='blue', s=40, label='Points')
        x_fit = np.linspace(np.min(x), np.max(x), 100)
        y_asymptote = fixed_slope * x_fit + b
        plt.plot(x_fit, y_asymptote, 'r-', linewidth=2, 
                    label=f'AsympLine: y = {fixed_slope}x + {b:.4f}')

        plt.xlabel('Log Intensity')
        plt.ylabel('Log Noise')
        plt.title('Weighted Fit')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
        
        return b
    
    else:
        raise ValueError("方法必须是 'tail_fit' 或 'weighted'")
    
    

