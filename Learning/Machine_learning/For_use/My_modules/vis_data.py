import matplotlib as mpl
import matplotlib.pyplot as plt
"""该模块包含常用作图相关函数"""

# 对数据中数值项按50个等同区间作分布直方图
def show_hist(data):
    mpl.rc('axes', labelsize=14) # axis font size
    mpl.rc('xtick', labelsize=12)# x,y-label font size
    mpl.rc('ytick', labelsize=12) # 全局图片生成样式设置
    data.hist(bins=50, figsize=(20,15))
    plt.show()

