import matplotlib as mpl
import matplotlib.pyplot as plt
def show_hist(data):
    mpl.rc('axes', labelsize=14) # axis font size
    mpl.rc('xtick', labelsize=12)# x,y-label font size
    mpl.rc('ytick', labelsize=12) # 全局图片生成样式设置
    data.hist(bins=50, figsize=(20,15))
    plt.show()