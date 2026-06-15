### 提供作图需要的一些函数

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Heiti TC'] # 用来正常显示中文标签
import os # os模块是与操作系统交互的模块，用来操作文件
from sklearn.metrics import confusion_matrix
# 本项目目录config.py文件中的ORDER变量(序号数组，confusion_matrix作图时使用)与STRAINS变量(菌株名称的字典)
from config import ORDER, STRAINS
import seaborn as sns

def plot_spectra(x, ys, title, labels=None, dy=0):  # labels, dy两个参数有缺省值
  """
  x: x轴坐标
  ys: y轴坐标，支持一条光谱数据(1D数组)与多条数据(2D数组)
  title: 图的标题
  labels：每条数据的标签
  dy: 多条数光谱数据偏移高度，主要目的是不使数据的谱线重合
  """
  if isinstance(ys, np.ndarray) and len(ys.shape) == 1:
      ys = np.array([ys])  # 使ys为2D数组
      
  plt.figure(figsize=(10,6))
  if x is None:  # 如果没传x轴坐标，则使用最后一个维度上的数组序号作为x轴坐标
      x = np.array(list(range(ys.shape[-1])))      
      
  for i, y in enumerate(ys):
      if labels is not None and len(labels) > 0:  # 如果有传labels则作图时使用label
          plt.plot(x, y + dy * i, label=labels[i])
      else:
          plt.plot(x, y + dy * i)
  plt.xlabel('位移(cm$^{-1}$)')
  plt.ylabel('强度(a.u.)')
  plt.title(title)
  if labels is not None and len(labels) > 0:  # 如果有传labels则作图时显示出label
    plt.legend()
  plt.show()

def plot_cm(cm, labels):
    """
    根据名字绘制混淆矩阵图。

    参数:
    cm : 2D array-like, 混淆矩阵
    labels : list, 混淆矩阵中对应的类别标签
    """
    plt.figure(figsize=(8, 6)) # 创建一个颜色图，并设置图形尺寸
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=labels, yticklabels=labels, cbar=False)# 使用seaborn的heatmap来绘制混淆矩阵，cmap选择颜色风格，annot显示数字
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Confusion Matrix')
    plt.show()
    

def plot_confusion_matrix(y_test, y_pred): 
    """
    计算confution matrix，并根据SRRAINS中的名字作为名字绘制热图
    """
    label = [STRAINS[i] for i in ORDER] # 字符串列表
    cm = confusion_matrix(y_test, y_pred, labels=ORDER)      
    # Plot confusion matrix
    plot_cm(cm, label)