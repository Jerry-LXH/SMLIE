from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
from torch.utils import data
import pandas as pd
import random
import torch

def simple_split(data,test_size,target_name):
    """简单分割dataframe为测试训练和X，y"""
    tr_raw,ts_raw=train_test_split(data, test_size=test_size, random_state=44) # 拆分20%
    tr_target=tr_raw[target_name]
    tr_input=tr_raw.drop(target_name,axis=1)
    ts_target=ts_raw[target_name]
    ts_input=ts_raw.drop(target_name,axis=1)
    return tr_target,tr_input,ts_target,ts_input

def shuffle_split(data,test_size,target_name,shuffle_col,bins):
    """按目标列和bins分层分割dataframe为测试训练和X，y"""
    data["new_cat"] = pd.cut(data[shuffle_col], #按income分割并建立新列
                               bins=bins,
                               ) 
    split = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=44) # 分层抽样
    for train_index, test_index in split.split(data, data["new_cat"]):
        tr = data.loc[train_index]
        ts = data.loc[test_index]
    tr.drop("new_cat",axis=1,inplace=True)
    ts.drop("new_cat",axis=1,inplace=True)
    tr_target=tr[target_name]
    tr_input=tr.drop(target_name,axis=1)
    ts_target=ts[target_name]
    ts_input=ts.drop(target_name,axis=1)
    return tr_target,tr_input,ts_target,ts_input

def batches_iter(batch_size,X,y):
    """Input Tensors X, y to become generator for Batches"""
    num=len(X)
    indices=list(range(num))
    random.shuffle(indices) # 生成随机指标序列
    for i in range(0,num,batch_size):
        batch_index=torch.tensor(indices[i:min(i+batch_size,num)]) # 一个batch的所有指标
        yield X[batch_index],y[batch_index] # yelid 暂停函数，返回一个生成器

def load_array(data_arrays, batch_size, is_train=True):
    """Convert Tensor array (X, y) to a generator for Batches"""
    dataset = data.TensorDataset(*data_arrays) # 打包，数据沿第一个轴维度一致
    return data.DataLoader(dataset, batch_size, shuffle=is_train) # 返回负责随机拆分(X,y)的对象的实例

class Accumulator:
    """Accumulate on n variables"""
    def __init__(self, n):
        self.data = [0.0] * n # 列表，储存n个累计量
    
    def add(self, *args): # *args传递tuple给函数
        """Add to data new updates"""
        self.data = [a + float(b) for a, b in zip(self.data, args)] # zip组合列表中同编号元素，返回元组列表

    def reset(self):
        self.data = [0.0] * len(self.data)

    def __getitem__(self, idx): # 可用name[idx]
        return self.data[idx]