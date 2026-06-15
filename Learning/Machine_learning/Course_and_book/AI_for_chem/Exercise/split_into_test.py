from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
import pandas as pd

def simple_split(data,test_size,target_name):
    tr_raw,ts_raw=train_test_split(data, test_size=test_size, random_state=44) # 拆分20%
    tr_target=tr_raw[target_name]
    tr_input=tr_raw.drop(target_name,axis=1)
    ts_target=ts_raw[target_name]
    ts_input=ts_raw.drop(target_name,axis=1)
    return tr_target,tr_input,ts_target,ts_input

def shuffle_split(data,test_size,target_name,shuffle_col,bins):
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