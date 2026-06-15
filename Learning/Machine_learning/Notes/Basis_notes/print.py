def print_inputs_copy(*names): # 特殊用法：用空元组实现无限封装
    print(names)
    return names

def save_all_inputs_copy(first,last,**info): # **传入实参，然后用字典封装
    profile={}
    profile["first_name"]=first
    profile["last_name"]=last
    for k,v in info.items(): # info是传入的字典
        profile[k]=v
    return profile