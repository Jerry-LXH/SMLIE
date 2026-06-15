import os
import matplotlib.pyplot as plt
import tarfile
import urllib.request
import joblib
"""该模块含有和操作系统相关的函数"""

def create_dir(PROJECT_ROOT_DIR,data_class,second_class=None): # data_class 为图片/数据/输出/模型等
    """建立根目录下名为dataclass的路径并返回之"""
    if second_class==None:
        aim_path = os.path.join(PROJECT_ROOT_DIR,data_class)
    else:
        aim_path = os.path.join(PROJECT_ROOT_DIR,data_class,second_class)
    if not os.path.isdir(aim_path):
        os.makedirs(aim_path)
        print("Dir ["+aim_path+"] has been made!\n")
    else:
        print("Dir ["+aim_path+"] already exists!\n")
    return aim_path 

def save_fig(fig_name,second_class=None, tight_layout=True, fig_extension="png", resolution=300,PROJECT_ROOT_DIR = "."):
    """保存matplot对象中的图片到image路径"""
    if second_class==None:
        IMAGES_PATH = create_dir(PROJECT_ROOT_DIR,"image")
    else:
        IMAGES_PATH = create_dir(PROJECT_ROOT_DIR,"image",second_class)
    full_path = os.path.join(IMAGES_PATH, fig_name + "." + fig_extension) # 图片文件完整路径
    print("Saving figure ["+fig_name+"] in ["+IMAGES_PATH+"]...\n")
    if tight_layout:
        plt.tight_layout()
    plt.savefig(full_path, format=fig_extension, dpi=resolution) # pyplot的保存功能，需要先生成图像
    print("Yah! Figure ["+full_path+"] has been saved!\n")
    return full_path

def fetch_data_from_url(url,data_name,extension="tgz",second_class=None,PROJECT_ROOT_DIR = "."):
    """下载url中的tgz数据到data路径"""
    DATA_PATH = create_dir(PROJECT_ROOT_DIR,"data",second_class)
    full_path = os.path.join(DATA_PATH, data_name + "." + extension) # 图片文件完整路径
    print("Saving data ["+data_name+"] in ["+DATA_PATH+"]...\n")
    urllib.request.urlretrieve(url, full_path) # 下载
    print("Yah! ["+url+"] has been downloaded to ["+full_path+"]!\n")
    return full_path

def de_compress_tgz(tgz_full_path,aim_path=None):
    """解压缩tgz文件，默认到压缩文件所在目录"""
    housing_tgz = tarfile.open(tgz_full_path) # 打开压缩文件
    tgz_path=os.path.dirname(tgz_full_path)
    if aim_path==None:
        housing_tgz.extractall(path=tgz_path) # 解压到当前目录
        print("["+tgz_path+"] has been extracted to its dir!\n")
    else:
        housing_tgz.extractall(path=aim_path)
        print("["+tgz_path+"] has been extracted to ["+aim_path+"]!\n")

def save_sklearn_model(model,model_name,second_class=None,PROJECT_ROOT_DIR = "."):
    """保存模型到model路径"""
    if second_class==None:
        MODEL_PATH = create_dir(PROJECT_ROOT_DIR,"model")
    else:
        MODEL_PATH = create_dir(PROJECT_ROOT_DIR,"model",second_class)

    full_path = os.path.join(MODEL_PATH, model_name + "." + "pkl") # 图片文件完整路径
    joblib.dump(model, full_path)
    print("Yah! ["+model_name+"] has been saved to ["+full_path+"]!\n")
    return full_path