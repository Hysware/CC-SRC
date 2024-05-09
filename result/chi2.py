###########################################################################################
# -*- coding: utf-8 -*-
################# 运行环境 python3
###这个程序可以打开原始数据文件与计算数据，并计算两者的chi^2
######C：cal  E:exp
# 从shell 中读取5个参数，分别为C和E文件的存储位置以及要读取的列数,以及输出的拟合图片的位置路径
import sys
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
C_file_pth = sys.argv[1]  # cal文件路径
C_column  =sys.argv[2]    # 读取的文件列数cal
E_file_pth=sys.argv[3]
E_column  =sys.argv[4]
fig_info=sys.argv[5]

#### 这个函数用于读取理论计算结果和实验结果，并将其返回为4个list类型的数据,分别是实验和理论的有效能量以及截面比

def get_data(C_file_pth,C_column,E_file_pth,E_column):
    C_file_path = C_file_pth  # 打开并读取计算文件
    try:
        with open(C_file_path, 'r') as file:
        # 读取文件的每一行
            lines = file.readlines()

        # 指定要提取的列数（从0开始计数）
            target_column = int(C_column)  # 这里假设你想提取第二列数据
            C_array=[]
            Cx_array=[]                    
        # 提取指定列的数据并打印出来
            for line in lines:
            # 分割每一行的数据（假设数据之间使用空格或制表符分隔）
               columns = line.strip().split()  # 使用split('\t')来处理制表符分隔的文件
               if len(columns) > target_column:
                # 检查指定列是否存在
                 # 将能量截面比值都写入数组中
                   data = columns[target_column]
                   C_array.append(data)
                   data_c=columns[target_column-1]
                   Cx_array.append(data_c)
               else:
                   print(f"行 '{line.strip()}' 中没有第{target_column + 1}列的数据")
                   print(f"0")
    except FileNotFoundError:
        print(f"找不到文件: {file_path}")
    except Exception as e:
        print(f"发生错误: {e}")


    E_file_path = E_file_pth  # 打开并读取实验结果
    try:
        with open(E_file_path, 'r') as file:
        # 读取文件的每一行
            lines = file.readlines()

        # 指定要提取的列数（从0开始计数）
            target_column = int(E_column)  # 这里假设你想提取第二列数据
            E_array=[]
            Ex_array=[]
        # 提取指定列的数据并打印出来
            for line in lines:
            # 分割每一行的数据（假设数据之间使用空格或制表符分隔）
               columns = line.strip().split()  # 使用split('\t')来处理制表符分隔的文件
               if len(columns) > target_column:
                # 检查指定列是否存在
                   data = columns[target_column]
                   E_array.append(data)
                   data_e=columns[target_column-1]
                   Ex_array.append(data_e)

               else:
                   print(f"行 '{line.strip()}' 中没有第{target_column + 1}列的数据")
                   print(f"e")
    except FileNotFoundError:
        print(f"找不到文件: {file_path}")
    except Exception as e:
        print(f"发生错误: {e}")
    # 由于上面的C_array,E_array,Cx_array,Ex_array 都是str类型的，需要将其转化为float 类型
    C_array =list(map(float,C_array))
    Cx_array=list(map(float,Cx_array))
    E_array =list(map(float,E_array))
    Ex_array=list(map(float,Ex_array))
    return C_array,Cx_array,E_array,Ex_array

C,Cx,E,Ex=get_data(C_file_pth,C_column,E_file_pth,E_column)  # 得到一组理论计算结果和实验值
# print(E)
# print(C)


## 此时 已经有了一组实验值和理论值，下面开始计算chi^2

chi2,p=stats.chisquare(C,E,) # 将chi^2存入一个变量中
#print("chi^2:",chi2)
print(chi2)

# x=np.arange(36,72,1)
plt.plot(Cx, C,label="Cal",color='b')
plt.plot(Ex,E,label="Exp",color='k',linestyle=':',marker='D',markersize=2)  
plt.legend()  #添加图例
plt.text(60,0.7, f"$ \chi^2$:{round(chi2,5)}", size = 15, alpha = 1)
plt.savefig(fig_info)
