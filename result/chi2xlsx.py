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
import openpyxl

C_file_pth = sys.argv[1]  # cal文件路径
E_file_pth=sys.argv[2]
cross_fig=sys.argv[3]
barri_fig=sys.argv[4]
#### 这个函数用于读取理论计算结果和实验结果，并将其返回为4个list类型的数据,分别是实验和理论的有效能量以及截面比

def get_data(C_file_pth):
    C_file_path = C_file_pth  # 打开并读取计算文件
    try:
       with open(C_file_path, 'r') as file:
        # 读取文件的每一行
            lines = file.readlines()

        # 指定要提取的列数（从0开始计数）
            target_column = int(2)  # 这里假设你想提取第二列数据
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
            
            # 定义变量,读取了贾老师beta4=-0.01,0.03,0.07的数据，可以一起画出来
            E_array=[] #B
            Ex_array=[]#A
            Jia_x=[]   #D
            Jia_001=[] #E
            Jia_003=[] #F
            Jia_007=[] #G
            barrierX_array=[] #H
            barrier_array=[]  #I
            barrier_Jia_x=[]  #K
            barrier_Jia_001=[]#L
            barrier_Jia_003=[]#M
            barrier_Jia_007=[]#N
        # 提取指定列的数据并打印出来，重复了两次，因为截面和势垒分布的行数不一样，放一起会出现问题
          
        #下面开始读取数据
            workbook = openpyxl.load_workbook(E_file_pth)
            sheet = workbook.active
            for cell in sheet['A']:
                Ex_array.append(cell.value)
        
            for cell in sheet['B']:
                E_array.append(cell.value)

            for cell in sheet['D']:
                Jia_x.append(cell.value)

            for cell in sheet['E']:
                Jia_001.append(cell.value)

            for cell in sheet['F']:
                Jia_003.append(cell.value)

            for cell in sheet['G']:
                Jia_007.append(cell.value)

            for cell in sheet['H']:
                barrierX_array.append(cell.value)

            for cell in sheet['I']:
                barrier_array.append(cell.value)

            for cell in sheet['K']:
                barrier_Jia_x.append(cell.value)

            for cell in sheet['L']:
                barrier_Jia_001.append(cell.value)

            for cell in sheet['M']:
                barrier_Jia_003.append(cell.value)

            for cell in sheet['N']:
                barrier_Jia_007.append(cell.value)

    except FileNotFoundError:
        print(f"找不到文件: {file_path}")
    except Exception as e:
        print(f"发生错误: {e}")
    # 由于上面的C_array,E_array,Cx_array,Ex_array 都是str类型的，需要将其转化为float 类型
    C_array =list(map(float,C_array))
    Cx_array=list(map(float,Cx_array))

    return C_array,Cx_array,Jia_x,Jia_001,Jia_003,Jia_007,E_array,Ex_array,barrierX_array,barrier_array,barrier_Jia_x,barrier_Jia_001,barrier_Jia_003,barrier_Jia_007


C,Cx,Jia_x,Jia_001,Jia_003,Jia_007,E,Ex,Bx,B,barrier_Jia_x,barrier_Jia_001,barrier_Jia_003,barrier_Jia_007=get_data(C_file_pth)  # 得到一组理论计算>结果和实验值
# print(len(E))
# print(len(C))


#####################################
#chi2,p=stats.chisquare(C,E,) # 将chi^2存入一个变量中

#print("chi^2:",chi2)
#print(chi2)
#########################################

# 截面分布
'''
######################################################################
#plt.text(60,0.7, f"$ \chi^2$:{round(chi2,5)}", size = 15, alpha = 1)
#plt.text(60,0.7, f"$ \chi^2$:{round(chi2,5)}", size = 15, alpha = 1)
fig0=plt.figure()
plt.plot(Cx, C,label="Cal",color='r')
plt.plot(Jia_x,Jia_001,label='-0.01')
plt.plot(Jia_x,Jia_003,label='0.03')
plt.plot(Jia_x,Jia_007,label='0.07')
plt.plot(Ex,E,label="Exp",color='k',linestyle=':',marker='D',markersize=2)
plt.legend()  #添加图例
#plt.text(60,0.6, f"$ \chi^2$:{round(chi2,5)}", size = 15, alpha = 1)
plt.savefig(cross_fig)
###################################################################
'''


# 下面用中心差分方法计算势垒分布,实验值已经获得，直接读取就可
barrier_C = []
barrier_Cx = []
# barrier_E = []

for i in range(0,len(C)-2):
    barrier_C.append(-(C[i+2]-C[i])/(Cx[i+2]-Cx[i]))
    barrier_Cx.append(Cx[i]+0.5*(Cx[i+2]-Cx[i])) ##attiention!!!!! here need to be change
# for j in range(0,len(E)-2):
#    barrier_E.append(-(E[j+2]-E[j])/(Ex[j+2]-Ex[j]))

#print(len(B))
#print(len(Bx))

fig = plt.figure()
#设置坐标轴范围
plt.xlim((52, 70))
plt.ylim((0, 0.2))

#设置坐标轴刻度
my_x_ticks = np.arange(52,70,2)
my_y_ticks = np.arange(0,0.2,0.04)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)

#生成网格
plt.grid()

plt.plot(barrier_Cx, barrier_C,label="Cal",color='r')
# plt.plot(barrier_Jia_x,barrier_Jia_001,label='-0.01')
# plt.plot(barrier_Jia_x,barrier_Jia_003,label='0.03')
#plt.plot(barrier_Jia_x,barrier_Jia_007,label='0.07')
plt.plot(Bx,B,label="Exp",color='k',linestyle=':',marker='D',markersize=2)
plt.legend()  #添加图例
#plt.text(60,0.7, f"$ \chi^2$:{round(chi2,5)}", size = 15, alpha = 1)
plt.savefig(barri_fig)














