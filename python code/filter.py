import pandas as pd
import numpy as np
from allkindsoffunctions import *
import matplotlib.pyplot as plt

c0102=pd.read_csv('/Users/dingding/Desktop/4c+01.02/0102filter.csv')
c2807=pd.read_csv('/Users/dingding/Desktop/4c+28.07/2807filter.csv')
time0102=c0102['time']
fluence0102=c0102['fluence']
time2807=c2807['time']
fluence2807=c2807['fluence']
time2807=np.array(time2807)
time0102=np.array(time0102)

# 筛选
i=0
timeselected0102=[]
fluenceselected0102=[]
for i in range(len(time0102)):
    if time0102[i] in time2807:
        timeselected0102=np.append(timeselected0102,time0102[i])
        fluenceselected0102=np.append(fluenceselected0102,fluence0102[i])
    i=i+1

# 重组原数组格式
i=0
timeselected2807=[]
fluenceselected2807=[]
for i in range(len(time2807)):
    if time2807[i] in timeselected0102:
        timeselected2807=np.append(timeselected2807,time2807[i])
        fluenceselected2807=np.append(fluenceselected2807,fluence2807[i])
    i=i+1

# 剪裁数据 控制大小为500
time1=timeselected0102[:500]
time2=timeselected2807[:500]
flu1=fluenceselected0102[:500]
flu2=fluenceselected2807[:500]

# 求根据时间分布所筛选的两个新样本数据的互相关函数,并画出图像：
cf=correlationothersfunction(flu1,flu2)
x=np.arange(0,len(cf),1)
plt.figure(figsize=(20,10))
plt.scatter(x, cf,s=10,alpha=1,marker='*',c='green',label='Interrelation-func of 4c0102&4c2807')
plt.axhline(0.0, color= 'black')
plt.grid(True, linestyle = "-.", color = "grey", linewidth = "0.2")
plt.title("4C+01.02&4C+29.07 Interrelation Function distribution")
plt.legend(loc='best')
plt.xlabel("τ")
plt.ylabel("Interrelation Functionn")
plt.show()
