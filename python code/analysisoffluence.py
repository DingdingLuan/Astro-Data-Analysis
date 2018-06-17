import numpy as np
import matplotlib.pyplot as plt
from allkindsoffunctions import *
import random
import pandas as pd
from scipy.stats import kstest
from scipy.stats import ks_2samp


# 读取两个源的原始流量数据
c0102=pd.read_excel('/Users/dingding/Desktop/4c+01.02/0102selected.xlsx')
fluence1=c0102['flux_0p1_300_gev']
start_time1=c0102['time']
c2807=pd.read_excel('/Users/dingding/Desktop/4c+28.07/2807selected.xlsx')
fluence2=c2807['flux_0p1_300_gev']
start_time2=c2807['time']
f3=fluence2


# 对两个源进行随机抽样 然后进行期望与方差估计 抽样标准为完备样本的一半
i=0
mean1array=[]
mean2array=[]
var1array=[]
var2array=[]
k=50
# 取随机操作数为k次，然后通过ADF方法检测序列的平稳性 从而来检验期望与方差
for i in range(k):
    random_fluence1=random.sample(list(fluence1),int(len(fluence1)/2))
    random_fluence2=random.sample(list(fluence2),int(len(fluence2)/2))
    estimted_mean1=np.mean(random_fluence1)
    estimted_mean2=np.mean(random_fluence2)
    estimted_var1=np.var(random_fluence1)
    estimted_var2=np.var(random_fluence2)
    mean1array=np.append(mean1array,estimted_mean1)
    mean2array=np.append(mean2array,estimted_mean2)
    var1array=np.append(var1array,estimted_var1)
    var2array=np.append(var2array,estimted_var2)
    i=i+1
Finaltest(mean1array,'mean1array')
Finaltest(mean2array,'mean2array')
Finaltest(var1array,'var1array')
Finaltest(var2array,'var2array')


# 对两个源的流量分别进行K-S检验，判断其分布的正太性
kstest_0102=kstest(fluence1, 'norm')
kstest_2807=kstest(fluence2,'norm')
if kstest_0102[1]<=0.05:
    print('4c+0102 is not normal distribution\npvalue=%a' %(kstest_0102[1]))
else:
    print('4c+0102 is normal distribution\npvalue=%a' %(kstest_0102[1]))

if kstest_2807[1]<=0.05:
    print('4c+2807 is not normal distribution\npvalue=%a' %(kstest_2807[1]))
else:
    print('4c+2807 is normal distribution\npvalue=%a' %(kstest_2807[1]))


# 对两个源进行z-score标准化处理，并生成新的流量序列
new0102fluence=z_score(fluence1)
new2807fluence=z_score(fluence2)
a=ks_2samp(new0102fluence,new2807fluence)
if a[1]<=0.05:
    print('After z-score standardization,4c+0102 and 4c+2807 are not in same distribution,pvalue=%a' %(a[1]))
else:
    print('After z-score standardization,4c+0102 and 4c+2807 are  in same distribution,pvalue=%a' % (a[1]))


# 绘制p-p和q-q图来检测两个源原始数据的正太性
# 求累积分布的数据值,并画出p-p图：
# 4c+0102:
total_fluence=np.sum(fluence1)
fluence1=fluence1[np.argsort(fluence1)]
fluence1=np.array(fluence1)
i=1
for i in range(len(start_time1)):
    fluence1[i]=fluence1[i]+fluence1[i-1]
    i=i+1
y1=fluence1
y1=fluence1/total_fluence
x1=start_time1-np.min(start_time1)
x2=np.arange(0, len(start_time1), 1)
y2=np.random.normal(0,1,len(start_time1))
y2=np.abs(y2)
total_normal=np.sum(y2)
i=1
for i in range(len(y2)):
    y2[i]=y2[i]+y2[i-1]
    i=i+1
y2=y2/total_normal
a=np.arange(0,1,0.01)
b=np.arange(0,1,0.01)
plt.plot(y2,y1)
plt.title('p-p plot for 4C+0102')
plt.xlabel('Theoritical normal Data')
plt.ylabel('Observed Data')
plt.plot(a,b,linewidth=2,color='r')
plt.legend(labels=['Actual result','Perfect result'],loc='best')
plt.show()
# 4c+2807:
total_fluence=np.sum(fluence2)
fluence2=fluence2[np.argsort(fluence2)]
fluence2=np.array(fluence2)
i=1
for i in range(len(start_time2)):
    fluence2[i]=fluence2[i]+fluence2[i-1]
    i=i+1
y1=fluence2
y1=fluence2/total_fluence
x1=start_time2-np.min(start_time2)
x2=np.arange(0, len(start_time2), 1)
y2=np.random.normal(0,1,len(start_time2))
y2=np.abs(y2)
total_normal=np.sum(y2)
i=1
for i in range(len(y2)):
    y2[i]=y2[i]+y2[i-1]
    i=i+1
y2=y2/total_normal
a=np.arange(0,1,0.01)
b=np.arange(0,1,0.01)
plt.plot(y2,y1)
plt.title('p-p plot for 4C+2807')
plt.xlabel('Theoritical normal Data')
plt.ylabel('Observed Data')
plt.plot(a,b,linewidth=2,color='r')
plt.legend(labels=['Actual result','Perfect result'],loc='best')
plt.show()


# 求累积分布的数据值,并画出q-q图：
# 4c+0102:
total_fluence=np.sum(fluence1)
fluence1=fluence1[np.argsort(fluence1)]
fluence1=np.array(fluence1)
i=1
for i in range(len(start_time1)):
    fluence1[i]=fluence1[i]+fluence1[i-1]
    i=i+1
y1=fluence1
x1=start_time1-np.min(start_time1)
x2=np.arange(0, len(start_time1), 1)
y2=np.random.normal(0,1,len(start_time1))
y2=np.abs(y2)
total_normal=np.sum(y2)
i=1
for i in range(len(y2)):
    y2[i]=y2[i]+y2[i-1]
    i=i+1
a=np.arange(0,np.max(y2),np.max(y2)/len(y2))
b=np.arange(0,np.max(y1),np.max(y1)/len(y2))
plt.plot(y2,y1)
plt.title('q-q plot for 4C+0102')
plt.xlabel('Theoritical random normal data')
plt.ylabel('Observed fluence Data')
plt.plot(a,b,linewidth=2,color='r')
plt.legend(labels=['Actual result','Perfect result'],loc='best')
plt.show()

# 4c+2807:
total_fluence=np.sum(f3)
fluence2=f3[np.argsort(f3)]
fluence2=np.array(fluence2)
i=1
for i in range(len(start_time2)):
    fluence2[i]=fluence2[i]+fluence2[i-1]
    i=i+1
y1=fluence2
x1=start_time2-np.min(start_time2)
x2=np.arange(0, len(start_time2), 1)
y2=np.random.normal(0,1,len(start_time2))
y2=np.abs(y2)
total_normal=np.sum(y2)
i=1
for i in range(len(y2)):
    y2[i]=y2[i]+y2[i-1]
    i=i+1
a=np.arange(0,np.max(y2),np.max(y2)/len(y2))
b=np.arange(0,np.max(y1),np.max(y1)/len(y2))
plt.plot(y2,y1)
plt.title('q-q plot for 4C+2807')
plt.xlabel('Theoritical random normal data')
plt.ylabel('Observed fluence Data')
plt.plot(a,b,linewidth=2,color='r')
plt.legend(labels=['Actual result','Perfect result'],loc='best')
plt.show()