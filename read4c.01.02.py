import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from allkindsoffunctions import *

# 读取4C+01.02头文件
c0102=pd.read_excel('/Users/dingding/Desktop/4c+01.02/0102selected.xlsx')
start_time1=c0102['time']
end_time1=c0102['end_time']
fluence1=c0102['flux_0p1_300_gev']
error1=c0102['flux_0p1_300_gev_error']

# 做出4C+01.02原始的光变曲线
plt.figure(figsize=(20,10))
plt.plot(start_time1,np.log10(fluence1),label='4c0102',linewidth=0.2,color='r',marker='o',markerfacecolor='blue',markersize=3)
plt.xlabel('time[day]')
plt.ylabel('log10_fluence[erg s$^{-1}$ cm$^{-2}$]')
plt.title('4C+0102 Light Curve Origin')
plt.legend(loc='best')
plt.savefig('/users/dingding/desktop/2.jpg')

# 判断4C+01.02的时间序列的平稳性
Finaltest(fluence1,'fluence1')


# 画出4C+01.02的自相关函数的分布图
self=correlationselffunction(fluence1)
x=np.arange(0,len(self),1)
# ax=plt.figure(figsize=(20,10)).add_subplot(111)
# ax.spines['bottom'].set_position(('data', 0))        #改变坐标轴位置
plt.figure(figsize=(20,10))
plt.scatter(x, self,s=5,alpha=1,marker='*',c='blue',label='Correlation-func of 4c0102')
plt.axhline(0.0, color= 'black')
plt.grid(True, linestyle = "-.", color = "grey", linewidth = "0.2")
# plt.plot(x,np.zeros(len(self)),linewidth=0.2,color='black')
plt.legend(loc='best')
plt.title("4C+0102 Self-correlative Function distribution\nSample number=%a" %(len(fluence1)))
plt.xlabel("τ")
plt.ylabel("Correlation Function")
plt.savefig('/users/dingding/desktop/1.jpg')
plt.show()




