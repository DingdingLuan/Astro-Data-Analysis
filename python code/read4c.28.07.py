import numpy as np
import matplotlib.pyplot as plt
from allkindsoffunctions import *
import pandas as pd

# 读取4C+28.07头文件
c2807=pd.read_excel('/Users/dingding/Desktop/4c+28.07/2807selected.xlsx')
start_time2=c2807['time']
end_time2=c2807['end_time']
fluence2=c2807['flux_0p1_300_gev']
error2=c2807['flux_0p1_300_gev_error']

# 做出4C+28.07原始的光变曲线
plt.figure(figsize=(20,10))
plt.plot(start_time2,np.log10(fluence2),label='4c2807',linewidth=0.2,color='r',marker='o',markerfacecolor='blue',markersize=3)
plt.xlabel('time[day]')
plt.ylabel('log10_fluence[erg s$^{-1}$ cm$^{-2}$]')
plt.title('4C+2807 Light Curve Origin')
plt.legend(loc='best')

# 判断4C+28.07的时间序列的平稳性
Finaltest(fluence2,'fluence2')


# 画出4C+28.07的自相关函数的分布图
self=correlationselffunction(fluence2)
x=np.arange(0,len(self),1)
plt.figure(figsize=(20,10))
plt.scatter(x, self,s=5,alpha=1,marker='*',c='blue',label='Correlation-func of 4c2807')
plt.axhline(0.0, color= 'black')
plt.grid(True, linestyle = "-.", color = "grey", linewidth = "0.2")
plt.title("4C+2807 Self-correlative Function distribution\nSample number=%a" %(len(fluence2)))
plt.xlabel('τ')
plt.ylabel("Correlation Function")
plt.legend(loc='best')
plt.show()

