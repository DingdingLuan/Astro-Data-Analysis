import numpy as np
import pandas as pd
from allkindsoffunctions import *
import matplotlib.pyplot as plt

# 对4C+0102进行高斯平滑
c0102=pd.read_excel('/Users/dingding/Desktop/4c+01.02/0102selected.xlsx')
start_time1=c0102['time']
end_time1=c0102['end_time']
fluence1=c0102['flux_0p1_300_gev']
error1=c0102['flux_0p1_300_gev_error']
fluence1=np.array(fluence1)
gaussfluence=Gauss_multiFlat(fluence1,a=1,part=50,n=6)     #倒数第二个是窗口数目，最后一个是平滑阶shu
plt.figure(figsize=(20,10))
plt.plot(start_time1,np.log10(gaussfluence),label='Gauss_4c0102',linewidth=1,color='blue')
plt.plot(start_time1,np.log10(fluence1),label='4c0102',linewidth=0.2,color='r')
plt.xlabel('time[day]')
plt.ylabel('log10_gaussfluence[erg s$^{-1}$ cm$^{-2}$]')
plt.title('4C+0102 Gauss-Flat Light Curve,n=6')
plt.legend(loc='best')
# plt.savefig('/users/dingding/desktop/123.jpg')
# plt.show()

# 对4C+2807进行高斯平滑
c2807=pd.read_excel('/Users/dingding/Desktop/4c+28.07/2807selected.xlsx')
start_time2=c2807['time']
end_time2=c2807['end_time']
fluence2=c2807['flux_0p1_300_gev']
error2=c2807['flux_0p1_300_gev_error']
fluence2=np.array(fluence2)
gaussfluence=Gauss_multiFlat(fluence2,a=1,part=50,n=6)

plt.figure(figsize=(20,10))
plt.plot(start_time2,np.log10(gaussfluence),label='Gauss_4c2807',linewidth=1,color='blue')
plt.plot(start_time2,np.log10(fluence2),label='4c2807',linewidth=0.2,color='r')
plt.xlabel('time[day]')
plt.ylabel('log10_gaussfluence[erg s$^{-1}$ cm$^{-2}$]')
plt.title('4C+2807 Gauss-Flat Light Curve,n=6')
plt.legend(loc='best')
plt.show()