import pandas as pd
import numpy as np
from scipy import  stats
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.graphics.api import qqplot
from allkindsoffunctions import *


df = pd.read_csv('/users/dingding/desktop/4c+28.07/2807arimaselectedtest.csv', encoding='utf-8', index_col='date')
df.index = pd.to_datetime(df.index)  # 将字符串索引转换成时间索引
ts = df['flux']  # 生成pd.Series对象


dta=ts
dta.index=df.index
# dta.plot(figsize=(12,8))
# plt.show()

# fig = plt.figure(figsize=(12,8))
# ax1= fig.add_subplot(111)
# diff1 = dta.diff(1)
# diff1.plot(ax=ax1)
# plt.show()

# diff1= dta.diff(1)
# fig = plt.figure(figsize=(12,8))
# ax1=fig.add_subplot(211)
# fig = sm.graphics.tsa.plot_acf(dta,lags=40,ax=ax1)
# ax2 = fig.add_subplot(212)
# fig = sm.graphics.tsa.plot_pacf(dta,lags=40,ax=ax2)
# # plt.show()

arma_mod80 = sm.tsa.ARMA(dta,(8,0)).fit()
# print(arma_mod80.aic,arma_mod80.bic,arma_mod80.hqic)
resid = arma_mod80.resid
# fig = plt.figure(figsize=(12,8))
# ax1 = fig.add_subplot(211)
# fig = sm.graphics.tsa.plot_acf(resid.values.squeeze(), lags=40, ax=ax1)
# ax2 = fig.add_subplot(212)
# fig = sm.graphics.tsa.plot_pacf(resid, lags=40, ax=ax2)
# plt.show()


# print(stats.normaltest(resid))
# fig = plt.figure(figsize=(12,8))
# ax = fig.add_subplot(111)
# fig = qqplot(resid, line='q', ax=ax, fit=True)
# plt.show()


r,q,p = sm.tsa.acf(resid.values.squeeze(), qstat=True)
data = np.c_[range(1,41), r[1:], q, p]
table = pd.DataFrame(data, columns=['lag', "AC", "Q", "Prob(>Q)"])
# print(table.set_index('lag'))
#
predict_dta = arma_mod80.predict('2016-12-27', '2017-1-7', dynamic=True)
print(predict_dta)

fig, ax = plt.subplots(figsize=(12, 8))
ax = dta.ix['2015':].plot(ax=ax)
fig = arma_mod80.plot_predict('2016-12-27', '2017-1-7', dynamic=True, ax=ax, plot_insample=False)
plt.show()