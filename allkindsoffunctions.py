# 运用ADF方法检测时间序列平稳性(原假设为非平稳),若其p值可以拒绝原假设 即其为平稳时间序列
def testStationarity(ts):
    import pandas as pd
    from statsmodels.tsa.stattools import adfuller
    dftest = adfuller(ts)
    dfoutput = pd.Series(dftest[0:4], index=['Test Statistic','p-value','#Lags Used','Number of Observations Used'])
    for key,value in dftest[4].items():
        dfoutput['Critical Value (%s)'%key] = value
    return dfoutput



# 求序列的自相关函数（输出为相关系数数列）
def correlationselffunction(ts):
    import numpy as np
    k=0
    tao=[]
    miu=np.mean(ts)
    for k in range(len(ts)):
        up=[]
        down=[]
        t=0
        for t in range(len(ts)-k):
            up=np.append(up,(ts[t]-miu)*(ts[t+k]-miu))
            t=t+1
        t=0
        for t in range(len(ts)):
            down=np.append(down,(ts[t]-miu)**2)
            t=t+1
        tao=np.append(tao,np.sum(up)/np.sum(down))
        k=k+1
    return tao



# 求序列的互相关函数（输出为相关系数数列）
def correlationothersfunction(ts1,ts2):     #这里要求ts1与ts2是用维度的时间序列(需要预先筛选同一时间段的两个源流量的数列值进行计算)
    import numpy as np
    if len(ts1)!=len(ts2):
        print('Warning:Array Dimensions are not matching!!!Please check your data..')
        exit()
    k=0
    tao=[]
    for k in range(len(ts1)):
        up=[]
        down1=[]
        down2=[]
        miu1 = np.mean(ts1[k:])
        miu2 = np.mean(ts2[k:])
        t=0
        for t in range(len(ts1)-k):
            up=np.append(up,(ts1[t+k]-miu1)*(ts2[t+k]-miu2))
            down1=np.append(down1,(ts1[t+k]-miu1)**2)
            down2=np.append(down2,(ts2[t+k]-miu2)**2)
        ch1=np.sum(up)
        ch2=np.sum(down1)
        ch3=np.sum(down2)
        stao=ch1/np.sqrt(ch2*ch3)
        tao=np.append(tao,stao)
        k=k+1
    return tao



# 判断时间序列的平稳性
def Finaltest(ts,tsname):
    a=testStationarity(ts)
    if a[1]<0.99:
        print('Time Series %s is Stable，P-value=%a' %(tsname,a[1]))
    elif a[1]>=0.99:
        print('Time series %s is Unstable，p-value=%a' %(tsname,a[1]))
    return



# z-score标准化处理方法,并生成新的数据
def z_score(ts):
    import numpy as np
    mean=np.mean(ts)
    std=np.std(ts)
    new=np.zeros_like(ts)
    i=0
    for i in range(len(ts)):
        new[i]=(ts[i]-mean)/std
        i=i+1
    return new



# 定义一个高斯概率密度分布函数p(x)
def p(x,miu,a):
    import numpy as np
    p=1/(np.sqrt(2*np.pi)*a)*np.exp(-(x-miu)**2/(2*a**2))
    return p



# 定义一个高斯权函数平滑过程计算函数
def gaussflat(ts,a,l):
    import numpy as np
    unitnumber=int(len(ts)/l)                   #将高斯过平滑程定义为l个窗口
    leftseriesnumber=len(ts)-l*unitnumber
    i=unitnumber
    n=int(leftseriesnumber/i)
    c=n*i
    if (i%2)==0:                                  #判断采样值的奇偶性
        miu=int(i/2)
    else:
        miu=int(i/2+1)
    j=0
    window=ts[:i]
    x_prime=np.zeros_like(ts)
    for j in range((l-1)*unitnumber+c):
        k=0
        m=[]
        s=[]
        for k in range(unitnumber):
            mm=np.abs(p(k,miu,a)-p(k+1,miu,a))*np.exp(-(k-miu)**2/(2*a**2))
            ss=mm*window[k]
            m=np.append(m,mm)
            s=np.append(s,ss)
            k=k+1
        M=np.sum(m)
        S=np.sum(s)
        x_prime[miu+j]=S/M
        j=j+1
        window=ts[j:i+j]
    return x_prime



# 定义优化后的高斯平滑函数，使得其首尾的值可以进行相连从而让全序列都进行优化
def optimizedgaussflat(ts, a, part):
    import numpy as np
    unitnumber=int(len(ts)/part)+1
    difference=unitnumber*part-len(ts)
    newts=np.append(ts[len(ts)-int(unitnumber/2):len(ts)],ts)
    newts=np.append(newts,ts[:difference])                            #构造新的ghost序列第一步：补全整体性
    newts=np.append(newts,ts[difference:int(unitnumber/2)])           #补全ghost新数列的尾巴 添加一个半宽高斯通道
    if unitnumber%2==0:
        windows=newts[:unitnumber+1]                                  #类似upwind算法嘿嘿嘿
        lens=unitnumber+1
    else:
        lens=unitnumber
        windows=newts[:unitnumber]                                    #对窗口数目的奇偶性质进行讨论
    miu=int(unitnumber/2+1)
    aft=miu+len(ts)

    x_prime=newts
    for j in range(len(ts)):
        m=[]
        s=[]
        for k in range(lens):
            mm=np.abs(p(k,miu,a)-p(k+1,miu,a))*np.exp(-(k-miu)**2/(2*a**2))
            ss=mm*windows[k]
            m=np.append(m,mm)
            s=np.append(s,ss)
        M=np.sum(m)
        S=np.sum(s)
        x_prime[miu+j]=S/M
        j=j+1
        if unitnumber%2==0:
            windows=newts[j:unitnumber+1+j]
        else:
            windows=newts[j:unitnumber+j]
    return x_prime[miu:aft]



#添加一个参数n，来控制平滑的次数
def Gauss_multiFlat(ts,a,part,n):
    new=optimizedgaussflat(ts,a,part)
    for i in range(n):
        new=optimizedgaussflat(new,a,part)
    return new







import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf

# 移动平均图
def draw_trend(timeSeries, size):
    f = plt.figure(facecolor='white')
    # 对size个数据进行移动平均
    rol_mean = timeSeries.rolling(window=size).mean()
    # 对size个数据进行加权移动平均
    rol_weighted_mean = pd.ewma(timeSeries, span=size)

    timeSeries.plot(color='blue', label='Original')
    rol_mean.plot(color='red', label='Rolling Mean')
    rol_weighted_mean.plot(color='black', label='Weighted Rolling Mean')
    plt.legend(loc='best')
    plt.title('Rolling Mean')
    plt.show()



def draw_ts(timeSeries):
    f = plt.figure(facecolor='white')
    timeSeries.plot(color='blue')
    plt.show()



# 自相关和偏相关图，默认阶数为31阶
def draw_acf_pacf(ts, lags=31):
    f = plt.figure(facecolor='white')
    ax1 = f.add_subplot(211)
    plot_acf(ts, lags=31, ax=ax1)
    ax2 = f.add_subplot(212)
    plot_pacf(ts, lags=31, ax=ax2)
    plt.show()



# 差分操作
def diff_ts(ts, d):
    global shift_ts_list
    #  动态预测第二日的值时所需要的差分序列
    global last_data_shift_list
    shift_ts_list = []
    last_data_shift_list = []
    tmp_ts = ts
    for i in d:
        last_data_shift_list.append(tmp_ts[-i])
        print(last_data_shift_list)
        shift_ts = tmp_ts.shift(i)
        shift_ts_list.append(shift_ts)
        tmp_ts = tmp_ts - shift_ts
    tmp_ts.dropna(inplace=True)
    return tmp_ts

# 还原操作
def predict_diff_recover(predict_value, d):
    if isinstance(predict_value, float):
        tmp_data = predict_value
        for i in range(len(d)):
            tmp_data = tmp_data + last_data_shift_list[-i-1]
    elif isinstance(predict_value, np.ndarray):
        tmp_data = predict_value[0]
        for i in range(len(d)):
            tmp_data = tmp_data + last_data_shift_list[-i-1]
    else:
        tmp_data = predict_value
        for i in range(len(d)):
            try:
                tmp_data = tmp_data.add(shift_ts_list[-i-1])
            except:
                raise ValueError('What you input is not pd.Series type!')
        tmp_data.dropna(inplace=True)
    return tmp_data



from dateutil.relativedelta import relativedelta
def _add_new_data(ts, dat, type='day'):
    if type == 'day':
        new_index = ts.index[-1] + relativedelta(days=1)
    elif type == 'month':
        new_index = ts.index[-1] + relativedelta(months=1)
        ts[new_index] = dat

def add_today_data(model, ts,  data, d, type='day'):
    _add_new_data(ts, data, type)  # 为原始序列添加数据
    # 为滞后序列添加新值
    d_ts = diff_ts(ts, d)
    model.add_today_data(d_ts[-1], type)

def forecast_next_day_data(model, type='day'):
    if model == None:
        raise ValueError('No model fit before')
    fc = model.forecast_next_day_value(type)
    return predict_diff_recover(fc, [12, 1])

def proper_model(data_ts, maxLag):
    import numpy as np
    from statsmodels.tsa.arima_model import ARMA
    import sys
    init_bic = sys.maxint
    init_p = 0
    init_q = 0
    init_properModel = None
    for p in np.arange(maxLag):
        for q in np.arange(maxLag):
            model = ARMA(data_ts, order=(p, q))
            try:
                results_ARMA = model.fit(disp=-1, method='css')
            except:
                continue
            bic = results_ARMA.bic
            if bic < init_bic:
                init_p = p
                init_q = q
                init_properModel = results_ARMA
                init_bic = bic
    return init_bic, init_p, init_q, init_properModel


