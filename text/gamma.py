import matplotlib.pyplot as plt
import scipy.special as sps
import numpy as np
from pylab import mpl
mpl.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号

shape, scale = 2, 0.01 # mean=4, std=2*sqrt(2)
s = np.random.gamma(shape, scale, 1000)
print(np.mean(s))
# print(np.median(s))
plt.subplot(211)
[count, bins, ignored] = plt.hist(s, 50, density='True')
plt.title('a=2,b=0.01')
plt.ylabel(u'抽样结果统计')
# [count, bins, ignored] = plt.hist(s, 50)
y = bins**(shape-1)*(np.exp(-bins/scale) /
                      (sps.gamma(shape)*scale**shape))
plt.subplot(212)
plt.plot(bins[:], y[:], linewidth=2, color='r')
# plt.axis([0.0, 0.2, 0, 2])
plt.ylabel(u'概率密度分布函数')
plt.xlabel(u't/s')
plt.show()