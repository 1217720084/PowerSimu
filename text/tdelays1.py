# -*- coding: utf-8 -*-
from devices.delays import delays
import numpy as np
import matplotlib.pyplot as plt
from pylab import mpl
mpl.rcParams['font.sans-serif'] = ['SimHei']

Delays = delays()
Delays.setup()
Delays.setx0()

ti_1 = 0.0
t = 0.0
tend = 2
tp1 = []
ts1 = []
tstart = 0.0
ta = [0.0]
T1 = 0.02 #仿真周期
while t < tend:

    [tp, ts, t0] = Delays.calt(ti_1, t)
    print(tp, ts, t0)
    tp1.append(tp[0])
    ts1.append(ts[0])
    ti_1 = t
    t = t + T1
    ta.append(t)

n = len(Delays.u2)
plt.figure(1)
plt.subplot(211)
for i in range(n):
    if Delays.u2[i] == 1:
        x = np.linspace(i*Delays.T[0], (i+1)*Delays.T[0], endpoint='True')
        y = x-i*Delays.T[0]
        tstart = i*Delays.T[0]
        plt.plot(x, y, 'k')

        if i < n-1:
            if Delays.u2[i+1] == 1:
                ymax = max(y)
                y = np.linspace(0.0, ymax, 100, endpoint='True')
                x = [(i+1)*Delays.T[0]]*100
                plt.plot(x, y, 'k')

    else:
        x = np.linspace(i*Delays.T[0], (i+1)*Delays.T[0], endpoint='True')
        y = x - tstart
        plt.plot(x, y, 'k')
        if i < n-1:
            if Delays.u2[i+1] == 1:
                ymax = max(y)
                y = np.linspace(0.0, ymax, 100, endpoint='True')
                x = [(i+1)*Delays.T[0]]*100
                plt.plot(x, y, 'k')
plt.ylabel(u'准周期时延tp(t)/s')
plt.xlabel(u'时间t/s')
plt.axis([0, 2, 0, 0.2])

# plt.figure(2)
plt.subplot(212)
n = len(Delays.u2)
print(Delays.u2)
print(Delays.ts1)
for i in range(n-1):
    if Delays.u2[i+1] == 1:
        x = [(i+1)*Delays.T[0]]*100
        y = np.linspace(Delays.ts1[i], Delays.ts1[i+1], 100, endpoint='True')
        plt.plot(x, y, 'k')

        x = np.linspace((i+1)*Delays.T[0], (i+2)*Delays.T[0], 100, endpoint='True')
        y = [Delays.ts1[i+1]]*100
        plt.plot(x, y, 'k')
    else:
        x = np.linspace((i + 1) * Delays.T[0], (i + 2) * Delays.T[0], 100, endpoint='True')
        y = [Delays.ts1[i + 1]] * 100
        plt.plot(x, y, 'k')
plt.axis([0, 2, 0, 0.08])
plt.ylabel(u'随机时延ts(t)/s')
plt.xlabel(u'时间t/s')
plt.show()

