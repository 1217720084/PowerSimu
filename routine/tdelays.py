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
while t < tend:

    [tp, ts, t0] = Delays.calt(ti_1, t)
    print(tp, ts, t0)
    tp1.append(tp[0])
    ts1.append(ts[0])
    ti_1 = t
    t = t + 0.02

n = len(Delays.u2)
plt.figure(1)
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

plt.figure(2)
t = np.linspace(0.0, 2, 100)
nx = len(t)
for i in range(nx):
    x = np.linspace(i*Delays.T[0], (i+1)*Delays.T[0], 100)
    y = [ts1[i]]*100
    plt.plot(x, y, 'k')
    if i < nx - 1:
        if ts1[i] >= ts1[i+1]:
            y = np.linspace(ts1[i+1], ts1[i], 100, endpoint='True')
            x = [(i + 1) * Delays.T[0]] * 100
            plt.plot(x, y, 'k')
        else:
            y = np.linspace(ts1[i], ts1[i+1], 100, endpoint='True')
            x = [(i + 1) * Delays.T[0]] * 100
            plt.plot(x, y, 'k')

plt.ylabel(u'随机时延ts(t)/s')
plt.xlabel(u'时间t/s')
plt.axis([0, 2, 0, 0.07])
plt.show()
