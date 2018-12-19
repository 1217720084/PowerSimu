"""


"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, spmatrix
import random
import numpy as np

class delays():
    def __init__(self):

        self.n1 = 0
        self._type = 'Delays'
        self._name = 'Delays'
        self._data = {'t0': 50.0, 'tp': 50, 'ts': 0.0, 'delt': 0.0, 'ts_max': 100, 'a': 2, 'b': 0.01,
                      'p': 0.3, 'n': 7, 'm': 10, 'errt': 1e-6, 'T': 50, 'u': 1}
        for arg in self._data:
            self.__dict__[arg] = []
        self.u2 = [1] # 记录数据包是否成功接收
        self.ts1 = [0.0] # 记录随机延时

    def setup(self):

        self.n1 = 1
        for key, value in self._data.items():
            self.__dict__[key].append(value)


    def setx0(self):

        if self.n1 == 0:
            return
        self.T[0] = self.T[0] / 1000
        self.tp[0] = self.tp[0] / 1000
        self.t0[0] = self.t0[0] / 1000
        self.t_drop = [0.0]*self.n1
        self.ts[0] = 0.0
        # self.t_drop = matrix(0.0, (self.n, 1))

    def calt(self, ti_1, ti):

        if self.n1 == 0:
            return

        tp1 = ti % self.T[0]
        self.tp[0] = tp1+self.t_drop[0]
        # # 随机时延
        # for i in range(self.n1):
        #     # self.ts[i] = random.gammavariate(self.a[i], self.b[i])
        #     self.ts[i] = np.random.gamma(self.a[i], self.b[i])
        #     if self.ts[i] > self.ts_max[i]:
        #         self.ts[i] = self.ts_max[i]

        # 判断是否有新的数据包到达
        if int(ti_1/self.T[0]) < int(ti/self.T[0]):
            # self.delt = self.delt - ti
            for i in range(self.n1):
                # 判断是否成功接收到数据包
                q = random.randint(1, self.m[i])
                if q <= self.n[i]:
                    self.t_drop[i] = 0.0
                    self.tp[0] = tp1 + self.t_drop[0]
                    print('成功接收数据包')
                    self.u2.append(1)

                    # 随机时延 （成功接收数据包才算）
                    for i in range(self.n1):
                        # self.ts[i] = random.gammavariate(self.a[i], self.b[i])
                        self.ts[i] = np.random.gamma(self.a[i], self.b[i])
                        if self.ts[i] > self.ts_max[i]:
                            self.ts[i] = self.ts_max[i]
                    self.ts1.append(self.ts[i])
                else:
                    self.t_drop[i] = self.t_drop[i] + self.T[i]
                    self.tp[i] = tp1 + self.t_drop[i]
                    print('不成功接收数据包')
                    self.u2.append(0)
                    self.ts1.append(self.ts[i])
            return (self.tp, self.ts, self.t0)
        else:
            return(self.tp, self.ts, self.t0)