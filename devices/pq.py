"""

"""
from devices.base_device import base_device
import system
import random
import numpy as np
from cvxopt.base import matrix, spmatrix, cos, sin, sparse, mul,exp
import time

class pq(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pl': 1, 'bus': None, 'Ql': 1,  'Vmax': 1.2, 'Vmin': 0.95, 'z': 1, 'gen': 0, 'shunt': 0})
        self._type = 'PQ'
        self._name = 'PQ'
        self._bus = {'bus': ['a', 'v']}
        self.n = 0
        self._algebs = ['V0', 'Va']
        self._params.extend(['Pl', 'Ql',  'Vmax', 'Vmin'])
        self._powers = ['Pl', 'Ql']
        self._voltages = ['Vmax', 'Vmin']
        self.properties.update({'gcall': True, 'Gycall': True})

    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.Ql):
            dae.g[key] += value
        for key, value in zip(self.a, self.Pl):
            dae.g[key] += value

    def gcall(self):

        if self.n == 0:
            return
        for i in range(system.PQ.n):
            system.DAE.g[self.a[i]] += self.Pl[i]
            system.DAE.g[self.v[i]] += self.Ql[i]
        # 判断是否PQ负荷转为恒阻抗模型

        if system.Settings.forcepq:
            return

        for i in range(self.n):
            if system.DAE.y[self.v[i]] < self.Vmin[i] and self.z[i] and self.u[i] or self.shunt[i]:
                system.DAE.g[self.a[i]] = system.DAE.g[self.a[i]] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2\
                                                                                 / self.Vmin[i] / self.Vmin[i]
                system.DAE.g[self.v[i]] = system.DAE.g[self.v[i]] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]] ** 2 / \
                                                                         self.Vmin[i] / self.Vmin[i]
            if system.DAE.y[self.v[i]] > self.Vmax[i] and self.z[i] and self.u[i]:

                # print('pq达到vmax')
                system.DAE.g[self.a[i]] = system.DAE.g[self.a[i]] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2 /self.Vmax[i] / self.Vmax[i]

                system.DAE.g[self.v[i]] = system.DAE.g[self.v[i]] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]]**2 / self.Vmax[i] / self.Vmax[i]

        # try:
        #     for i in range(self.n):
        #         if system.DAE.y[self.v[i]] < self.Vmin[i] and self.z[i] and self.u[i] or self.shunt[i]:
        #             system.DAE.g[self.a[i]] = system.DAE.g[self.a[i]] - self.Pl[i] + self.Pl[i] * system.DAE.y[
        #                                                                                               self.v[i]] ** 2 \
        #                                                                              / self.Vmin[i] / self.Vmin[i]
        #             system.DAE.g[self.v[i]] = system.DAE.g[self.v[i]] - self.Ql[i] + self.Ql[i] * system.DAE.y[
        #                                                                                               self.v[i]] ** 2 / \
        #                                                                              self.Vmin[i] / self.Vmin[i]
        #         if system.DAE.y[self.v[i]] > self.Vmax[i] and self.z[i] and self.u[i]:
        #             # print('pq达到vmax')
        #             system.DAE.g[self.a[i]] = system.DAE.g[self.a[i]] - self.Pl[i] + self.Pl[i] * system.DAE.y[
        #                                                                                               self.v[i]] ** 2 / \
        #                                                                              self.Vmax[i] / self.Vmax[i]
        #
        #             system.DAE.g[self.v[i]] = system.DAE.g[self.v[i]] - self.Ql[i] + self.Ql[i] * system.DAE.y[
        #                                                                                               self.v[i]] ** 2 / \
        #                                                                              self.Vmax[i] / self.Vmax[i]
        # except:
        #
        #     print( system.DAE.y[self.v])

                    # a = []
        # b = []
        # for i in range(self.n):
        #     if system.DAE.y[self.v[i]] < self.Vmin[i]:
        #
        #         print('pq达到vmin')
        #         system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]
        #
        #         system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]
        #
        #     if system.DAE.y[self.v[i]] > self.Vmax[i]:
        #
        #         print('pq达到vmax')
        #         system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2 /self.Vmin[i] / self.Vmin[i]
        #
        #         system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]



    def Gycall(self):

        if self.n == 0:
            return

        for i in range(self.n):
            if system.DAE.y[self.v[i]] < self.Vmin[i] and self.z[i] and self.u[i] or self.shunt[i]:
                system.DAE.Gy[self.a[i], self.v[i]] = system.DAE.Gy[self.a[i], self.v[i]] + 2*self.Pl[i]*system.DAE.y[self.v[i]]/self.Vmin[i]/self.Vmin[i]
                system.DAE.Gy[self.v[i], self.v[i]] = system.DAE.Gy[self.v[i], self.v[i]] + 2 * self.Ql[i] * system.DAE.y[self.v[i]] / \
                                                                                            self.Vmin[i] / self.Vmin[i]
            if system.DAE.y[self.v[i]] > self.Vmax[i] and self.z[i] and self.u[i]:

                # print('pq达到vmax')
                system.DAE.Gy[self.a[i], self.v[i]] = system.DAE.Gy[self.a[i], self.v[i]] + 2 * self.Pl[i] * \
                                                                                            system.DAE.y[self.v[i]] / \
                                                                                            self.Vmax[i] / self.Vmax[i]
                system.DAE.Gy[self.v[i], self.v[i]] = system.DAE.Gy[self.v[i], self.v[i]] + 2 * self.Ql[i] * \
                                                                                            system.DAE.y[self.v[i]] / \
                                                                                            self.Vmax[i] / self.Vmax[i]
    def pqshunt(self):

        if self.n == 0:
            return

        if not system.Settings.pq2z:
            return

        for i in range(self.n):
            if self.gen[i] == 0:
                self.shunt[i] = 1
            else:
                self.shunt[i] = 0
            if self.shunt[i] == 1:
                self.Vmin[i] = system.DAE.y[self.v[i]]
                self.z[i] = 0

    def suiji(self, mup, muq, sigmap, sigmaq, suijih, syst=-1.0, Pl=[], Ql=[], type='1'):

        # print('负荷随机波动模块启动:%s' % time.time())
        k = 0
        t = -1.0
        if type == '1':
            while True:
                t1 = time.time()
                try:
                    # tnew = syst.get()
                    tnew = syst.value
                except EOFError:
                    break
                if round(tnew, 2) != round(t, 2):

                    k = k + 1
                    t = tnew
                    # print('负荷随机波动模拟模块')
                    if round(round(t, 2) % suijih, 1) == 0.0:  # 是否有更好的表达形式

                        t2 = time.time()
                        # print('负荷随机波动模拟耗时：%s' % (t2 - t1))
                        for i in range(self.n):
                            # ttot = 0.0
                            # t3 = time.time()
                            a = random.normalvariate(mup[i], sigmap)
                            b = random.normalvariate(muq[i], sigmaq)
                            Pl[i] = a
                            Ql[i] = b
                            # t4 = time.time()
                            # ttot = ttot + t4 - t3
                            # print('一个负荷波动模拟耗时：%s' % (t4-t3))
                        t3 = time.time()
                        # print('所有负荷波动模拟总耗时：%s' % (t3-t2))
                    # time.sleep(0.05)

                    t2 = time.time()
                    # print('负荷随机波动模拟耗时：%s' % (t2-t1))
                    # print('负荷随机波动模块：第%s次，仿真时间%ss' % (k, t))
        if type == '2':

            t = syst
            t1 = time.time()
            if round(round(t, 2) % suijih, 1) == 0.0:  # 是否有更好的表达形式

                t3 = time.time()
                for i in range(self.n):
                    a = random.normalvariate(mup[i], sigmap)
                    b = random.normalvariate(muq[i], sigmaq)
                    system.PQ.Pl[i] = a
                    system.PQ.Ql[i] = b
                    # Pl.append(a)
                    # Ql.append(b)
                # system.PQ.Pl = Pl
                # system.PQ.Ql = Ql
                t4 = time.time()
                # print('所有负荷波动模拟总耗时：%s' % (t4 - t3))
            t2 = time.time()
            # print('负荷随机波动模拟耗时：%s' % (t2-t1))
            # print('负荷随机波动模块：仿真时间%ss' % t)

