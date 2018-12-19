"""

"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, sin, cos, spmatrix
import numpy as np

class tg1(base_device):
    def __init__(self):

        base_device.__init__(self)

        # self.u = mul(matrix(self.u), matrix(system.Syn6.u))
        self._data.update({'bus': None, 'Type': 1, 'wref0': 1, 'R': 0.04, 'Pmax': 1.1, 'Pmin': 0, 'Ts': 0.03,
                           'Tc': 0.5, 'T3': 0, 'T4': 0.02, 'T5': 1})
        self._name = 'Tg1'
        self._type = 'Tg1'
        self.n = 0
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['wref']
        self._states = ['tg1', 'tg2', 'tg3']
        self._params.extend(['wref0', 'R', 'Pmax', 'Pmin', 'Ts', 'u',
                            'Tc', 'T3', 'T4', 'T5'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改
        self.ba = []
        self.bv = []
        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def getbus(self):
        for i in range(self.n):
            self.ba.append(system.Syn6.a[self.a[i]])
            self.bv.append(system.Syn6.v[self.a[i]])
    def _dxy_index(self):

        base_device._dxy_index(self)
        system.Syn6.pm = matrix(system.Syn6.pm)
        self.pm = system.Syn6.pm[self.a]   # 需要改进
        self.pm = list(self.pm)

    def base(self):

        self.R = div(mul(self.R, system.Settings.mva), system.Syn6.Sn[self.a])
        self.Pmax = div(mul(self.Pmax, system.Syn6.Sn[self.a]), system.Settings.mva)
        self.Pmin = div(mul(self.Pmin, system.Syn6.Sn[self.a]), system.Settings.mva)

    def setx0(self):

        if self.n == 0:
            return
        system.Syn6._list2matrix()
        self.Porder = system.Syn6.pm0[self.a]  # 需要改进，若不是Syn6
        gain = div(matrix(1.0, (self.n, 1)), self.R)

        a_s = div(matrix(1.0, (self.n, 1)), self.Ts)
        ac = div(matrix(1.0, (self.n, 1)), self.Tc)
        a5 = div(matrix(1.0, (self.n, 1)), self.T5)
        K1 = mul(self.T3, ac)
        K2 = 1 - K1
        K3 = mul(self.T4, a5)
        K4 = 1 - K3

        system.DAE.x[self.tg1] = mul(self.u, self.Porder)
        system.DAE.x[self.tg2] = mul(self.u, K2, self.Porder)
        system.DAE.x[self.tg3] = mul(self.u, K4, self.Porder)
        system.DAE.f[self.tg1] = 0
        system.DAE.f[self.tg2] = 0
        system.DAE.f[self.tg3] = 0

        for i in range(self.n):
            if self.u[i] == 1:
                system.Syn6.pm0[self.a[i]] = 0
        # system.Syn6.pm0[self.a] = 0
        system.DAE.y[self.wref] = mul(self.u, self.wref0)


        for i in range(self.n):
            if self.Porder[i] > self.Pmax[i]:
                print('第%i 台Tg1机械功率超过最大值'% (i+1))
            elif self.Porder[i] < self.Pmin[i]:
                print('第%i 台Tg1机械功率低于最小值' % (i+1))
            else:
                print('初始化Tg1完成')

        # self.pmech = system.DAE.x[self.tg3] + mul(K3, system.DAE.x[self.tg2]) + mul(K1, system.DAE.x[self.tg1])


    def gcall(self):

        if self.n == 0:
            return

        a_s = div(matrix(1.0, (self.n, 1)), self.Ts)
        ac = div(matrix(1.0, (self.n, 1)), self.Tc)
        a5 = div(matrix(1.0, (self.n, 1)), self.T5)
        K1 = mul(self.T3, ac)
        K2 = 1 - K1
        K3 = mul(self.T4, a5)
        K4 = 1 - K3
        pmech = system.DAE.x[self.tg3] + mul(K3, system.DAE.x[self.tg2]) + mul(K1, system.DAE.x[self.tg1])
        system.DAE.g = system.DAE.g \
                       + spmatrix(mul(self.u, pmech), self.pm, matrix(0, (self.n, 1)), (system.DAE.ny, 1)) \
                       + spmatrix(mul(self.u, self.wref0) - system.DAE.y[self.wref], self.wref, matrix(0, (self.n, 1)), (system.DAE.ny, 1))

    def Gycall(self):

        if self.n == 0:
            return
        system.DAE.Gy = system.DAE.Gy \
                        - spmatrix(matrix(1, (self.n, 1)), self.wref, self.wref, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        tg1 = system.DAE.x[self.tg1]
        tg2 = system.DAE.x[self.tg2]
        tg3 = system.DAE.x[self.tg3]
        wref = system.DAE.y[self.wref]


        gain = div(matrix(1.0, (self.n, 1)), self.R)

        a_s = div(matrix(1.0, (self.n, 1)), self.Ts)
        ac = div(matrix(1.0, (self.n, 1)), self.Tc)
        a5 = div(matrix(1.0, (self.n, 1)), self.T5)
        K1 = mul(self.T3, ac)
        K2 = 1 - K1
        K3 = mul(self.T4, a5)
        K4 = 1 - K3
        system.Syn6.omega = matrix(system.Syn6.omega)
        omega = system.Syn6.omega[self.a]

        tin = self.Porder + mul(gain, wref - system.DAE.x[omega])
        for i in range(self.n):
            tin[i] = max(tin[i], self.Pmin[i])
            tin[i] = min(tin[i], self.Pmax[i])

        system.DAE.f[self.tg1] = mul(self.u, a_s, -tg1+tin)
        system.DAE.f[self.tg2] = mul(self.u, ac, -tg2 + mul(K2, tg1))
        system.DAE.f[self.tg3] = mul(self.u, a5, -tg3 + mul(K4, tg2+mul(K1, tg1)))

    def Fxcall(self):

        if self.n == 0:
            return


        wref = system.DAE.y[self.wref]



        gain = div(matrix(1.0, (self.n, 1)), self.R)

        a_s = div(matrix(1.0, (self.n, 1)), self.Ts)
        ac = div(matrix(1.0, (self.n, 1)), self.Tc)
        a5 = div(matrix(1.0, (self.n, 1)), self.T5)
        K1 = mul(self.T3, ac)
        K2 = 1 - K1
        K3 = mul(self.T4, a5)
        K4 = 1 - K3
        system.Syn6.omega = matrix(system.Syn6.omega)
        omega = system.Syn6.omega[self.a]

        tin = self.Porder + mul(gain, wref - system.DAE.x[omega])
        for i in range(self.n):
            tin[i] = max(tin[i], self.Pmin[i])
            tin[i] = min(tin[i], self.Pmax[i])

        u = matrix(1, (self.n, 1))

        # windup limit
        for i in range(self.n):
            if (tin[i] < self.Pmax[i]) and (tin[i] > self.Pmin[i]):
                u[i] = 1
            else:
                u[i] = 0

        system.DAE.Fx = system.DAE.Fx \
                        - spmatrix(mul(u, self.u, a_s, gain), self.tg1, omega, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(a_s, self.tg1, self.tg1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(ac, self.tg2, self.tg2, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(a5, self.tg3, self.tg3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(ac, self.u, K2), self.tg2, self.tg1, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(a5, self.u, K4), self.tg3, self.tg2, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(a5, self.u, K1, K4), self.tg3, self.tg1, (system.DAE.nx, system.DAE.nx))

        system.DAE.Fy = system.DAE.Fy \
                        + spmatrix(mul(u, self.u, a_s, gain), self.tg1, self.wref, (system.DAE.nx, system.DAE.ny))

        system.DAE.Gx = system.DAE.Gx \
                        + spmatrix(self.u, self.pm, self.tg3, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(self.u, K3), self.pm, self.tg2, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(self.u, K1, K3), self.pm, self.tg1, (system.DAE.ny, system.DAE.nx))


class tg2(base_device):
    def __init__(self):
        base_device.__init__(self)
        return