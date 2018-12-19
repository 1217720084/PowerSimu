"""

"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, sin, cos, spmatrix
import copy
import numpy as np

class pss1(base_device):
    def __init__(self):
        base_device.__init__(self)
        return


class pss2(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'avr': None, 'm_model': 2, 'ni': 1, 'vsmax': 0.1, 'vsmin': -0.1, 'Kw': 8, 'Tw': 10, 'T1': 0.8,
                           'T2': 0.05, 'T3': 0.8, 'T4': 0.05, 'Ka': 0, 'Ta': 0, 'Kp': 0, 'Kv': 0, 'vamax': 0, 'va1max': 0,
                           'vs1max': 0, 'vs1min': 0, 'ethr': 0, 'wthr': 0, 's2': 0, 'u': 1, 'syn': 0})
        self._type = 'Pss2'
        self._name = 'Pss2'
        # self._bus = {'bus': ['a', 'v']}
        self.n = 0
        self._algebs = ['vss']
        self._states = ['v1', 'v2', 'v3']
        self._params.extend(['ni', 'vsmax', 'vsmin', 'Kw', 'Tw', 'T1', 'T2', 'T3',
                             'T4', 'Ka', 'Ta', 'Kp', 'Kv', 'vamax', 'va1max',
                           'vs1max', 'vs1min', 'ethr', 'wthr', 's2'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改
        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def _bus_index(self):


        # 后续得修改
        system.Avr2.a = matrix(system.Avr2.a)
        for i in range(self.n):
            self.syn[i] = system.Avr2.a[self.avr[i]]
        system.Syn6.a = matrix(system.Syn6.a)
        system.Syn6.v = matrix(system.Syn6.v)
        self.a = [0]*self.n
        self.v = [0]*self.n
        self.a = system.Syn6.a[self.syn]
        self.v = system.Syn6.v[self.syn]
        # the PSS is inactive if the AVR is off-line
        for i in range(self.n):
            self.u[i] = self.u[i]*system.Avr2.u[self.avr[i]]

    def _dxy_index(self):

        base_device._dxy_index(self)

        self.s1 = [0]*self.n
        self.omega = [0] * self.n
        self.p = [0] * self.n
        self.vf = [0] * self.n
        self.vref = [0] * self.n
        system.Syn6.omega = matrix(system.Syn6.omega)
        system.Syn6.p = matrix(system.Syn6.p)
        system.Syn6.vf = matrix(system.Syn6.vf)
        system.Avr2.vref = matrix(system.Avr2.vref)
        self.omega = system.Syn6.omega[self.syn]
        self.p = system.Syn6.p[self.syn]
        self.vf = system.Syn6.vf[self.syn]
        self.vref = system.Avr2.vref[self.avr]
        print('1')




    def setx0(self):

        if self.n == 0:
            return

        VSI = [0]*self.n
        SIw = [0]*self.n
        SIp = [0] * self.n
        SIv = [0] * self.n
        for i in range(self.n):
            if self.ni[i] == 1:
                SIw[i] = i
                VSI[i] = system.DAE.x[self.omega[SIw[i]]]
            elif self.ni[i] == 2:
                SIp[i] = i
                VSI[i] = system.DAE.y[self.p[SIp[i]]]
                self.Kp[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kv[i] = 0.0
            elif self.ni[i] == 3:
                SIv[i] = i
                VSI[i] = system.DAE.y[self.v[SIv[i]]]
                self.Kv[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kp[i] = 0.0
        VSI = matrix(VSI)
        # VSI = system.DAE.x[self.omega[SIw]]



        Kw = mul(self.u, self.Kw)
        Kp = mul(self.u, self.Kp)
        Kv = mul(self.u, self.Kv)

        Tw = self.Tw
        T2 = self.T2
        T4 = self.T4
        Ta = self.Ta

        for i in range(self.n):
            if Tw[i] == 0.0:
                Tw[i] = 0.01
                print(' Tw cannot be zero. Default value Tw = 0.01 will be used')

            if T2[i] == 0.0:
                T2[i] = 0.01
                print(' T2 cannot be zero. Default value Tw = 0.01 will be used')
            if T4[i] == 0.0:
                T4[i] = 0.01
                print(' T2 cannot be zero. Default value Tw = 0.01 will be used')

        system.DAE.x[self.v1] = mul(-(Kw+Kp+Kv), VSI)
        system.DAE.x[self.v2] = 0
        system.DAE.x[self.v3] = 0

        for i in range(self.n):
            if self.s2[i] == 0.0:
                self.s2[i] = -1.0
        system.DAE.y[self.vss] = 0.0

    def gcall(self):

        if self.n == 0:
            return

        VSI = [0] * self.n
        SIw = [0] * self.n
        SIp = [0] * self.n
        SIv = [0] * self.n
        for i in range(self.n):
            if self.ni[i] == 1:
                SIw[i] = i
                VSI[i] = system.DAE.x[self.omega[SIw[i]]]
            elif self.ni[i] == 2:
                SIp[i] = i
                VSI[i] = system.DAE.y[self.p[SIp[i]]]
                self.Kp[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kv[i] = 0.0
            elif self.ni[i] == 3:
                SIv[i] = i
                VSI[i] = system.DAE.y[self.v[SIv[i]]]
                self.Kv[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kp[i] = 0.0
        VSI = matrix(VSI)


        T1 = self.T1
        T2 = self.T2
        T3 = self.T3
        T4 = self.T4

        Kw = self.Kw
        Kp = self.Kp
        Kv = self.Kv

        vsmax = self.vsmax
        vsmin = self.vsmin
        vathr = self.va1max
        v3max = self.vs1max
        v3min = self.vs1min
        S2 = copy.deepcopy(self.s2)

        for i in range(self.n):
            if system.DAE.x[self.omega[i]]-1 < 0 or S2[i]:
                if S2[i] >= 0:
                    S2[i] = 1
                else:
                    S2[i] = 0
            else:
                S2[i] = 0

        Vs = [0]*self.n
        y = mul(Kw+Kp+Kv, VSI) + system.DAE.x[self.v1]

        A = div(T1, T2)
        C = div(T3, T4)
        Vs = system.DAE.x[self.v3] + mul(C, system.DAE.x[self.v2] + mul(A, y))

        for i in range(self.n):
            if Vs[i] >= vsmin[i]:
                Vs[i] = Vs[i]
            else:
                Vs[i] = vsmin[i]

            if Vs[i] >= vsmax[i]:
                Vs[i] = vsmax[i]
            else:
                Vs[i] = Vs[i]

        zero = [0]*self.n
        system.DAE.g = system.DAE.g + spmatrix(mul(self.u, Vs-system.DAE.y[self.vss]), self.vss, zero, (system.DAE.ny, 1))\
                            + spmatrix(mul(self.u, system.DAE.y[self.vss]), self.vref, zero, (system.DAE.ny, 1))


    def Gycall(self):

        if self.n == 0:
            return

        T1 = self.T1
        T2 = self.T2
        T3 = self.T3
        T4 = self.T4

        Kp = self.Kp
        Kv = self.Kv

        vsmax = self.vsmax
        vsmin = self.vsmin

        vss = system.DAE.y[self.vss]
        z = [0]*self.n
        z = matrix(z)
        for i in range(self.n):
            if vss[i] < vsmax[i] and vss[i] > vsmin[i]:
                if self.u[i]:
                    z[i] = 1

        system.DAE.Gy = system.DAE.Gy - spmatrix(1, self.vss, self.vss, (system.DAE.ny, system.DAE.ny))\
                             + spmatrix(z, self.vref, self.vss, (system.DAE.ny, system.DAE.ny))

        A = div(T1, T2)
        C = div(T3, T4)
        E = mul(C, A)
        system.DAE.Gy = system.DAE.Gy + spmatrix(mul(z, mul(Kv, E)), self.vss, self.v, (system.DAE.ny, system.DAE.ny)) \
                        + spmatrix(mul(z, mul(Kp, E)), self.vss, self.p, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        VSI = [0] * self.n
        SIw = [0] * self.n
        SIp = [0] * self.n
        SIv = [0] * self.n
        for i in range(self.n):
            if self.ni[i] == 1:
                SIw[i] = i
                VSI[i] = system.DAE.x[self.omega[SIw[i]]]
            elif self.ni[i] == 2:
                SIp[i] = i
                VSI[i] = system.DAE.y[self.p[SIp[i]]]
                self.Kp[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kv[i] = 0.0
            elif self.ni[i] == 3:
                SIv[i] = i
                VSI[i] = system.DAE.y[self.v[SIv[i]]]
                self.Kv[i] = self.Kw[i]
                self.Kw[i] = 0.0
                self.Kp[i] = 0.0
        VSI = matrix(VSI)

        Tw = self.Tw
        T1 = self.T1
        T2 = self.T2
        T3 = self.T3
        T4 = self.T4
        Ta = self.Ta


        Kw = self.Kw
        Ka = self.Ka
        Kp = self.Kp
        Kv = self.Kv

        ETHR = self.ethr
        WTHR = self.wthr
        S2 = copy.deepcopy(self.s2)

        for i in range(self.n):
            if system.DAE.x[self.omega[i]] - 1 < 0 or S2[i]:
                if S2[i] >= 0:
                    S2[i] = 1
                else:
                    S2[i] = 0
            else:
                S2[i] = 0

        y = mul(Kw+Kp+Kv, VSI) + system.DAE.x[self.v1]
        system.DAE.f[self.v1] = -mul(self.u, div(y, Tw))

        A = div(T1, T2)
        B = 1 - A
        C = div(T3, T4)
        D = 1 - C

        system.DAE.f[self.v2] = div(mul(self.u, mul(B, y)-system.DAE.x[self.v2]), T2)
        system.DAE.f[self.v3] = div(mul(self.u, mul(D, system.DAE.x[self.v2]+mul(A, y))-system.DAE.x[self.v3]), T4)




    def Fxcall(self):

        if self.n == 0:
            return

        Tw = self.Tw
        T1 = self.T1
        T2 = self.T2
        T3 = self.T3
        T4 = self.T4
        Ta = self.Ta

        Kw = self.Kw
        Ka = self.Ka
        Kp = self.Kp
        Kv = self.Kv

        vsmax = self.vsmax
        vsmin = self.vsmin
        vamax = self.vamax
        vathr = self.va1max
        S2 = copy.deepcopy(self.s2)

        for i in range(self.n):
            if system.DAE.x[self.omega[i]] - 1 < 0 or S2[i]:
                if S2[i] >= 0:
                    S2[i] = 1
                else:
                    S2[i] = 0
            else:
                S2[i] = 0

        vss = system.DAE.y[self.vss]
        z = [0] * self.n
        z = matrix(z)
        for i in range(self.n):
            if vss[i] < vsmax[i] and vss[i] > vsmin[i]:
                if self.u[i]:
                    z[i] = 1

        # common Jacobians elements

        system.DAE.Fx = system.DAE.Fx \
                        - spmatrix(div(1, Tw), self.v1, self.v1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(mul(self.u, div(Kw, Tw)), self.v1, self.omega, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy \
                        - spmatrix(mul(self.u, div(Kv, Tw)), self.v1, self.v, (system.DAE.nx, system.DAE.ny)) \
                        - spmatrix(mul(self.u, div(Kp, Tw)), self.v1, self.p, (system.DAE.nx, system.DAE.ny))

        A = div(T1, T2)
        B = 1 - A
        C = div(T3, T4)
        D = 1 - C
        E = mul(C, A)
        F = mul(self.u, div(D, T4))
        G = mul(self.u, div(B, T2))

        system.DAE.Fx = system.DAE.Fx \
                        + spmatrix(G, self.v2, self.v1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(1, T2), self.v2, self.v2, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(G, Kw), self.v2, self.omega, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(F, A), self.v3, self.v1, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(F, self.v3, self.v2, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(1, T4), self.v3, self.v3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(F, mul(A, Kw)), self.v3, self.omega, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy \
                        + spmatrix(mul(G, Kv), self.v2, self.v, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(mul(G, Kp), self.v2, self.p, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(mul(F, mul(A, Kv)), self.v3, self.v, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(mul(F, mul(A, Kp)), self.v3, self.p, (system.DAE.nx, system.DAE.ny))
        system.DAE.Gx = system.DAE.Gx \
                        + spmatrix(z, self.vss, self.v3, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(z, C), self.vss, self.v2, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(z, E), self.vss, self.v1, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(z, mul(Kw, E)), self.vss, self.omega, (system.DAE.ny, system.DAE.nx))

    def windup(self, type):

        if type == 'td':

            return

    def mosig(self, busangle, avr, ni=2):

        # ni 指pss输入信号，1为omega,2为p,3为v
        # syn, avr为列表形式，如：[1,2]，定义pss所修改的信号
        self.ni = ni  # 整数
        self.avr = avr
        system.Avr2.vref = matrix(system.Avr2.vref)
        system.Bus.a = matrix(system.Bus.a)
        if self.ni == 1:
            self.omega = system.Bus.a[busangle]
        elif self.ni == 2:
            self.p = system.Bus.a[busangle]
        elif self.ni == 3:
            self.v = system.Bus.a[busangle]

        self.vref = system.Avr2.vref[self.avr]






class pss3(base_device):
    def __init__(self):
        base_device.__init__(self)
        return