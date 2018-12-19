"""

"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, spmatrix
import random

class avr1(base_device):
    def __init__(self):
        base_device.__init__(self)
        self.u = mul(matrix(self.u), matrix(system.PV.u))
        self._data.update(
            {'bus': None, 'Type': 1, 'vrmax': 3.3, 'vrmin': -2.6, 'K0': 0, 'T1': 0, 'T2': 0, 'T3': 0,
             'T4': 0, 'Te': 0, 'Tr': 0, 'Ae': 0, 'Be': 0})
        self._type = 'Avr1'
        self._name = 'Avr1'
        self.n = 0
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['vref']  # 后续得修改
        self._states = ['vm', 'vr1', 'vr2', 'vf']  # 后续得修改
        self._params.extend(['u', 'vrmax', 'vrmin', 'K0', 'T1', 'T2', 'T3', 'T4', 'Te', 'Tr', 'Ae', 'Be'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改
        self.ba = []
        self.bv = []

        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def setx0(self):

        if self.n == 0:
            return

        vg = system.DAE.y[self.bv]
        vf = mul(self.u, system.Syn6.vf0[self.a])
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)

        # 检查参数

        for i in range(self.n):
            if self.Te[i] <= 0:
                self.Te[i] = 1
                print('<%i> Te不能小于等于0，设置Te = 1 [s]' % i)
            if self.Tr[i] <= 0:
                self.Tr[i] = 0.001
                print('<%i> Tr不能小于等于0，设置Te = 0.001 [s]' % i)
            if self.K0[i] <= 0:
                self.K0[i] = 400
                print('<%i> K0不能小于等于0，设置K0 = 400 ' % i)
            if self.T1[i] <= 0:
                self.T1[i] = 0.1
                print('<%i> T1不能小于等于0，设置T1 = 0.1 [s]' % i)
            if self.T4[i] <= 0:
                self.T4[i] = 0.01
                print('<%i> T4不能小于等于0，设置T4 = 0.01 [s]' % i)

        #
        Ce = -vf - mul(mul(self.Ae, exp(mul(self.Be, vf))), vf)
        K1 = mul(self.K0, div(self.T2, self.T1))
        K2 = self.K0 - K1
        K3 = div(self.T3, self.T4)
        K4 = 1 - K3

        system.DAE.x[self.vm] = mul(self.u, vg)
        self.vref0 = vg - div(Ce, self.K0)

        system.DAE.x[self.vr1] = mul(self.u, mul(K2, self.vref0-vg))
        system.DAE.x[self.vr2] = mul(self.u, mul(K4, self.vref0 - vg))
        system.DAE.x[self.vf] = vf

        system.DAE.y[self.vref] = mul(self.u, self.vref0)

        vr = mul(self.u, mul(self.K0, system.DAE.x[self.vr2])) + mul(K3, mul(K1, self.vref0-vg)+system.DAE.x[self.vr1])
        for i in range(self.n):
            if vr[i] > self.vrmax[i]:
                print('Warn: vr1超出最大值vrmax')
            if vr[i] < self.vrmin[i]:
                print('Warn: vr1小于最小值vrmin')

        system.Syn6.vf0[self.a] = 0
    def getbus(self):
        for i in range(self.n):
            self.ba.append(system.Syn6.a[self.a[i]])
            self.bv.append(system.Syn6.v[self.a[i]])

    def gcall(self):

        if self.n == 0:
            return

        self.vfd = system.Syn6.vf[self.a]



        system.DAE.g = system.DAE.g + spmatrix(system.DAE.x[self.vf], self.vfd, [0]*self.n, (system.DAE.ny, 1)) \
                       + spmatrix(mul(self.u, self.vref0)-system.DAE.y[self.vref], self.vref, [0]*self.n, (system.DAE.ny, 1))


    def Gycall(self):

        if self.n == 0:
            return

        system.DAE.Gy = system.DAE.Gy - spmatrix([1]*self.n, self.vref, self.vref, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        vg = system.DAE.y[self.bv]
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)
        vm = system.DAE.x[self.vm]
        vr1 = system.DAE.x[self.vr1]
        vr2 = system.DAE.x[self.vr2]
        vf = system.DAE.x[self.vf]
        vref = system.DAE.y[self.vref]

        system.DAE.f[self.vm] = div(mul(self.u, vg-system.DAE.x[self.vm]), self.Tr)

        K1 = div(mul(self.K0, self.T2), self.T1)
        K2 = self.K0 - K1
        K3 = div(self.T4, self.T3)
        K4 = matrix(1, (self.n, 1)) - K3

        vr = mul(self.K0, vr2) + mul(K3, mul(K1, vref-vm)+vr1)

        system.DAE.f[self.vr1] = div(mul(self.u, mul(K2, vref-vm)-vr1), self.T1)
        system.DAE.f[self.vr2] = div(mul(self.u, mul(K4, vr1+mul(K1, vref - vm))-mul(self.K0, vr2)), mul(self.T3, self.K0))

        # hard limit

        for i in range(self.n):
            vr[i] = min(vr[i], vrmax[i])
            vr[i] = max(vr[i], vrmin[i])



        Se = mul(self.Ae, exp(mul(self.Be, abs(vf))))

        system.DAE.f[self.vf] = div(mul(self.u, -vf+vr-mul(Se, vf)), self.Te)

    def Fxcall(self):

        if self.n == 0:
            return

        vg = system.DAE.y[self.bv]
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)

        system.DAE.Gx = system.DAE.Gx + spmatrix(self.u, self.vfd, self.vf, (system.DAE.ny, system.DAE.nx))
        system.DAE.Fx = system.DAE.Fx - spmatrix(div(matrix(1.0, (self.n, 1)), self.Tr), self.vm, self.vm, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy + spmatrix(div(self.u, self.Tr), self.vm, self.bv, (system.DAE.nx, system.DAE.ny))

        vm = system.DAE.x[self.vm]
        vr1 = system.DAE.x[self.vr1]
        vr2 = system.DAE.x[self.vr2]
        vf = system.DAE.x[self.vf]
        vref = system.DAE.y[self.vref]

        K1 = div(mul(self.K0, self.T2), self.T1)
        K2 = self.K0 - K1
        K3 = div(self.T4, self.T3)
        K4 = matrix(1, (self.n, 1)) - K3

        vr = mul(self.K0, vr2) + mul(K3, mul(K1, vref - vm) + vr1)

        Se = mul(self.Ae, exp(mul(self.Be, abs(vf))))

        Se = Se + mul(mul(Se, self.Ae), vf)

        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if vr[i] < vrmax[i] and vr[i] > vrmin[i]:
                z[i] = 1

        system.DAE.Fx = system.DAE.Fx \
                        - spmatrix(div(K2, self.T1), self.vr1, self.vm, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(matrix(1.0, (self.n, 1)), self.T1), self.vr1, self.vr1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(mul(K4, K1), mul(self.K0, self.T3)), self.vr2, self.vm, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(div(K4, mul(self.T3, self.K0)), self.vr2, self.vr1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(matrix(1.0, (self.n, 1)), self.T3), self.vr2, self.vr2, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(1+Se, self.Te), self.vf, self.vf, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(div(mul(z, K3), self.Te), self.vf, self.vr1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(mul(mul(z, K3), K1), self.Te), self.vf, self.vm, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(div(mul(mul(z, self.u), self.K0), self.Te), self.vf, self.vr2, (system.DAE.nx, system.DAE.nx))

        system.DAE.Fy = system.DAE.Fy \
                        + spmatrix(div(K2, self.T1), self.vr1, self.vref, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(div(mul(K4, K1), mul(self.K0, self.T3)), self.vr2, self.vref, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(div(mul(mul(z, K3), K1), self.Te), self.vf, self.vref, (system.DAE.nx, system.DAE.ny))

    def suiji(self, muvref01, sigmavref):

        if self.n == 0:
            return

        for i in range(self.n):
            a = random.normalvariate(muvref01[i], sigmavref)

            self.vref0[i] = a


class avr2(base_device):
    def __init__(self):
        base_device.__init__(self)
        self.u = mul(matrix(self.u), matrix(system.PV.u))
        self._data.update(
            {'bus': None, 'Type': 2, 'vrmax': 5, 'vrmin': 0, 'Ka': 0, 'Ta': 0, 'Kf': 0, 'Tf': 0,
             'Ke': 0, 'Te': 0, 'Tr': 0, 'Ae': 0, 'Be': 0})
        self._type = 'Avr2'
        self._name = 'Avr2'
        self.n = 0
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['vref']  # 后续得修改
        self._states = ['vm', 'vr1', 'vr2', 'vf']  # 后续得修改
        self._params.extend(['u', 'vrmax', 'vrmin', 'Ka', 'Ta', 'Kf', 'Tf', 'Ke', 'Te', 'Tr', 'Ae', 'Be'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改
        self.ba = []
        self.bv = []

        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})


    def setx0(self):

        vg = system.DAE.y[self.bv]
        vf = mul(self.u, system.Syn6.vf0[self.a])
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)

        # 检查参数

        for i in range(self.n):
            if self.Te[i] == 0:
                self.Te[i] = 1
                print('<%i> Te不能小于等于0，设置Te = 1 [s]' % i)
            if self.Tr[i] <= 0:
                self.Tr[i] = 0.001
                print('<%i> Tr不能小于等于0，设置Te = 0.001 [s]' % i)
            if self.Ke[i] <= 0:
                self.Ke[i] = 1
                print('<%i> Ke不能小于等于0，设置Ke = 1 ' % i)
            if self.Tf[i] <= 0:
                self.Tf[i] = 0.1
                print('<%i> T1不能小于等于0，设置T1 = 0.1 [s]' % i)
            if self.Ta[i] <= 0:
                self.Ta[i] = 0.1
                print('<%i> T4不能小于等于0，设置T4 = 0.1 [s]' % i)

        #
        Ce = mul(self.Ke, vf) + mul(mul(self.Ae, exp(mul(self.Be, abs(vf)))), vf)


        system.DAE.x[self.vm] = mul(self.u, vg)
        self.vref0 = div(Ce, self.Ka) + vg

        system.DAE.x[self.vr1] = Ce
        system.DAE.x[self.vr2] = div(mul(-self.Kf, vf), self.Tf)
        system.DAE.x[self.vf] = vf

        system.DAE.y[self.vref] = mul(self.u, self.vref0)

        for i in range(self.n):
            if Ce[i] > self.vrmax[i]:
                print('Warn: vr1超出最大值vrmax')
            if Ce[i] < self.vrmin[i]:
                print('Warn: vr1小于最小值vrmin')

        for i in range(self.n):
            if self.u[i] == 1:
                system.Syn6.vf0[self.a[i]] = 0
        # system.Syn6.vf0[self.a] = 0

    def getbus(self):
        for i in range(self.n):
            self.ba.append(system.Syn6.a[self.a[i]])
            self.bv.append(system.Syn6.v[self.a[i]])



    def gcall(self):

        self.vfd = system.Syn6.vf[self.a]

        system.DAE.g = system.DAE.g + spmatrix(system.DAE.x[self.vf], self.vfd, [0] * self.n, (system.DAE.ny, 1)) \
                       + spmatrix(mul(self.u, self.vref0) - system.DAE.y[self.vref], self.vref, [0] * self.n, (system.DAE.ny, 1))


    def Gycall(self):

        system.DAE.Gy = system.DAE.Gy - spmatrix([1]*self.n, self.vref, self.vref, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        vg = system.DAE.y[self.bv]
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)
        vm = system.DAE.x[self.vm]
        vr1 = system.DAE.x[self.vr1]
        vr2 = system.DAE.x[self.vr2]
        vf = system.DAE.x[self.vf]
        vref = system.DAE.y[self.vref]

        system.DAE.f[self.vm] = div(mul(self.u, vg - system.DAE.x[self.vm]), self.Tr)  #

        K5 = div(self.Kf, self.Tf)

        system.DAE.f[self.vr1] = div(mul(self.u, mul(self.Ka, vref - vm-vr2-mul(K5, vf)) - vr1), self.Ta)
        system.DAE.f[self.vr2] = div(mul(-self.u, mul(K5, vf) + vr2), self.Tf)

        # non-windup limit

        for i in range(self.n):
            if vr1[i] >= vrmax[i] and system.DAE.f[self.vr1[i]] > 0:
                system.DAE.f[self.vr1[i]] = 0
            if vr1[i] <= vrmin[i] and system.DAE.f[self.vr1[i]] < 0:
                system.DAE.f[self.vr1[i]] = 0

        for i in range(self.n):
            vr1[i] = min(vr1[i], vrmax[i])
            vr1[i] = max(vr1[i], vrmin[i])

        system.DAE.x[self.vr1] = mul(self.u, vr1)

        Se = mul(self.Ae, exp(mul(self.Be, abs(vf))))

        system.DAE.f[self.vf] = div(mul(self.u, vr1 - mul(self.Ke, vf) - mul(Se, vf)), self.Te)

    def Fxcall(self):

        vg = system.DAE.y[self.bv]
        vrmax = mul(self.u, self.vrmax)
        vrmin = mul(self.u, self.vrmin)

        system.DAE.Gx = system.DAE.Gx + spmatrix(self.u, self.vfd, self.vf, (system.DAE.ny, system.DAE.nx))
        system.DAE.Fx = system.DAE.Fx - spmatrix(div(matrix(1.0, (self.n, 1)), self.Tr), self.vm, self.vm, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy + spmatrix(div(self.u, self.Tr), self.vm, self.bv, (system.DAE.nx, system.DAE.ny))


        vr1 = system.DAE.x[self.vr1]
        vf = system.DAE.x[self.vf]
        Ka = mul(self.u, self.Ka)
        Kf = mul(self.u, self.Kf)
        K5 = div(Kf, self.Tf)

        Se = mul(self.Ae, exp(mul(self.Be, abs(vf))))

        Se = Se + mul(mul(Se, self.Be), vf)

        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if vr1[i] < vrmax[i] and vr1[i] > vrmin[i]:
                z[i] = 1

        system.DAE.Fx = system.DAE.Fx \
                        - spmatrix(div(matrix(1.0, (self.n, 1)), self.Tf), self.vr2, self.vr2, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(K5, self.Tf), self.vr2, self.vf, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(self.Ke+Se, self.Te), self.vf, self.vf, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(div(z, self.Te), self.vf, self.vr1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(mul(z, Ka), self.Ta), self.vr1, self.vm, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(matrix(1.0, (self.n, 1)), self.Ta), self.vr1, self.vr1, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(mul(z, Ka), self.Ta), self.vr1, self.vr2, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(mul(mul(z, K5), Ka), self.Ta), self.vr1, self.vf, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy \
                        + spmatrix(div(mul(z, Ka), self.Ta), self.vr1, self.vref, (system.DAE.nx, system.DAE.ny))

    def windup(self, type):
        idx = self.vr1
        xmax = self.vrmax
        xmin = self.vrmin

        if type == 'td':
            x = system.DAE.x[idx]
            for i in range(len(x)):
                if x[i] >= xmax[i] or x[i] <= xmin[i]:
                    if system.DAE.f[idx[i]] == 0:

                        system.DAE.tn[idx[i]] = 0
                        system.DAE.Ac[idx[i], :] = 0
                        system.DAE.Ac[:, idx[i]] = 0
                        system.DAE.Ac = system.DAE.Ac - spmatrix(1.0, [idx[i]], [idx[i]], ((system.DAE.nx+system.DAE.ny, system.DAE.nx+system.DAE.ny)))


    def suiji(self, muvref02, sigmavref):

        for i in range(self.n):
            a = random.normalvariate(muvref02[i], sigmavref)

            self.vref0[i] = a



class avr3(base_device):
    def __init__(self):
        base_device.__init__(self)
        return