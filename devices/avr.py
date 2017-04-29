"""

"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div

class avr1(base_device):
    def __init__(self):
        base_device.__init__(self)
        self.u = mul(matrix(self.u), matrix(system.PV.u))
        self._data.update(
            {'bus': None, 'Type': 1, 'vrmax': 3.3, 'vrmin': -2.6, 'K0': 0, 'T1': 0, 'T2': 0, 'T3': 0,
             'T4': 0, 'Te': 0, 'Tr': 0, 'Ae': 0, 'Be': 0})
        self._type = 'Avr1'
        self._name = 'Avr1'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['vref']  # 后续得修改
        self._states = ['vm', 'vr1', 'vr2', 'vf']  # 后续得修改
        self._params.extend(['u', 'vrmax', 'vrmin', 'K0', 'T1', 'T2', 'T3', 'T4', 'Te', 'Tr', 'Ae', 'Be'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改

    def setx0(self):

        vg = system.DAE.y[self.v]
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
            if vr[i] < self.vrmax[i]:
                print('Warn: vr1小于最小值vrmin')



class avr2(base_device):
    def __init__(self):
        base_device.__init__(self)
        self.u = mul(matrix(self.u), matrix(system.PV.u))
        self._data.update(
            {'bus': None, 'Type': 2, 'vrmax': 5, 'vrmin': 0, 'Ka': 0, 'Ta': 0, 'Kf': 0, 'Tf': 0,
             'Ke': 0, 'Te': 0, 'Tr': 0, 'Ae': 0, 'Be': 0})
        self._type = 'Avr1'
        self._name = 'Avr1'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['vref']  # 后续得修改
        self._states = ['vm', 'vr1', 'vr2', 'vf']  # 后续得修改
        self._params.extend(['u', 'vrmax', 'vrmin', 'Ka', 'Ta', 'Kf', 'Tf', 'Ke', 'Te', 'Tr', 'Ae', 'Be'])
        self._voltages = ['V0']  # 后续得修改
        self._powers = ['Pg']  # 后续得修改

    def setx0(self):

        vg = system.DAE.y[self.v]
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
        Ce = mul(self.Ke, vf) - mul(mul(self.Ae, exp(mul(self.Be, vf))), vf)


        system.DAE.x[self.vm] = mul(self.u, vg)
        self.vref0 = div(Ce, self.Ka) + vg

        system.DAE.x[self.vr1] = Ce
        system.DAE.x[self.vr2] = div(mul(-self.Kf, vf), self.Tf)
        system.DAE.x[self.vf] = vf

        system.DAE.y[self.vref] = mul(self.u, self.vref0)

        for i in range(self.n):
            if Ce[i] > self.vrmax[i]:
                print('Warn: vr2超出最大值vrmax')
            if Ce[i] < self.vrmax[i]:
                print('Warn: vr2小于最小值vrmin')

class avr3():
    def __init__(self):
        return