"""

"""
"""

"""
import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, sin, cos, spmatrix
import numpy as np

class fault(base_device):

    def __init__(self):

        base_device.__init__(self)
        self._data.update({'bus': None, 'Sn': 0, 'Vn': 0, 'fn': 50, 'tf': 0, 'tc': 0, 'rf': 1, 'xf': 1})
        self.u =[]
        self.properties.update({'gcall': True, 'Gycall': True})




        self._type = 'Fault'
        self._name = 'Fault'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._params.extend(['rf', 'xf'])       #要转成矩阵的元素
        self._z = ['rf', 'xf']
        self.y = 0+0*1j
        self.dat = [0, 0]
        self.n = 0




    def intervention(self,t):
        #

        for item in range(self.n):
            if t == self.tf[item]:      #发生故障
                print('故障发生在t = %s s' % self.tf[item])

                #启用故障

                self.u[item] = 1
                system.DAE.factorize = True

                #存储故障前
                self._ang = matrix(system.DAE.y[system.Bus.a])
                self._vol = matrix(system.DAE.y[system.Bus.n:len(system.DAE.y)])

            if t == self.tc[item]:    #故障清除

                print('在t = %s s清除故障' % self.tc[item])

                # 禁用故障
                self.u[item] = 0
                system.DAE.factorize = True

                #恢复电压
                if system.Settings.resetangles:
                    system.DAE.y[system.Bus.a] = matrix(self._ang)
                    system.DAE.y[system.Bus.n:len(system.DAE.y)] = matrix(self._vol)



    def gettimes(self):

        t = []
        if self.n == 0:
            return t

        # u = np.unique(self.tf,self.tc)    #此处省略同时发生故障的情况，待完善
        # print(self.tf[0])
        tf = self.tf[0] - 0.000001
        tc = self.tc[0] - 1e-6

        a = [tf, self.tf[0]]
        # print(a)
        b = [tc, self.tc[0]]
        a.extend(b)
        # print(a)
        t = a
        return t
        # t = [u,u]
        # print(t)


    def setup(self):
        self.u = matrix(0, (self.n, 1))
        z = self.rf + self.xf*1j
        h = div(1, z)
        self.y = h
        # self.y = h.H.T
        self.dat = [self.y.real(), self.y.imag()]

    def gcall(self):

        if self.n == 0:
            return

        V = system.DAE.y[self.v]
        V2 = mul(self.u, V, V)
        V2 = matrix(V2)




        system.DAE.g = system.DAE.g + spmatrix(mul(self.dat[0],V2),self.a,[0],(system.DAE.ny,1)) -spmatrix(mul(self.dat[1],V2),self.v,[0],(system.DAE.ny,1))

    def Gycall(self):

        if self.n == 0:
            return

        V = mul(2*self.u, system.DAE.y[self.v])


        system.DAE.Gy = system.DAE.Gy +spmatrix(mul(self.dat[0],V),self.a,self.v,(system.DAE.ny,system.DAE.ny)) -spmatrix(mul(self.dat[1],V),self.v,self.v,(system.DAE.ny,system.DAE.ny))




    def istime(self,t):

        if self.n == 0:
            return


        u = False
        for i in range(self.n):
            if self.tf[i] == t or self.tc[i] == t:
                u = True

        return u