"""

"""
from devices.base_device import base_device
import system
import numpy as np
from cvxopt.base import matrix, spmatrix, cos, sin, sparse, mul,exp

class pq(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pl': 1, 'bus': None, 'Ql': 1,  'Vmax': 1.1, 'Vmin': 0.95})
        self._type = 'PQ'
        self._name = 'PQ'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._params.extend(['Pl', 'Ql',  'Vmax', 'Vmin'])
        self._powers = ['Pl', 'Ql']

    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.Ql):
            dae.g[key] += value
        for key, value in zip(self.a, self.Pl):
            dae.g[key] += value

    def gcall(self):
        for i in range(system.PQ.n):
            system.DAE.g[self.a[i]] += self.Pl[i]
            system.DAE.g[self.v[i]] += self.Ql[i]
        # 判断是否PQ负荷转为恒阻抗模型

        a = []
        b = []
        for i in range(self.n):
            if system.DAE.y[self.v[i]] < self.Vmin[i]:

                print('pq达到vmin')
                system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]

                system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]

            if system.DAE.y[self.v[i]] > self.Vmax[i]:

                print('pq达到vmax')
                system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] * system.DAE.y[self.v[i]]**2 /self.Vmin[i] / self.Vmin[i]

                system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] * system.DAE.y[self.v[i]]**2 / self.Vmin[i] / self.Vmin[i]



    # def Gycall(self):
    #
    #
    #
    #
    #     system.DAE.y = matrix(system.DAE.y)
    #     U = exp(system.DAE.y[self.a] * 1j)
    #     print(U)

        # system.DAE.Gy = np.zeros((2*system.PQ.n,2*system.PQ.n))
        #
        # for i in range(system.PQ.n):
        #     for j in range(system.PQ.n):
        #
        #         system.DAE.Gy[i,j] = -system.DAE.y[self.v[i]]*system.DAE.y[self.v[j]]*(system.DAE.Y_G[i,j]*sin(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]])- system.DAE.Y_B[i,j]*cos(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]]))
        #
        #         system.DAE.Gy[i,j+system.PQ.n] = -system.DAE.y[self.v[i]]*system.DAE.y[self.v[j]]*(system.DAE.Y_G[i,j]*cos(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]])+ system.DAE.Y_B[i,j]*sin(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]]))
        #         system.DAE.Gy[i+system.PQ.n,j] = system.DAE.y[self.v[i]]*system.DAE.y[self.v[j]]*(system.DAE.Y_G[i,j]*cos(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]])+ system.DAE.Y_B[i,j]*sin(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]]))
        #         system.DAE.Gy[i + system.PQ.n,j+system.PQ.n] = -system.DAE.y[self.v[i]]*system.DAE.y[self.v[j]]*(system.DAE.Y_G[i,j]*sin(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]])- system.DAE.Y_B[i,j]*cos(system.DAE.y[self.a[i]]-system.DAE.y[self.a[j]]))










