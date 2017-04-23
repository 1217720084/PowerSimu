"""

"""
from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix
import numpy as np


class shunt(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'fn': 50, 'bus': None, 'g': 1, 'b': 1.1})
        self._type = 'Shunt'
        self._name = 'Shunt'
        self._bus = {'bus': ['a', 'v']}
        self._params.extend(['fn', 'g', 'b'])
        self._y = ['g', 'b']



    def gcall(self):

        #V = matrix(0, (system.Shunt.n, 1))
        V = [1.0] * system.Shunt.n
        print(V)
        for i in range(system.Shunt.n):
            V[i] = system.DAE.y[self.v[i]]
        print(V)
        V = np.array(V)
        V2 = V * V
        print(V2)
        J = [1] * system.Shunt.n
        print(self.a)
        print(J)
        p = [0] * system.Shunt.n
        q = [0] * system.Shunt.n

        for i in range(system.Shunt.n):
             p[i] = self.g[i] * V2[i]
             q[i] = self.b[i] * V2[i]

        for key, value in zip(self.a, p):
            system.DAE.g[key] += value
        for key, value in zip(self.v, q):
            system.DAE.g[key] += value
