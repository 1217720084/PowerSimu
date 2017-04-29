"""

"""
from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, mul,sparse
from numpy import multiply
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

        system.DAE.y = matrix(system.DAE.y)
        V = system.DAE.y[self.v]
        V2 = mul(V, V)
        p = mul(matrix(self.g), V2)
        q = mul(matrix(self.b), V2)



        # #V = matrix(0, (system.Shunt.n, 1))
        # V = [1.0] * system.Shunt.n
        # print(V)
        # for i in range(system.Shunt.n):
        #     V[i] = system.DAE.y[self.v[i]]
        # print(V)
        # V = np.array(V)
        # V2 = V * V
        # print(V2)
        # J = [1] * system.Shunt.n
        # print(self.a)
        # print(J)
        # p = [0] * system.Shunt.n
        # q = [0] * system.Shunt.n

        # for i in range(system.Shunt.n):
        #      p[i] = self.g[i] * V2[i]
        #      q[i] = self.b[i] * V2[i]

        for key, value in zip(self.a, p):
            system.DAE.g[key] += value
        for key, value in zip(self.v, q):
            system.DAE.g[key] -= value
    def Gycall(self):
        system.DAE.y=matrix(system.DAE.y)
        V=2*system.DAE.y[self.v]        #shunt bus 的索引？
        conductance = self.g   #网络中每条母线的电导列向量
        conductance = matrix(conductance)
        susceptance = self.b    #电纳列向量
        susceptance = matrix(susceptance)

        m=len(system.DAE.y)        #代数变量的个数
        conducv=mul(conductance,V)

        suscepv=mul(susceptance,V)


        spconducv = spmatrix(conducv, self.a, self.v, (m, m), 'z')

        spcsuscepv = spmatrix(suscepv, self.v,self.v, (m, m), 'z')

        system.DAE.Gy += spmatrix(conducv, self.a, self.v, (m, m), 'z')- spmatrix(suscepv, self.v,self.v, (m, m), 'z')
        system.DAE.Gy = sparse(system.DAE.Gy)
        system.DAE.Gy = matrix(system.DAE.Gy.real())
