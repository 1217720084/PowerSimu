"""

"""
from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, cos, sin, sparse, mul
import cvxopt.blas
import math
import numpy as np


class line(base_device):
    def __init__(self):
        base_device.__init__(self)
        self._data.update({'fn': 50, 'f_bus': None, 'to_bus': None, 'l': 0,  'kT': 1, 'r': 2, 'x': 1.1, 'b': 1, 'tap_ratio': 1.1, 'theta': 0.1,
                           'Imax': 0, 'Pmax': 0, 'Smax': 0})
        self._type = 'Line'
        self._name = 'Line'
        self._bus = {'f_bus': ['af', 'vf'], 'to_bus': ['at', 'vt']}
        self._params.extend(['fn', 'kT', 'r', 'x', 'tap_ratio', 'theta', 'Imax', 'Pmax', 'Smax'])
        self.z = ['r', 'x']

    def _bus_index(self):

        idx = []

        for index in self._bus.keys():

            for item in self.__dict__[index]:
                if not item in system.Bus.int:
                    continue
                    # self.message('Bus index <%s> does not exist', data_tuple = item, level = self.ERROR)
                else:
                    idx = system.Bus.int[item]
                    self.__dict__[self._bus[index][0]].append(system.Bus.a[idx])
                    self.__dict__[self._bus[index][1]].append(system.Bus.a[idx]+system.Bus.n)

    def build_y(self):

        for i in range(len(self.tap_ratio)):
            if self.tap_ratio[i] == 0:
                self.tap_ratio[i] = 1
        print(self.tap_ratio)

        # chrg = np.mat(self.b) * 0.5
        # chrg = chrg.T
        #
        # print(chrg)
        chrg = np.array(np.zeros((len(self.b), 1), complex))
        for i in range(len(self.x)):
            chrg[i] = complex(0, self.b[i] * 0.5)

        #print(chrg)
        #zeros = [0] * len(self.x)
        #y = matrix(zeros, (len(self.x), 1, complex))
        y = np.array(np.zeros((len(self.x), 1), complex))
        for i in range(len(self.x)):
            y[i] = 1 / complex(self.r[i], self.x[i])
        #print(y)

        ts = np.array(np.zeros((len(self.theta), 1), complex))
        for i in range(len(self.theta)):
            ts[i] = complex(self.tap_ratio[i]*cos(self.theta[i]*math.pi/180), self.tap_ratio[i]*sin(self.theta[i]*math.pi/180))
        #print(ts)

        ts2 = ts * ts.conj()
        #print(ts2)

        y1 = -y / ts.conj()
        #print(y1)

        self.Y = spmatrix(y1, self.af, self.at, (system.Bus.n, system.Bus.n)) +\
            spmatrix(y1, self.at, self.af, (system.Bus.n, system.Bus.n)) +\
            spmatrix(y/ts2 + chrg, self.af, self.af, (system.Bus.n, system.Bus.n)) +\
            spmatrix(y + chrg, self.at, self.at, (system.Bus.n, system.Bus.n))
        system.DAE.Y = self.Y

        system.DAE.Y_G = self.Y.real()
        print(system.DAE.Y_G)
        system.DAE.Y_B = self.Y.imag()
        print(system.DAE.Y_B)


        print(self.Y)
        print(self.Y.V)

    def gcall(self):
        # V0 = np.array(np.ones((system.Bus.n, 1), complex))
        # for i in system.Bus.v:
        #     V0[i] = system.DAE.y[i]
        # Va = np.array(np.ones((system.Bus.n, 1), complex))
        # for i in system.DAE.a:
        #     Va[i] = complex(0, system.DAE.y[i])
        Vc = np.array(np.ones((system.Bus.n, 1), complex))
        for i in range(system.Bus.n):
            Vc[i] = complex(system.DAE.y[i + system.Bus.n] * cos(system.DAE.y[i]), system.DAE.y[i + system.Bus.n] * sin(system.DAE.y[i]))

        Vc = matrix(Vc)
        print(Vc)
        print(Vc.size)
        print(self.Y.size)
        print(type(Vc))
        self.Y = matrix(self.Y)
        print(type(self.Y))

        I = np.matmul(self.Y, Vc)
        print(I)
        I = I.conj() #?
        #print(I.ndim)

        S = Vc * I
        S = matrix(S)
        print(S)
        self.p = S.real()
        print(self.p)
        self.q = S.imag()
        print(self.q)
        for i in range(system.Bus.n):
             system.DAE.g[i] = self.p[i]
             system.DAE.g[i+system.Bus.n] = self.q[i]


