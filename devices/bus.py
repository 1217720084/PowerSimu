from devices.base_device import base_device
import system
from cvxopt.base import matrix, sin, cos
import math
class bus(base_device):

    def __init__(self):

        base_device.__init__(self)
        self._data.update({'bus': None, 'Vb': 10.5, 'Pl': 0, 'Ql': 0, 'Pg': 0, 'Qg': 0})
        self._type = 'Bus'
        self._name = 'Bus'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._params.extend(['Vb', 'Pl', 'Ql', 'Pg', 'Qg'])

    def _bus_index(self,):

        for index in self._bus.keys():
            for item in self.__dict__[index]:
                self.__dict__[self._bus[index][0]].append(self.int[item])
                self.__dict__[self._bus[index][1]].append(self.int[item] + self.n)

    def yinit(self, dae):
        zeros = [0] * (2*self.n + 2*system.Statcom.n + 2*system.Sssc.n)
        dae.y = zeros[:]
        dae.g = zeros[:]
        #dae.Gy = sparse(m*m)
        for item in self.Va:
            dae.y[item] = 0.0
        for item in self.V0:
            dae.y[item] = 1.0
        for item in range(system.Statcom.n):
            dae.y[2*self.n+2*item] = 0.0
            dae.y[2*self.n+2*item+1] = 1.0
        for item in range(system.Sssc.n):
            dae.y[2*self.n+2*system.Statcom.n+2*item] = -math.pi/2
            dae.y[2*self.n+2*system.Statcom.n+2*item+1] = 0.1







