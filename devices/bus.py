from devices.base_device import base_device
import system
from cvxopt.base import matrix, sin, cos



class bus(base_device):

    def __init__(self):

        base_device.__init__(self)
        self._data.update({'bus': None, 'Vb': 10.5, 'V_0': 1, 'theta0': 0, 'Pl': 0, 'Ql': 0, 'Pg': 0, 'Qg': 0})
        self._type = 'Bus'
        self._name = 'Bus'
        self.n = 0
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._params.extend(['Vb', 'Pl', 'Ql', 'Pg', 'Qg'])



    def _bus_index(self,):

        for index in self._bus.keys():
            for item in self.__dict__[index]:
                self.__dict__[self._bus[index][0]].append(self.int[item])
                self.__dict__[self._bus[index][1]].append(self.int[item] + self.n)

    def yinit(self, dae):
        zeros = [0.0] * (2*self.n)
        dae.y = zeros[:]
        dae.g = zeros[:]

        for i in range(self.n):
            if self.V_0[i] < 0.5:
                print('Warning: Bus %i initial guess voltage magnitudes are too low.' % i)
            if self.V_0[i] > 1.5:
                print('Warning: Bus %i initial guess voltage magnitudes are too high.' % i)
        system.DAE.y = matrix(system.DAE.y)
        system.DAE.y[self.v] = self.V_0

        aref = min(abs(matrix(self.theta0)))
        for i in range(self.n):
            if self.theta0[i] - aref < -1.5708:
                print('Warning: Bus %i initial guess voltage phases are too low.' % i)
            if self.theta0[i] - aref > 1.5708:
                print('Warning: Bus %i initial guess voltage phases are too high.' % i)
        system.DAE.y[self.a] = self.theta0








