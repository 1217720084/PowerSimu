"""

"""
from devices.base_device import base_device
import system


class pv(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pg': 1, 'bus': None, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05, 'Vmax': 1.1, 'Vmin': 0.95,})
        self._type = 'PV'
        self._name = 'PV'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['Va', 'V0']
        self._params.extend(['Pg', 'qgmax', 'qgmin', 'V0', 'Vmax', 'Vmin',])
        self.n_PV = 0
        self._voltages = ['V0']
        self._powers = ['Pg']


    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        #dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.V0):
            dae.y[key] = value
        for key, value in zip(self.a, self.Pg):
            dae.g[key] += value


class slack(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pg': 1, 'bus': None, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05, 'Vmax': 1.1, 'Vmin': 0.95, 'Va': 0})
        self._type = 'SW'
        self._name = 'SW'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._params.extend(['Pg', 'qgmax', 'qgmin', 'V0', 'Vmax', 'Vmin', 'Va'])
        self._voltages = ['V0']
        self._powers = ['Pg']

    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.V0):
            dae.y[key] = value
        for key, value in zip(self.a, self.Va):
            dae.y[key] = value
        for key, value in zip(self.a, self.Pg):
            dae.g[key] = 0
        for key, value in zip(self.v, self.Pg):
            dae.g[key] = 0

    def gcall(self,dae):
        dae.g[self.v] = self.P1 + self.Pg
        i = 0
        while i < self.n: #PV,PQ节点索引号排列顺序？
            dae.g[self.v] -= dae.y[v] * dae.y[i] * (G * cos(dae.y[a] - dae.y[i]) + B * sin(dae.y[a] - dae.y[i]))
            i += 1
        return