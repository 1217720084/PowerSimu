"""

"""
from devices.base_device import base_device
import system


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

    def gcall(self,dae):
        dae.g[self.a] = self.P1 + self.Pg
        dae.g[self.v] = self.Q1 + self.Qg
        i = 0
        while i < dae.n_bus:#system.Bus.n
            dae.g[self.a] -= dae.y[v] * dae.y[i + dae.n_bus] * (system.DAE.Y_G[self.a][i] * cos(dae.y[a] - dae.y[i]) + system.DAE.Y_B[self.a][i] * sin(dae.y[a] - dae.y[i]))
            dae.g[self.v] -= dae.y[v] * dae.y[i + dae.n_bus] * (system.DAE.Y_G[self.a][i] * cos(dae.y[a] - dae.y[i]) - system.DAE.Y_B[self.a][i] * sin(dae.y[a] - dae.y[i]))
            i += 1