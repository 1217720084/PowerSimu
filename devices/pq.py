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
        self._algebs = ['Va', 'V0']
        self._params.extend(['Pl', 'Ql',  'Vmax', 'Vmin'])
        self._powers = ['Pl', 'Ql']

    def yinit(self, a1):

        a1.y = system.DAE.y
        a1.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.Ql):
            a1.g[key] += value
        for key, value in zip(self.a, self.Pl):
            a1.g[key] += value