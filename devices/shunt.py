"""

"""
from devices.base_device import base_device
import system


class shunt(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'fn': 50, 'bus': None, 'g': 1, 'b': 1.1})
        self._type = 'Shunt'
        self._name = 'Shunt'
        self._bus = {'bus': ['a', 'v']}
        self._params.extend(['fn', 'g', 'b'])
        self._y = ['g', 'b']

    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        #for key, value in zip(self.a, self.Pg):
        #    dae.g[key] += 0
        #for key, value in zip(self.v, self.Pg):
         #   dae.g[key] += 0