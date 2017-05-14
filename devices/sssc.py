from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, mul,sparse
from numpy import multiply
import numpy as np

class sssc(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Sn': 1, 'bus': None, 'Vn': 1, 'Imax': 1.2, 'Imin': 0.95, 'fn': 50, 'Kr': 1, 'Tr': 1})
        self._type = 'Sssc'
        self._name = 'Sssc'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._states = ['Vsh', 'Vsha']
        self._params.extend(['Sn', 'bus', 'Vn', 'Imax', 'Imin', 'fn', 'Kr', 'Tr'])