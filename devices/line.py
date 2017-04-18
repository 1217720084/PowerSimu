"""

"""
from devices.base_device import base_device
import system


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
                    self.__dict__[self._bus[index][1]].append(system.Bus.a[idx]+system.DAE.n_bus)
