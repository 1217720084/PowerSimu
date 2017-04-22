"""

"""
from devices.base_device import base_device
import system
from cvxopt.base import matrix


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

        zero = [0] * system.Bus.n
        zeros = [zero] * system.Bus.n
        system.DAE.Y = zeros
        Line_n = len(self.r)
        i = 0
        while i < Line_n:
            if system.Line.kT[i] == 0:
                system.DAE.Y[self.af[i]][self.at[i]] += -complex(self.r[i], self.x[i])
                print(system.DAE.Y[1])
                system.DAE.Y[self.at[i]][self.af[i]] += -complex(self.r[i], self.x[i])
                system.DAE.Y[self.af[i]][self.af[i]] += complex(self.r[i], self.x[i]) + complex(0, self.b[i] / 2)
                system.DAE.Y[self.at[i]][self.at[i]] += complex(self.r[i], self.x[i]) + complex(0, self.b[i] / 2)
            else:
                system.DAE.Y[self.af[i]][self.at[i]] += -complex(self.r[i], self.x[i]) * \
                                                         self.kT[i]
                system.DAE.Y[self.af[i]][self.af[i]] += complex(self.r[i], self.x[i]) * self.kT[i] + complex(self.r[i], self.x[i]) * \
                                                        self.kT[i] / (self.kT[i] - 1)
                system.DAE.Y[self.at[i]][self.at[i]] += complex(self.r[i], self.x[i]) * self.kT[i] + complex(self.r[i], self.x[i]) * \
                                                          self.kT[i] ** 2 / (1 - self.kT[i])
            i += 1
            # print(matrix(system.DAE.Y))
        # for i in range(Line_n):
        #     if system.Line.kT[i] == 0:
        #       system.DAE.Y[self.af[i]][self.at[i]] += - 1/(complex(self.r[i], self.x[i]))
        #       print(system.DAE.Y[10])
        #       system.DAE.Y[self.at[i]][self.af[i]] += -complex(self.r[i], self.x[i])
        #       print(system.DAE.Y[24])
        #       system.DAE.Y[self.af[i]][self.af[i]] += complex(self.r[i], self.x[i]) + complex(0, self.b[i] / 2)
        #       system.DAE.Y[self.at[i]][self.at[i]] += complex(self.r[i], self.x[i]) + complex(0, self.b[i] / 2)




        system.DAE.Y_G = matrix.real(matrix(system.DAE.Y))
        system.DAE.Y_B = matrix.imag(matrix(system.DAE.Y))

        system.DAE.Y = matrix(system.DAE.Y)

        # for i in range(system.Bus.n):
        #     system.DAE.Y = matrix([zero, system.DAE.Y[i]])
