from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, cos, sin, sparse, mul,exp,sparse
from numpy import multiply
import numpy as np


class statcom(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Sn': 1, 'bus': None, 'Vn': 1, 'Imax': 1.2, 'Imin': 0.95, 'fn': 50, 'Kr':1,'Tr':1 ,'gsh':0.1,'bsh':0.1,'Vref':1.0})
        self._type = 'Statcom'
        self._name = 'Statcom'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._states = ['Vsh','Vsha']
        # self._voltages = ['Vref']
        # self.y = ['gsh', 'bsh']
        self._params.extend(['Sn', 'Vn', 'Imax', 'Imin', 'fn', 'Kr','Tr','gsh','bsh','Vref'])

    def gcall(self):
        for i in range(system.Statcom.n):
            system.DAE.g[self.a[i]] += system.DAE.y[self.v[i]]*system.DAE.y[self.v[i]]*self.gsh[i] - \
                                       system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                       (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                        self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            system.DAE.g[self.v[i]] += -system.DAE.y[self.v[i]]*system.DAE.y[self.v[i]]*self.bsh[i] - \
                                       system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                       (self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                        self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            system.DAE.g[2*system.Bus.n+2*i] = system.DAE.y[2*system.Bus.n+2*i+1]*system.DAE.y[2*system.Bus.n+2*i+1]*self.gsh[i] - \
                                               system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                               (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            system.DAE.g[2*system.Bus.n+2*i+1] = system.DAE.y[self.v[i]] - self.Vref[i]
            # print(system.DAE.g)

    def Gycall(self):
        # for i in range(system.Bus.n):
        #     for j in range(2*system.Statcom.n):
        #         Gy[i][system.Bus.n+j] = 0
        for i in range(system.Statcom.n):
            system.DAE.Gy[2*system.Bus.n+2*i, :] = 0
            system.DAE.Gy[:, 2*system.Bus.n+2*i] = 0
            system.DAE.Gy[2*system.Bus.n+2*i+1, :] = 0
            system.DAE.Gy[:, 2*system.Bus.n+2*i+1] = 0
            # system.DAE.Gy[2 * system.Bus.n + 2 * i, 2 * system.Bus.n + 2 * i] = 1
            # system.DAE.Gy[2 * system.Bus.n + 2 * i + 1, 2 * system.Bus.n + 2 * i + 1] = 1
        for i in range(system.Statcom.n):
            #Pi/vi
            system.DAE.Gy[self.a[i],self.v[i]] += 2.0 * system.DAE.y[self.v[i]]*self.gsh[i] - system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                  (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                                   self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[self.a[i], self.v[i]] = system.DAE.Gy[self.a[i],self.v[i]]*system.DAE.y[self.v[i]]
            #Qi/vi
            system.DAE.Gy[self.v[i],self.v[i]] += -2.0 * system.DAE.y[self.v[i]] * self.bsh[i] - system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                  (self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                   self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[self.v[i], self.v[i]] = system.DAE.Gy[self.v[i],self.v[i]]*system.DAE.y[self.v[i]]
            #PE/vi
            system.DAE.Gy[2*system.Bus.n+2*i,self.v[i]] += -system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                         (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                          self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[2 * system.Bus.n + 2 * i, self.v[i]] = system.DAE.Gy[2*system.Bus.n+2*i,self.v[i]]*system.DAE.y[self.v[i]]
            #Vi/vi
            system.DAE.Gy[2*system.Bus.n+2*i+1,self.v[i]] += 1
            #system.DAE.Gy[2 * system.Bus.n + 2 * i + 1, self.v[i]] = system.DAE.Gy[2*system.Bus.n+2*i+1,self.v[i]]*system.DAE.y[self.v[i]]
            #Pi/ai
            system.DAE.Gy[self.a[i],self.a[i]] -= system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                   (-self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                                    self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #Qi/ai
            system.DAE.Gy[self.v[i],self.a[i]] -= system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                  (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                                   self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            # PE/ai
            system.DAE.Gy[2 * system.Bus.n + 2*i,self.a[i]] -= system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                               (-self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                               self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            # Vi/ai
            system.DAE.Gy[2 * system.Bus.n + 2*i+1,self.a[i]] += 0
            #Pi/ash
            system.DAE.Gy[self.a[i],2 * system.Bus.n + 2 * i] -= system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                                 (self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                                  self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #Pi/vsh
            system.DAE.Gy[self.a[i],2 * system.Bus.n + 2 * i + 1] -= system.DAE.y[self.v[i]] * \
                                                                    (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                                                     self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[self.a[i], 2 * system.Bus.n + 2 * i + 1] = system.DAE.Gy[self.a[i],2 * system.Bus.n + 2 * i + 1]*system.DAE.y[2*system.Bus.n+2*i+1]
            #Qi/ash
            system.DAE.Gy[self.v[i],2 * system.Bus.n + 2 * i] -= system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                                 ((-self.gsh[i])*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                                  self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #Qi/vsh
            system.DAE.Gy[self.v[i],2 * system.Bus.n + 2 * i + 1] -= system.DAE.y[self.v[i]] * \
                                                                     (self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                                     self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[self.v[i], 2 * system.Bus.n + 2 * i + 1] = system.DAE.Gy[self.v[i],2 * system.Bus.n + 2 * i + 1]*system.DAE.y[2 * system.Bus.n + 2 * i + 1]
            #PE/ash
            system.DAE.Gy[2 * system.Bus.n + 2 * i,2 * system.Bus.n + 2 * i] -= system.DAE.y[self.v[i]] * system.DAE.y[2*system.Bus.n+2*i+1] * \
                                                                                (self.gsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) + \
                                                                                 self.bsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #PE/vsh
            system.DAE.Gy[2 * system.Bus.n + 2 * i,2 * system.Bus.n + 2 * i + 1] += 2.0 * system.DAE.y[2*system.Bus.n+2*i+1]*self.gsh[i] - system.DAE.y[self.v[i]] * \
                                                                                   (self.gsh[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]) - \
                                                                                    self.bsh[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[2*system.Bus.n+2*i]))
            #system.DAE.Gy[2 * system.Bus.n + 2 * i, 2 * system.Bus.n + 2 * i + 1] = system.DAE.Gy[2 * system.Bus.n + 2 * i,2 * system.Bus.n + 2 * i + 1]*system.DAE.y[2*system.Bus.n+2*i+1]
            #Vi/ash
            system.DAE.Gy[2 * system.Bus.n + 2 * i + 1,2 * system.Bus.n + 2 * i] += 0
            #Vi/vsh
            system.DAE.Gy[2 * system.Bus.n + 2 * i + 1,2 * system.Bus.n + 2 * i + 1] += 0
            #system.DAE.Gy[2 * system.Bus.n + 2 * i + 1, 2 * system.Bus.n + 2 * i + 1] = system.DAE.Gy[2 * system.Bus.n + 2 * i + 1,2 * system.Bus.n + 2 * i + 1]*system.DAE.y[2 * system.Bus.n + 2 * i + 1]







