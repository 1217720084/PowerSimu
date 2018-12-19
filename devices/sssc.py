from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, mul, sparse, cos, sin, exp
import cvxopt.blas
from numpy import multiply
import numpy as np
import math


class sssc(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({ 'L': 0,'bus': None, 'gse':0.1,'bse':-0.1,'Pref':1.0})
        self._type = 'Sssc'
        self._name = 'Sssc'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._states = ['Vse', 'Vsea']
        self._params.extend(['L','bus','gse','bse','Pref'])

    def adbus(self):

        case = open('text_d_004_sssc.txt')
        for i in range(self.n):
            for each_line in case:
                data = each_line.split()
                if data[0] == 'Bus':
                    if 'Bus_'+str(data[1]) == system.Sssc.bus[i]:
                        bus = data[0] + '_' + str(system.Bus.n+i+1)
                        Vb = float(data[2])
                        bus_case = {'bus': bus, 'Vb': Vb}
                        system.Bus.add(idx=bus, **bus_case)
                        break
    def gcall(self):
        for i in range(system.Sssc.n):
            #Pi
            system.DAE.g[self.a[i]] -= system.DAE.y[self.v[i]]*system.DAE.y[self.v[i]]*self.gse[i] - \
                                       system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) + \
                                        self.bse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) + \
                                        self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Qi
            system.DAE.g[self.v[i]] -= -system.DAE.y[self.v[i]]*system.DAE.y[self.v[i]]*self.bse[i] - \
                                       system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) - \
                                        self.bse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) - \
                                        self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Pj
            system.DAE.g[system.Bus.n-system.Sssc.n+i] -= system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.gse[i] - \
                                                          system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                          system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Qj
            system.DAE.g[2*system.Bus.n-system.Sssc.n+i] -= -system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.bse[i] - \
                                                          system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                          system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                          (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #PE
            system.DAE.g[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] = -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                            (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                             self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i])) + \
                                                                             system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                             (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                              self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] -system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F
            system.DAE.g[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] = system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.gse[i] - \
                                                                                system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                                                (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                                                 self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                                                system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                                 self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i])) - \
                                                                                 self.Pref[i]

    def Gycall(self):
        for i in range(system.Sssc.n):
            system.DAE.Gy[2*system.Bus.n+2*system.Statcom.n+2*i, :] = 0
            system.DAE.Gy[:, 2*system.Bus.n+2*system.Statcom.n+2*i] = 0
            system.DAE.Gy[2*system.Bus.n+2*system.Statcom.n+2*i+1, :] = 0
            system.DAE.Gy[:, 2*system.Bus.n+2*system.Statcom.n+2*i+1] = 0
        for i in range(system.Sssc.n):
            #Pi/vi
            system.DAE.Gy[self.a[i],self.v[i]] -= 2*system.DAE.y[self.v[i]]*self.gse[i] - system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) + \
                                        self.bse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) + \
                                        self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Qi/vi
            system.DAE.Gy[self.v[i],self.v[i]] -= -2*system.DAE.y[self.v[i]]*self.bse[i] - system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) - \
                                        self.bse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) - \
                                        self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Pj/vi
            system.DAE.Gy[2*system.Bus.n-system.Sssc.n+i, self.v[i]] -= -system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            #Qi/vj
            system.DAE.Gy[system.Bus.n - system.Sssc.n + i, self.v[i]] -= -system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            #PE/vi
            system.DAE.Gy[2*system.Bus.n+2*system.Statcom.n+2*i,self.v[i]] += -system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                            (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                             self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F/vi
            system.DAE.Gy[2*system.Bus.n+2*system.Statcom.n+2*i+1,self.v[i]] += -system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                                                (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                                                 self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            #Pi/ai
            system.DAE.Gy[self.a[i],self.a[i]] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (-self.gse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) + \
                                        self.bse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (-self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) + \
                                        self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Qi/ai
            system.DAE.Gy[self.v[i],self.a[i]] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) + \
                                        self.bse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i])) - \
                                        system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                       (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) + \
                                        self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Pj/ai
            system.DAE.Gy[system.Bus.n-system.Sssc.n+i, self.a[i]] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            #Qj/ai
            system.DAE.Gy[2*system.Bus.n - system.Sssc.n + i, self.a[i]] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (-self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            # PE/ai
            system.DAE.Gy[2 * system.Bus.n + 2*system.Statcom.n+2*i,self.a[i]] += -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                            (-self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                             self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            # F/ai
            system.DAE.Gy[2 * system.Bus.n + 2*system.Statcom.n+2*i+1,self.a[i]] += -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                                                (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                                                 self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]))
            #Pi/vj
            system.DAE.Gy[self.a[i], 2*system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]] * \
                                                                           (self.gse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) + \
                                                                            self.bse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]))
            #Qi/vj
            system.DAE.Gy[self.v[i], 2 * system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]] * \
                                                                           (self.gse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) - \
                                                                            self.bse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]))
            #Pj/vj
            system.DAE.Gy[system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n - system.Sssc.n + i] -= 2*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.gse[i] - system.DAE.y[self.v[i]] * \
                                                                                                  (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                                                                   self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                                                                  system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                  (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                                                   self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Qj/vj
            system.DAE.Gy[2*system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n - system.Sssc.n + i] -= -2*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.bse[i] -system.DAE.y[self.v[i]] * \
                                                                                                      (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) - \
                                                                                                       self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                                                                      system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                      (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                                                       self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #PE/vj
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i, 2 * system.Bus.n - system.Sssc.n + i] += system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * \
                                     (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                      self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] -system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F/vj
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1, 2 * system.Bus.n - system.Sssc.n + i] += 2*system.DAE.y[2*system.Bus.n-system.Sssc.n+i]*self.gse[i] - system.DAE.y[self.v[i]] * \
                                                                                                                      (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                                                                                       self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                                                                                      system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                                      (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                                                                       self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Pi/aj
            system.DAE.Gy[self.a[i], system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (self.gse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) - \
                                        self.bse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]))
            #Qi/aj
            system.DAE.Gy[self.v[i], system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                       (-self.gse[i]*cos(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]) - \
                                        self.bse[i]*sin(system.DAE.y[self.a[i]]-system.DAE.y[system.Bus.n-system.Sssc.n+i]))
            #Pj/aj
            system.DAE.Gy[system.Bus.n - system.Sssc.n + i, system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (-self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                          system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                          (-self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                           self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Qj/aj
            system.DAE.Gy[2*system.Bus.n - system.Sssc.n + i, system.Bus.n - system.Sssc.n + i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                          system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #PE/aj
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i, system.Bus.n - system.Sssc.n + i] += system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                             (-self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                              self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] -system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F/aj
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1, system.Bus.n - system.Sssc.n + i] += -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n-system.Sssc.n+i] * \
                                                                                                                  (-self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i]) + \
                                                                                                                   self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i]-self.a[i])) + \
                                                                                                                  system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                                  (-self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                                                                   self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Pi/vse
            system.DAE.Gy[self.a[i], 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] -= -system.DAE.y[self.v[i]] * \
                                                                                           (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) + \
                                                                                            self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Qi/vse
            system.DAE.Gy[self.v[i], 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] -= -system.DAE.y[self.v[i]] * \
                                                                                           (self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) - \
                                                                                            self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Pj/vse
            system.DAE.Gy[system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] -= system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * \
                                      (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                       self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Qj/vse
            system.DAE.Gy[2*system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] -= system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * \
                                      (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                       self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #PE/vse
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] += -system.DAE.y[self.v[i]] * \
                                        (self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                         self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i])) + \
                                         system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * \
                                         (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                          self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] -system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F/vse
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] += system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * \
                                                          (self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                           self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))

            #Pi/ase
            system.DAE.Gy[self.a[i], 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                       (self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) - \
                                                                        self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Qi/ase
            system.DAE.Gy[self.v[i], 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] -= -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                       (-self.gse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]) - \
                                                                        self.bse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i]))
            #Pj/ase
            system.DAE.Gy[system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] -= system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                      (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                                                       self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #Qj/ase
            system.DAE.Gy[2*system.Bus.n - system.Sssc.n + i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] -=system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                      (-self.gse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                                                       self.bse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #PE/ase
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] += -system.DAE.y[self.v[i]]*system.DAE.y[2*system.Bus.n+2*system.Statcom.n+2*i+1] * \
                                                                            (self.gse[i] * sin(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                             self.bse[i] * cos(system.DAE.y[self.a[i]] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i])) + \
                                                                             system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                             (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) + \
                                                                              self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] -system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))
            #F/ase
            system.DAE.Gy[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1, 2 * system.Bus.n + 2 * system.Statcom.n + 2 * i] += system.DAE.y[2 * system.Bus.n - system.Sssc.n + i] * system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i + 1] * \
                                                                                                      (self.gse[i] * sin(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]) - \
                                                                                                       self.bse[i] * cos(system.DAE.y[system.Bus.n - system.Sssc.n + i] - system.DAE.y[2 * system.Bus.n + 2 * system.Statcom.n + 2 * i]))