"""

"""
from devices.base_device import base_device
import system


class pv(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pg': 1, 'qg': 0, 'bus': None, 'qgmax': 6, 'qgmin': -6, 'V0': 1.05, 'Vmax': 1.1, 'Vmin': 0.95,})
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

    def gcall(self):

        for i in range(system.PV.n):
            system.DAE.g[self.a[i]] -= self.Pg[i]
            #system.DAE.g[self.v[i]] -= 0

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

    def gcall(self):

        for i in range(system.PV.n):
            system.DAE.g[self.a[i]] -= self.Pg[i]
            if ~system.settings.pv2pq:
                system.DAE.g[self.v[i]] -= 0

        #判断PV是否要转为PQ

        if system.settings.iter < system.settings.pv2pqiter:
            return
        else:
            prev_err = 2*system.settings.error

        self.newpq = 0

        # Qmin
        err = [0] * self.n

        for i in range(len(self.n)):
            err[i] = self.qgmin[i] - system.DAE.g[self.v[i]] - prev_err
        max_err = max(err)
        if max_err > 0:
            for i in range(len(self.n)):
                if max_err[i] == err[i]:
                    self.qg[i] = self.qgmin[i]
                    self.pq[i] = 1
                    self.newpq = ~system.settings.multipvswitch

        # Qmax
        err = [0] * self.n

        for i in range(len(self.n)):
            err[i] = self.qgmax[i] - system.DAE.g[self.v[i]] + prev_err
        min_err = min(err)
        if max_err < 0 & ~self.newpq:
            for i in range(len(self.n)):
                if min_err[i] == err[i]:
                    self.qg[i] = self.qgmax[i]
                    self.pq[i] = 1
                    self.newpq = ~system.settings.multipvswitch

        for i in range(system.PV.n):
            system.DAE.g[self.v[i]] -= self.qg[i]
