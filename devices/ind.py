"""

"""

import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, sin, cos, spmatrix
import numpy as np
import math

class ind3(base_device):

    def __init__(self):
        base_device.__init__(self)
        self._data.update({'bus': None, 'fn': 60, 'type': 3, 'Sn': 100, 'Vn': 230, 'sup': 0, 'rs': 0, 'xs': 0,
                           'rr1': 0.18, 'xr1': 0.12, 'rr2': 0.12, 'xr2': 0.001, 'xm': 0.04, 'Hm': 3.5, 'a1': 1,
                           'b': 0.51, 'c': 0, 'tup': 0.51, 'allow': 0})
        self._type = 'Ind'
        self._name = 'Ind'
        self._bus = {'bus': ['a', 'v']}
        self.n = 0
        self._algebs = []
        self._states = ['slip', 'e1r', 'e1m']
        self._params.extend(['fn',  'Sn', 'Vn', 'sup', 'rs', 'xs',
                            'rr1', 'xr1', 'rr2', 'xr2', 'xm', 'Hm', 'a1',
                            'b', 'c', 'tup', 'allow'])

        self._z = ['rs', 'xs', 'rr1', 'xr1', 'rr2', 'xr2', 'xm']
        self._powers = ['Hm', 'a1', 'b', 'c']

        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def setdata(self):

        Wn = 2 * math.pi * system.Settings.freq

        # 转矩系数：A+B*slip+C*slip^2

        self.A = self.a1 + self.b + self.c
        self.B = -self.b - 2 * self.c
        self.C = self.c

        # 1/2*Hm
        self.i2Hm = div(1, 2 * self.Hm)

        # x0,x',x''
        self.x0 = self.xs + self.xm
        self.x1 = self.xs + div(mul(self.xr1, self.xm), self.xr1 + self.xm)
        self.x2 = self.xs + div(mul(self.xr1, self.xm, self.xr2), mul(self.xr1, self.xm)+mul(self.xr2, self.xm)+mul(self.xr1, self.xr2))

        # T'0,T''0
        for i in range(self.n):
            if self.rr1[i] == 0:
                self.rr1[i] = 1
            if self.rr2[i] == 0:
                self.rr2[i] = 1

        self.T10 = div(self.xr1 + self.xm, mul(self.rr1, Wn))
        self.T20 = div(self.xr2+div(mul(self.xr1, self.xm), self.xr1+self.xm), mul(self.rr2, Wn))

        # 1/xm, x's = xs+xr1

        for i in range(self.n):
            if self.xm[i] == 0:
                print('第%i 台感应电机的磁抗为0' % i)
                self.xm[i] = 1
            if self.tup[i] < 0:
                self.tup = 0

        self.ixm = div(1, self.xm)
        self.x1s = self.xs + self.xr1

        if self.sup == 0:
            self.tup = 0



    def gcall(self):

        if self.n == 0:
            return

        slip = system.DAE.x[self.slip]

        V = mul(self.u, system.DAE.y[self.v]) # 后续得修改
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        Vr = mul(V, st)
        Vm = mul(V, ct)

        e1r = system.DAE.x[self.e1r]
        e1m = system.DAE.x[self.e1m]

        a03 = mul(self.rs, self.rs) + mul(self.x1, self.x1)
        a13 = div(self.rs, a03)
        a23 = div(self.x1, a03)
        a33 = self.x0 - self.x1

        Im = mul(-a23, -Vr-e1r+mul(a13, (Vm-e1m)))
        Ir = mul(a13, -Vr-e1r) + mul(a23, Vm-e1m)

        system.DAE.g = system.DAE.g \
                       + spmatrix(mul(-Vr, Ir)+mul(Vm, Im), self.a, matrix(0, (self.n, 1)), (system.DAE.ny, 1)) \
                       + spmatrix(mul(Vr, Im) + mul(Vm, Ir), self.v, matrix(0, (self.n, 1)), (system.DAE.ny, 1))

    def Gycall(self):

        if self.n == 0:
            return

        slip = system.DAE.x[self.slip]

        V = mul(self.u, system.DAE.y[self.v])  # 后续得修改
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        Vr = mul(V, st)
        Vm = mul(V, ct)

        e1r = system.DAE.x[self.e1r]
        e1m = system.DAE.x[self.e1m]

        a03 = mul(self.rs, self.rs) + mul(self.x1, self.x1)
        a13 = div(self.rs, a03)
        a23 = div(self.x1, a03)
        a33 = self.x0 - self.x1

        Im = mul(-a23, -Vr - e1r + mul(a13, (Vm - e1m)))
        Ir = mul(a13, -Vr - e1r) + mul(a23, Vm - e1m)

        Am = mul(a23, Vm) - mul(a13, Vr)
        Ar = mul(-a13, Vm) - mul(a23, Vr)
        Bm = mul(a23, st) + mul(a13, ct)
        Br = mul(-a13, st) + mul(a23, ct)

        system.DAE.Gy = system.DAE.Gy \
                       + spmatrix(mul(-Vm, Ir)-mul(Vr, Ar)-mul(Vr, Im)+mul(Vm, Am), self.a, self.a, (system.DAE.ny, system.DAE.ny)) \
                       + spmatrix(mul(-st, Ir)-mul(Vr, Br)+mul(ct, Im)+mul(Vm, Bm), self.a, self.v, (system.DAE.ny, system.DAE.ny)) \
                       + spmatrix(mul(Vm, Im)+mul(Vr, Am)-mul(Vr, Ir)+mul(Vm, Ar), self.v, self.a,
                                   (system.DAE.ny, system.DAE.ny)) \
                        + spmatrix(mul(st, Im)+mul(Vr, Bm)+mul(ct, Ir)+mul(Vm, Br), self.v, self.v,
                                   (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        slip = system.DAE.x[self.slip]

        u = self.u  # 后续得修改

        Wn = 2 * math.pi * system.Settings.freq * u

        V = mul(self.u, system.DAE.y[self.v])  # 后续得修改
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        Vr = mul(V, st)
        Vm = mul(V, ct)

        i2Hm = mul(u, self.i2Hm)

        Tm = self.A + mul(slip, self.B+mul(slip, self.C))

        e1r = system.DAE.x[self.e1r]
        e1m = system.DAE.x[self.e1m]

        a03 = mul(self.rs, self.rs) + mul(self.x1, self.x1)
        a13 = div(self.rs, a03)
        a23 = div(self.x1, a03)
        a33 = self.x0 - self.x1

        Im = mul(-a23, -Vr - e1r + mul(a13, (Vm - e1m)))
        Ir = mul(a13, -Vr - e1r) + mul(a23, Vm - e1m)

        A = mul(-Wn, slip, e1r) - div(e1m-mul(a33, Ir), self.T10)
        system.DAE.f[self.slip] = mul(Tm-mul(e1r, Ir)-mul(e1m, Im), i2Hm)
        system.DAE.f[self.e1r] = mul(Wn, slip, e1m) - div(e1r+mul(a33, Im), self.T10)
        system.DAE.f[self.e1m] = mul(-Wn, slip, e1r) - div(e1m-mul(a33, Ir), self.T10)


    def Fxcall(self):

        if self.n == 0:
            return

        slip = system.DAE.x[self.slip]

        u = self.u  # 后续得修改
        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if slip[i] < 1 or self.allow[i]:
                z[i] = 1
            else:
                z[i] = 0

        Wn = 2 * math.pi * system.Settings.freq * u

        V = mul(self.u, system.DAE.y[self.v])  # 后续得修改
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        Vr = mul(V, st)
        Vm = mul(V, ct)

        i2Hm = mul(u, self.i2Hm, z)

        km = self.B + mul(2*self.C, slip)

        e1r = system.DAE.x[self.e1r]
        e1m = system.DAE.x[self.e1m]

        a03 = mul(self.rs, self.rs) + mul(self.x1, self.x1)
        a13 = div(self.rs, a03)
        a23 = div(self.x1, a03)
        a33 = self.x0 - self.x1

        Im = mul(-a23, -Vr - e1r + mul(a13, (Vm - e1m)))
        Ir = mul(a13, -Vr - e1r) + mul(a23, Vm - e1m)

        Am = mul(a23, Vm) - mul(a13, Vr)
        Ar = mul(-a13, Vm) - mul(a23, Vr)
        Bm = mul(a23, st) + mul(a13, ct)
        Br = mul(-a13, st) + mul(a23, ct)

        is3 = self.slip
        er3 = self.e1r
        em3 = self.e1m
        iu = matrix(0, (self.n, 1))
        for i in range(self.n):
            if self.u[i] == 1:
                iu[i] = 0
            else:
                iu[i] = 1


        system.DAE.Fx = system.DAE.Fx \
                        + spmatrix(mul(km, i2Hm)-iu, is3, is3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(-Ir+mul(e1r, a13)-mul(e1m, a23), i2Hm), is3, er3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(-Im+mul(e1r, a23)-mul(e1m, a13), i2Hm), is3, em3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(z, Wn, e1m), er3, is3, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(1+mul(a33, a23), self.T10), er3, er3, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(Wn, slip)+div(mul(a33, a13), self.T10), er3, em3, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(mul(z, Wn, e1r), em3, is3, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(mul(Wn, slip)+div(mul(a33, a13), self.T10), em3, er3, (system.DAE.nx, system.DAE.nx)) \
                        - spmatrix(div(1+mul(a33, a23), self.T10), em3, em3, (system.DAE.nx, system.DAE.nx))



        system.DAE.Gx = system.DAE.Gx \
                        + spmatrix(mul(a13, Vr)+mul(a23, Vm), self.a, er3, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(a23, Vr) - mul(a13, Vm), self.a, em3, (system.DAE.ny, system.DAE.nx)) \
                        + spmatrix(mul(a23, Vr) - mul(a13, Vm), self.v, er3, (system.DAE.ny, system.DAE.nx)) \
                        - spmatrix(mul(a13, Vr) + mul(a23, Vm), self.v, em3, (system.DAE.ny, system.DAE.nx))



        system.DAE.Fy = system.DAE.Fy \
                        - spmatrix(mul(mul(e1r, Ar)+mul(e1m, Am), i2Hm), is3, self.a, (system.DAE.nx, system.DAE.ny)) \
                        - spmatrix(div(mul(a33, Am), self.T10), er3, self.a, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(div(mul(a33, Ar), self.T10), em3, self.a, (system.DAE.nx, system.DAE.ny)) \
                        - spmatrix(mul(mul(e1r, Br) + mul(e1m, Bm), i2Hm), is3, self.v, (system.DAE.nx, system.DAE.ny)) \
                        - spmatrix(div(mul(a33, Bm), self.T10), er3, self.v, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(div(mul(a33, Br), self.T10), em3, self.v, (system.DAE.nx, system.DAE.ny))



    def xinit(self):

        if self.n == 0:
            return

        system.DAE.x[self.e1r] = 0.05
        system.DAE.x[self.slip] = 0
        system.DAE.x[self.e1m] = 0.9

    def windup(self, type):

        if type == 'td':
            idx = []
            for i in range(self.n-1):
                if self.allow[i] == 0:
                    idx[i] = self.slip[i]
            xmax = 1
            xmin = -1e3
            x = system.DAE.x[idx]

            for i in range(len(idx)):
                if x[idx[i]] >= xmax or x[idx[i]] <= xmin:
                    if system.DAE.f[idx[i]] == 0:

                        system.DAE.tn[idx[i]] = 0
                        system.DAE.Ac[idx[i], :] = 0
                        system.DAE.Ac[:, idx[i]] = 0
                        system.DAE.Ac = system.DAE.Ac - spmatrix(1.0, [idx[i]], [idx[i]], (
                        (system.DAE.nx + system.DAE.ny, system.DAE.nx + system.DAE.ny)))











