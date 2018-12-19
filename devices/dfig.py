#!/usr/bin/env Python
# coding=utf-8
from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, mul, div, sin, cos, sparse
from cvxopt import umfpack
import numpy
import math
import random
from scipy import interpolate
import decimal
from decimal import Decimal
import time



class dfig(base_device):

    def __init__(self):

        base_device.__init__(self)
        self._data.update(
            {'fn': 50, 'bus': None, 'wind': None, 'rs': 0.0101, 'xs': 0.0747, 'rr': 0.0031, 'xr': 0.0098, 'xm': 5.6903,
             'Hm': 3, 'Kp': 200, 'Tp': 2, 'Kv': 0.1, 'Te': 0.01, 'R': 37.5, 'np': 4, 'nb': 3, 'nGB': 1/104.5,
             'pmax': 1.2, 'pmin': -0.1, 'qmax': 0.8, 'qmin': -0.7, 'ng': 1, 'u': 1})
        self._type = 'Dfig'
        self._name = 'Dfig'
        self.n = 0

        self._data.update({'wind': 0})
        self._bus = {'bus': ['a', 'v']}

        self._algebs = ['pwa', 'vref']
        self._states = ['omega_m', 'theta_p', 'idr', 'iqr']
        self._z = ['rs', 'xs', 'rr', 'xr', 'xm']
        self._powers = ['Hm', 'pmax', 'pmin', 'qmax', 'qmin']
        self._params.extend(['fn', 'rs', 'xs', 'rr', 'xr', 'xm', 'Sn', 'Vn',
                             'Hm', 'Kp', 'Tp', 'Kv', 'Te', 'R', 'np', 'nb', 'nGB', 'pmax',
                             'pmin', 'qmax', 'qmin', 'ng'])
        # self._params.extend(['Pg', 'qgmax', 'qgmin', 'V0', 'Vmax', 'Vmin'])
        # self.n_PV = 0
        # self._voltages = ['V0']
        # self._powers = ['Pg']
        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def setx0(self):
        if self.n == 0:
            return
        check = 1

        Pc = system.Bus.Pg[self.a]
        Qc = system.Bus.Qg[self.a]
        Vc = system.DAE.y[self.v]
        ac = system.DAE.y[self.a]

        vds = mul(-Vc, sin(ac))   # 角度还是弧度
        vqs = mul(Vc, cos(ac))

        rho = system.Wind.rho[self.wind]

        ones = matrix(1, (self.n, 1))
        # 常数
        # xs + xm
        self.dat1 = self.xs + self.xm
        self.dat2 = self.xr + self.xm
        self.dat3 = div(ones, mul(2*ones, self.Hm))
        self.dat4 = div(mul(4*math.pi*system.Settings.freq*self.R, self.nGB), self.np)
        self.dat5 = math.pi*mul(self.R, self.R)
        self.dat6 = Vc
        self.dat8 = mul(div(-self.pmin, self.xm), self.dat1)
        self.dat9 = mul(div(-self.pmax, self.xm), self.dat1)
        self.dat10 = -mul(div(div(ones, self.xm) + self.qmin, self.xm), self.dat1)
        self.dat11 = -mul(div(div(ones, self.xm) + self.qmax, self.xm), self.dat1)

        # 初始化状态变量
        for i in range(self.n):
            # 参数
            Vds = vds[i]
            Vqs = vqs[i]
            Rs = self.rs[i]
            Rr = self.rr[i]
            Xm = self.xm[i]
            x1 = self.dat1[i]
            x2 = self.dat2[i]
            Pg = Pc[i]
            Qg = Qc[i]

            # 转子速度
            Pci = Pc[i]*system.Settings.mva
            if (Pci < self.Sn[i]) & (Pc[i] > 0):
                omega = 0.5*Pc[i]*system.Settings.mva/self.Sn + 0.5
            elif Pci >= self.Sn[i]:
                omega = 1
            else:
                omega = 0.5

            slip = 1 - omega

            iqr = -x1*self.Sn*(2*omega-1)/Vc[i]/Xm/system.Settings.mva/omega
            A = [[-Rs, Vqs], [x1, -Vds]]
            B = [Vds-Xm*iqr, Qg]
            A = sparse(A)
            B = matrix(B)
            umfpack.linsolve(A, B)

            # B = numpy.array(B)
            # Is = numpy.linalg.solve(A, B)
            ids = B[0]
            iqs = B[1]
            idr = -(Vqs+Rs*iqs+x1*ids)/Xm
            vdr = -Rr*idr+slip*(x2*iqr+Xm*iqs)
            vqr = -Rr*iqr-slip*(x2*idr+Xm*ids)

            jac = matrix(0.0, (6, 6))
            eqn = matrix(1.0, (6, 1))
            inc = matrix(1.0, (6, 1))

            x = matrix(0.0, (6, 1))

            x[0] = ids
            x[1] = iqs
            x[2] = idr
            x[3] = iqr
            x[4] = vdr
            x[5] = vqr

            rows = [0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5]
            cols = [0, 1, 3, 0, 1, 2, 2, 4, 3, 5, 0, 1, 2]
            vals = [-Rs, x1, Xm, -x1, -Rs, -Xm, -Rr, -1, -Rr, -1, Vds, Vqs, -Xm * Vc[i] / x1]
            jac0 = spmatrix(vals, rows, cols, (6, 6), 'd')

            # jac0 = jac + spmatrix(-Rs, [0], [0], (6, 6))
            # jac0 = jac + spmatrix(x1, [0], [1], (6, 6))
            # jac0 = jac + spmatrix(Xm, [0], [3], (6, 6))
            # jac0 = jac + spmatrix(-x1, [1], [0], (6, 6))
            # jac0 = jac + spmatrix(-Rs, [1], [1], (6, 6))
            # jac0 = jac + spmatrix(-Xm, [1], [2], (6, 6))
            # jac0 = jac + spmatrix(-Rr, [2], [2], (6, 6))
            # jac0 = jac + spmatrix(-1, [2], [4], (6, 6))
            # jac0 = jac + spmatrix(-Rr, [3], [3], (6, 6))
            # jac0 = jac + spmatrix(-1, [3], [5], (6, 6))
            # jac0 = jac + spmatrix(Vds, [4], [0], (6, 6))
            # jac0 = jac + spmatrix(Vqs, [4], [1], (6, 6))
            # jac0 = jac + spmatrix(-Xm*Vc[i]/x1, [5], [2], (6, 6))


            k = x1*self.Sn[i]/Vc[i]/Xm/system.Settings.mva

            iter = 0

            while max(abs(eqn)) > 1e-8:
                if iter > 20:
                    print('双馈风力发电机%i初始化失败' % i)
                    check = 0
                    break
                eqn[0] = -Rs * x[0] + x1 * x[1] + Xm * x[3] - Vds
                eqn[1] = -Rs * x[1] - x1 * x[0] - Xm * x[2] - Vqs
                eqn[2] = -Rr * x[2] + slip * (x2 * x[3] + Xm * x[1]) - x[4]
                eqn[3] = -Rr * x[3] - slip * (x2 * x[2] + Xm * x[0]) - x[5]
                eqn[4] = Vds * x[0] + Vqs * x[1] + x[4] * x[2] + x[5] * x[3] - Pg
                eqn[5] = -Xm * Vc[i] * x[2]/x1 - Vc[i] * Vc[i] / x1 - Qg

                rows = [2, 2, 3, 3, 4, 4, 4, 4]
                cols = [1, 3, 0, 2, 2, 3, 4, 5]
                vals = [slip[i] * Xm, slip[i] * x2, -slip[i] * Xm, -slip[i] * x2, x[4], x[5], x[2], x[3]]

                jac = jac0 + spmatrix(vals, rows, cols, (6, 6), 'd')

                # jac = jac + spmatrix(slip * Xm, [2], [1], (6, 6))
                # jac = jac + spmatrix(slip * x2, [2], [3], (6, 6))
                # jac = jac + spmatrix(-slip * Xm, [3], [0], (6, 6))
                # jac = jac + spmatrix(-slip * x2, [3], [2], (6, 6))
                # jac = jac + spmatrix(x[4], [4], [2], (6, 6))
                # jac = jac + spmatrix(x[5], [4], [3], (6, 6))
                # jac = jac + spmatrix(x[2], [4], [4], (6, 6))
                # jac = jac + spmatrix(x[3], [4], [5], (6, 6))

                jac = sparse(jac)
                umfpack.linsolve(jac, eqn)
                #inc = -numpy.linalg.solve(jac, eqn)
                x = x - eqn
                iter = iter + 1

            ids = x[0]
            iqs = x[1]
            idr = x[2]
            iqr = x[3]
            vdr = x[4]
            vqr = x[5]

            if iqr > self.dat8[i]:
                print('Warning: Dfig %i at Bus %s : iqr is over its max limit' % (i, self.a[i]))
                check = 0
            if iqr < self.dat9[i]:
                print('Warning: Dfig %i at Bus %s : iqr is under its min limit' % (i, self.a[i]))
                check = 0
            if idr > self.dat10[i]:
                print('Warning: Dfig %i at Bus %s : idr is over its max limit' % (i, self.a[i]))
                check = 0
            if idr < self.dat11[i]:
                print('Warning: Dfig %i at Bus %s : idr is under its min limit' % (i, self.a[i]))
                check = 0

            # theta_p
            contex = decimal.getcontext()
            contex.rounding = decimal.ROUND_05UP
            theta = self.Kp * round(Decimal(1000*(omega[i]-1)))/1000
            theta = max([theta[0], 0])

            # wind turbine state variables
            system.DAE.x[self.idr[i]] = idr
            system.DAE.x[self.iqr[i]] = iqr
            system.DAE.x[self.omega_m[i]] = omega
            system.DAE.x[self.theta_p[i]] = theta
            # Vref
            Kv = self.Kv[i]
            if Kv == 0:  # 没有电压控制
                self.dat6 = 0
            else:
                self.dat6 = Vc[i] - (idr+Vc[i]/Xm)/Kv

            self.dat7 = -k*max([min([2*omega[i]-1, 1]), 0])/omega[i] - iqr

            # electric torque
            Tel = Xm*(iqr*ids-idr*iqs)
            if Pg < 0:
                print('** Turbine power is negative at Bus %i' % self.a[i])
                print('Wind speed % i can not be initialized.' % self.wind[i])
                system.DAE.x[system.Wind.vw[self.wind[i]]] = 1
                continue
            # wind power [MW]

            Pw = Tel * omega * system.Settings.mva * 1e6 / self.ng

            # wind speed
            iter = 0
            incvw = matrix(1.0, (1, 1))
            eqnvw = [1]
            R = self.dat4[i]
            AA = self.dat5[i]
            # initial guess wind speed
            vw = 0.9 * system.Wind.vwn[self.wind[i]]
            while abs(incvw[0]) > 1e-7:

                if iter > 50:
                    print('**Initialization of wind speed %i failed(converge problem)' % self.wind[i])
                    print('Tip: Try increasing the nominal wind speed')
                    check = 0
                    break
                output = system.Dfig.windpower(rho[i], vw, AA, R, omega, theta, 1)
                eqnvw = output - Pw
                jacvw = system.Dfig.windpower(rho[i], vw, AA, R, omega, theta, 2)
                # incvw = -numpy.linalg.solve(jacvw[2], eqnvw)
                incvw = -eqnvw/jacvw[1]
                vw = vw + incvw[0]
                iter = iter + 1
            # average initial wind speed
            system.DAE.x[system.Wind.vw[self.wind[i]]] = vw/system.Wind.vwn[self.wind[i]]

            # find and delete static generators
            for bj in range(self.n):
                for bk in range(system.PV.n):
                    if self.u[bj] * self.a[bj] == system.PV.a[bk]:
                        idx = 'PV_' + str(bk + 1)
                        system.PV.remove(idx)
                for bk in range(system.SW.n):
                    if self.u[bj] * self.a[bj] == system.SW.a[bk]:
                        system.SW.remove(idx='SW_1')
        system.DAE.x[self.idr] = mul(self.u, system.DAE.x[self.idr])
        system.DAE.x[self.iqr] = mul(self.u, system.DAE.x[self.iqr])
        system.DAE.x[self.omega_m] = mul(self.u, system.DAE.x[self.omega_m])
        system.DAE.x[self.theta_p] = mul(self.u, system.DAE.x[self.theta_p])
        system.DAE.y[self.vref] = mul(self.u, self.dat6)
        xomega = system.DAE.x[self.omega_m]
        xomega1 = matrix(0.0, (self.n, 1))
        for i in range(self.n):
            xomega1[i] = max([min([2*xomega[i]-1, 1]), 0])
        system.DAE.y[self.pwa] = mul(self.Sn, div(xomega1, system.Settings.mva))

        if not check:
            print('双馈风力发电机初始化失败')
        else:
            print('双馈风力发电机初始化完成')

    def windpower(self, rho, vw, Ar, R, omega, theta, type):

        lambd = div(mul(omega, R), vw)
        lambdi = 1/(1/(lambd+0.08*theta)-0.035/(theta ** 3+1))

        if type == 1:
            cp = 0.22 * mul((116/lambdi-0.4*theta-5), math.exp(-12.5/lambdi[0]))
            output = 0.5 * mul(rho, mul(cp, mul(Ar, vw ** 3)))
        elif type == 2:
            # output = matrix(0.0, (len(omega), 3))
            output = matrix(0.0, (1, 3))
            a1 = math.exp(-12.5/lambdi[0])
            a2 = mul(lambd+0.08*theta, lambd+0.08*theta)
            a3 = 116/lambdi - 0.4*theta - 5
            a4 = -9.28/a2+12.180*div(mul(theta, theta), (theta ** 3+1) ** 2)-0.4
            a5 = 1.000/a2-1.3125*div(mul(theta, theta), (theta ** 3+1) ** 2)

            # d Pw/d omega_m
            output1 = mul(R, mul(a1, mul(rho, mul(vw, mul(vw, mul(Ar, div(-12.760+1.3750*a3, a2)))))))
            # d Pw/d vw
            output2 = mul(mul(omega, mul(R, div((12.760-1.3750*a3, a2))))+0.330*mul(a3, vw), mul(vw, mul(Ar,mul(rho,a1))))
            # d Pw/d theta_p
            output3 = 0.110*mul(rho, mul(a4+mul(a3, a5), mul(a1, mul(Ar, vw ** 3))))
            output[0] = output1
            output[1] = output2
            output[2] = output3
            # output = [output1, output2, output3]
        return output

    def base(self):

        if self.n == 0:
            return

        Vold = self.Vn
        Vbus = system.Bus.Vb[self.a]
        Verr = abs(div(Vbus-Vold, Vbus))

        for i in range(len(Verr)):
            if Verr[i] > 0.1:
                print('节点%i 的Dfig % i 电压额定值与节点电压额定值相差超过0.1' % (self.a[i]+1, i+1))

        base_device.base(self, Sb=system.Settings.mva, Vb=system.Bus.Vb[self.a])

    def gcall(self):

        if self.n == 0:
            return
        omega_m = system.DAE.x[self.omega_m]
        idr = system.DAE.x[self.idr]
        iqr = system.DAE.x[self.iqr]

        pwa = system.DAE.y[self.pwa]
        V = system.DAE.y[self.v]
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        rs = self.rs
        rr = self.rr
        xm = self.xm
        as1 = rs**2 + self.dat1**2
        a13 = div(rs, as1)
        a23 = div(self.dat1, as1)
        a33 = self.dat2

        vds = mul(-V, st)
        vqs = mul(V, ct)

        ids = mul(-a13, vds-mul(xm, iqr)) - mul(a23, vqs+mul(xm, idr))
        iqs = mul(a23, vds - mul(xm, iqr)) - mul(a13, vqs + mul(xm, idr))

        ones = matrix(1.0, (self.n, 1))

        vdr = mul(-rr, idr) + mul(ones-omega_m, mul(a33, iqr)+mul(xm, iqs))
        vqr = mul(-rr, iqr) - mul(ones - omega_m, mul(a33, idr) + mul(xm, ids))

        p = mul(vds, ids) + mul(vqs, iqs) + mul(vdr, idr) + mul(vqr, iqr)
        q = div(mul(-V, mul(xm, idr)+V), self.dat1)

        xomega = system.DAE.x[self.omega_m]
        xomega1 = matrix(0.0, (self.n, 1))
        for i in range(self.n):
            xomega1[i] = max([min([2 * xomega[i] - 1, 1]), 0])

        pwa1 = mul(self.u, mul(self.Sn, div(xomega1, system.Settings.mva)))-pwa
        ones = [0] * self.n
        system.DAE.g = system.DAE.g - \
                       spmatrix(mul(self.u, p), self.a, ones, (system.DAE.ny, 1)) - \
                       spmatrix(mul(self.u, q), self.v, ones, (system.DAE.ny, 1)) + \
                       spmatrix(mul(self.u, self.dat6)-system.DAE.y[self.vref], self.vref, ones, (system.DAE.ny, 1)) + \
                       spmatrix(pwa1, self.pwa, ones, (system.DAE.ny, 1))

    def Gycall(self):

        if self.n == 0:
            return

        omega_m = system.DAE.x[self.omega_m]
        idr = system.DAE.x[self.idr]
        iqr = system.DAE.x[self.iqr]

        pwa = system.DAE.y[self.pwa]
        V = system.DAE.y[self.v]
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        rs = self.rs
        rr = self.rr
        xm = self.xm
        as1 = rs ** 2 + self.dat1 ** 2
        a13 = div(rs, as1)
        a23 = div(self.dat1, as1)
        a33 = self.dat2

        vds = mul(-V, st)
        vqs = mul(V, ct)

        ids = mul(-a13, vds - mul(xm, iqr)) - mul(a23, vqs + mul(xm, idr))
        iqs = mul(a23, vds - mul(xm, iqr)) - mul(a13, vqs + mul(xm, idr))

        vdr = mul(-rr, idr) + mul(1 - omega_m, mul(a33, iqr) + mul(xm, iqs))
        vqr = mul(-rr, iqr) + mul(1 - omega_m, mul(a33, idr) + mul(xm, ids))

        pv = mul(-ids, st) + mul(iqs, ct)
        pt = mul(-ids, vqs) + mul(iqs, vds)

        idsv = mul(a13, st) - mul(a23, ct)
        idst = mul(a13, vqs) - mul(a23, vds)
        iqsv = mul(-a23, st) - mul(a13, ct)
        iqst = mul(-a23, vqs) - mul(a13, vds)

        k = mul(1-omega_m, xm)

        vdrv = mul(k, iqsv)
        vdrt = mul(k, iqst)
        vqrv = mul(-k, idsv)
        vqrt = mul(-k, idst)

        j11 = mul(vds, idst) + mul(vqs, iqst) + mul(vdrt, idr) + mul(vqrt, iqr) + pt
        j12 = mul(vds, idsv) + mul(vqs, iqsv) + mul(vdrv, idr) + mul(vqrv, iqr) + pv
        j22 = div(-(mul(xm, idr)+2*V), self.dat1)

        ones = [0] * self.n

        system.DAE.Gy = system.DAE.Gy - \
                        spmatrix(mul(self.u, j11), self.a, self.a, (system.DAE.ny, system.DAE.ny)) - \
                        spmatrix(mul(self.u, j12), self.a, self.v, (system.DAE.ny, system.DAE.ny)) - \
                        spmatrix(mul(self.u, j22), self.v, self.v, (system.DAE.ny, system.DAE.ny)) - \
                        spmatrix([1] * self.n, self.vref, self.vref, (system.DAE.ny, system.DAE.ny)) - \
                        spmatrix([1] * self.n, self.pwa, self.pwa, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        omega_m = system.DAE.x[self.omega_m]
        theta_p = system.DAE.x[self.theta_p]
        idr = system.DAE.x[self.idr]
        iqr = system.DAE.x[self.iqr]

        vw = system.DAE.x[system.Wind.vw[self.wind]]
        rho = system.Wind.rho[self.wind]

        pwa = system.DAE.y[self.pwa]
        V = system.DAE.y[self.v]
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        rs = self.rs
        xs = self.xs
        rr = self.rr
        xr = self.xr
        xm = self.xm

        i2Hm = mul(self.u, self.dat3)
        Kp = self.Kp
        Tp = self.Tp
        Kv = self.Kv
        Te = self.Te
        R = self.dat4
        A = self.dat5
        ng = self.ng

        as1 = rs ** 2 + self.dat1 ** 2
        a13 = div(rs, as1)
        a23 = div(self.dat1, as1)
        a33 = self.dat2

        vds = mul(-V, st)
        vqs = mul(V, ct)

        ids = mul(-a13, vds - mul(xm, iqr)) - mul(a23, vqs + mul(xm, idr))
        iqs = mul(a23, vds - mul(xm, iqr)) - mul(a13, vqs + mul(xm, idr))

        # wind speed in m/s
        notu = matrix(0, (self.n, 1))
        for i in range(self.n):
            if self.u[i] == 1:
                notu[i] = 0
            else:
                notu[i] = 1
        iomega = 1/(omega_m+notu)
        Vw = mul(vw, system.Wind.vwn[self.wind])
        Pw = div(mul(ng, system.Dfig.windpower(rho, Vw, A, R, notu+omega_m, theta_p, 1)), mul(system.Settings.mva, 1e6))

        # mechanical torque
        Tm = mul(Pw, iomega)

        # motion equation
        system.DAE.f[self.omega_m] = mul(Tm-mul(xm, mul(iqr, ids)-mul(idr, iqs)), i2Hm)

        # speed control equation

        system.DAE.f[self.iqr] = mul(self.u, div(-div(mul((xs+xm), mul(pwa, iomega)), mul(V, xm))-iqr-self.dat7, Te))

        # voltage control equation
        system.DAE.f[self.idr] = mul(self.u, mul(Kv, (V-system.DAE.y[self.vref]))-div(V, xm)-idr)

        # pitch control equation
        contex = decimal.getcontext()
        contex.rounding = decimal.ROUND_05UP
        phi = matrix(0.0, (self.n, 1))

        for i in range(self.n):
            phi[i] = round(Decimal(1000 * (notu[i] + omega_m[i] - 1))) / 1000
        system.DAE.f[self.theta_p] = div(mul(self.u, mul(Kp, phi) - theta_p), Tp)

        # anti-windup limiter
        system.Dfig.windup(idx=self.iqr, xmax=self.dat8, xmin=self.dat9, type='f')
        system.Dfig.windup(idx=self.idr, xmax=self.dat10, xmin=self.dat11, type='f')
        system.Dfig.windup(idx=self.theta_p, xmax=[float('inf')]*self.n, xmin=[0]*self.n, type='f')

    def windup(self, idx, xmax, xmin, type):

        x = system.DAE.x[idx]
        if type == 'f':
            for i in range(len(xmax)):
                if (x[i] >= xmax[i]) and (system.DAE.f[idx[i]] > 0):
                    system.DAE.f[idx[i]] = 0
                    system.DAE.x[idx[i]] = xmax[i]
                elif (x[i] <= xmin[i]) and (system.DAE.f[idx[i]] < 0):
                    system.DAE.f[idx[i]] = 0
                    system.DAE.x[idx[i]] = xmin[i]

        if type == 'td':
            for i in range(len(xmax)):
                if (x[i] >= xmax[i] or x[i] <= xmin[i]) and (system.DAE.f[idx[i]] == 0):
                    system.DAE.tn[idx[i]] = 0
                    system.DAE.Ac[idx[i], :] = 0
                    system.DAE.Ac[:, idx[i]] = 0
                    system.DAE.Ac = system.DAE.Ac - spmatrix(1.0, [idx[i]], [idx[i]], ((system.DAE.nx + system.DAE.ny, system.DAE.nx + system.DAE.ny)))

    def windup_final(self):
        if self.n == 0:
            return
        system.Dfig.windup(idx=self.iqr, xmax=self.dat8, xmin=self.dat9, type='td')
        system.Dfig.windup(idx=self.idr, xmax=self.dat10, xmin=self.dat11, type='td')
        system.Dfig.windup(idx=self.theta_p, xmax=[float('inf')] * self.n, xmin=[0] * self.n, type='td')

    def Fxcall(self):

        if self.n == 0:
            return

        omega_m = system.DAE.x[self.omega_m]
        theta_p = system.DAE.x[self.theta_p]
        idr = system.DAE.x[self.idr]
        iqr = system.DAE.x[self.iqr]

        vw = system.DAE.x[system.Wind.vw[self.wind]]
        rho = system.Wind.rho[self.wind]

        pwa = system.DAE.y[self.pwa]
        V = system.DAE.y[self.v]
        t = system.DAE.y[self.a]
        st = sin(t)
        ct = cos(t)

        rs = self.rs
        xs = self.xs
        rr = self.rr
        xr = self.xr
        xm = self.xm

        i2Hm = mul(self.u, self.dat3)
        Kp = self.Kp
        Tp = self.Tp
        Kv = self.Kv
        Te = self.Te
        R = self.dat4
        A = self.dat5
        ng = self.ng

        as1 = rs ** 2 + self.dat1 ** 2
        a13 = div(rs, as1)
        a23 = div(self.dat1, as1)
        a33 = self.dat2

        vds = mul(-V, st)
        vqs = mul(V, ct)

        ids = mul(-a13, vds - mul(xm, iqr)) - mul(a23, vqs + mul(xm, idr))
        iqs = mul(a23, vds - mul(xm, iqr)) - mul(a13, vqs + mul(xm, idr))

        vdr = mul(-rr, idr) + mul(1 - omega_m, mul(a33, iqr) + mul(xm, iqs))
        vqr = mul(-rr, iqr) - mul(1 - omega_m, mul(a33, idr) + mul(xm, ids))

        notu = matrix(0, (self.n, 1))
        for i in range(self.n):
            if self.u[i] == 1:
                notu[i] = 0
            else:
                notu[i] = 1

        iomega = 1 / (omega_m + notu)
        Vwrate = system.Wind.vwn[self.wind]
        Vw = mul(vw, Vwrate)

        dPwdx = div(system.Dfig.windpower(rho, Vw, A, R, notu + omega_m, theta_p, 2), mul(system.Settings.mva, 1e6))
        Pw = div(mul(ng, system.Dfig.windpower(rho, Vw, A, R, notu + omega_m, theta_p, 1)),
                 mul(system.Settings.mva, 1e6))

        # mechanical torque
        Tm = mul(Pw, iomega)

        iqrsign = [1.0] * self.n
        w21 = 2 * omega_m - 1
        if len(w21) == 1:
            if (w21[0] <= 0) or (w21[0] >= 1):
                iqrsign = -1
        else:
            for i in range(len(w21)):
                if (w21[i] <= 0) or (w21[i] >= 1):
                    iqrsign[i] = -1
        w211 = matrix(0.0, (len(w21), 1))
        for i in range(len(w21)):
            w211[i] = min(w21[i], 1)
        Tsp1 = div(mul(self.Sn, mul(w211, iomega)), system.Settings.mva)
        Tsp = matrix(0.0, (len(Tsp1), 1))
        for i in range(len(Tsp1)):
            Tsp[i] = max(Tsp1[i], 0)


        slip = 1 - omega_m
        iqr_min = div(-self.Sn, system.Settings.mva)

        # df / dy
        idsv = mul(a13, st) - mul(a23, ct)
        idst = mul(a13, vqs) - mul(a23, vds)
        iqsv = mul(-a23, st) - mul(a13, ct)
        iqst = mul(-a23, vqs) - mul(a13, vds)

        ot = mul(mul(xm, i2Hm), mul(idr, iqst)-mul(iqr, idst))
        ov = mul(mul(xm, i2Hm), mul(idr, iqsv) - mul(iqr, idsv))
        iqrv = div(mul(xs+xm, Tsp), mul(Te, mul(V, mul(V, xm))))

        system.DAE.Fy = system.DAE.Fy \
                        + spmatrix(ot, self.omega_m, self.a, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(ov, self.omega_m, self.v, (system.DAE.nx, system.DAE.ny))

        # dg / dx
        idsidr = mul(-a23, xm)
        idsiqr = mul(a13, xm)
        iqsidr = mul(-a13, xm)
        iqsiqr = mul(-a23, xm)

        vdridr = -rr + mul(slip, mul(xm, iqsidr))
        vdriqr = mul(slip, (a33 + mul(xm, iqsiqr)))
        vqriqr = -rr - mul(slip, mul(xm, idsiqr))
        vqridr = mul(-slip, (a33 + mul(xm, idsidr)))

        vdrom = -(mul(a33, iqr) + mul(xm, iqs))
        vqrom = mul(a33, idr) + mul(xm, ids)

        pidr = mul(vds, idsidr) + mul(vqs, iqsidr) + mul(idr, vdridr) + vdr + mul(iqr, vqridr)
        piqr = mul(vds, idsiqr) + mul(vqs, iqsiqr) + mul(idr, vdriqr) + vqr + mul(iqr, vqriqr)
        pom = mul(idr, vdrom) + mul(iqr, vqrom)
        qidr = div(mul(-xm, V), self.dat1)

        system.DAE.Gx = system.DAE.Gx - spmatrix(mul(self.u, pom), self.a, self.omega_m, (system.DAE.ny, system.DAE.nx))

        # df / dx

        oidr = mul(mul(xm, i2Hm), mul(idr, iqsidr) + iqs - mul(iqr, idsidr))
        oiqr = mul(mul(xm, i2Hm), mul(idr, iqsiqr) - ids - mul(iqr, idsiqr))

        # mechanical equation
        dPwdx1 = matrix(0.0, (self.n, 1))
        dPwdx2 = matrix(0.0, (self.n, 1))
        for i in range(self.n):
            dPwdx1[i] = dPwdx[i]
            dPwdx2[i] = dPwdx[self.n*1+i]
        omedome = -notu + mul(mul(ng, dPwdx1)-Tm, mul(i2Hm, iomega))
        omedwind = mul(Vwrate, mul(ng, mul(dPwdx2), mul(i2Hm, iomega)))
        system.DAE.Fx = system.DAE.Fx + spmatrix(omedome, self.omega_m, self.omega_m, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(omedwind, self.omega_m, system.Wind.vw[self.wind], (system.DAE.nx, system.DAE.nx))

        # pitch angle control equation
        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if (theta_p[i] > 0) and (self.u[i]):
                z[i] = 1

        dPwdx3 = matrix(0.0, (self.n, 1))
        for i in range(self.n):
            dPwdx3[i] = dPwdx[self.n*2+i]

        ones = matrix(1, (self.n, 1))

        omedtheta = mul(mul(z, ng), mul(dPwdx3, mul(i2Hm, iomega)))
        system.DAE.Fx = system.DAE.Fx - spmatrix(div(ones, Tp), self.theta_p, self.theta_p, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(div(mul(z, Kp), Tp), self.theta_p, self.omega_m, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(omedtheta, self.omega_m, self.theta_p, (system.DAE.nx, system.DAE.nx))

        # speed control equation

        kiqr = div(mul(-matrix(iqrsign), (xs + xm)), mul(V, mul(xm, Te)))
        tspo = 2 * div(self.Sn, system.Settings.mva)

        iqrp = div(mul(-(xs + xm), iomega), mul(V, mul(xm, Te)))
        iqrw = div(mul((xs + xm), mul(pwa, mul(iomega, iomega))), mul(V, mul(xm, Te)))

        for i in range(len(Tsp)):
            if (Tsp[i] == 0) & (w21[i] >= 1):
                tspo[i] = 0
        iqr_max = self.dat8
        iqr_min = self.dat9

        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if (iqr[i] < iqr_max[i]) and (iqr[i] > iqr_min[i]):
                if self.u[i]:
                    z[i] = 1

        system.DAE.Fx = system.DAE.Fx - spmatrix(div(ones, Te), self.iqr, self.iqr, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(z, iqrw), self.iqr, self.omega_m, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(z, oiqr), self.omega_m, self.iqr, (system.DAE.nx, system.DAE.nx))

        system.DAE.Fy = system.DAE.Fy + spmatrix(mul(z, iqrv), self.iqr, self.v, (system.DAE.nx, system.DAE.ny)) \
                        + spmatrix(mul(z, iqrp), self.iqr, self.pwa, (system.DAE.nx, system.DAE.ny))

        system.DAE.Gx = system.DAE.Gx - spmatrix(mul(z, piqr), self.a, self.iqr, (system.DAE.ny, system.DAE.nx))

        # voltage control equation
        idr_max = self.dat10
        idr_min = self.dat11
        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if (idr[i] < idr_max[i]) and (idr[i] > idr_min[i]):
                if self.u[i]:
                    z[i] = 1

        system.DAE.Fx = system.DAE.Fx + spmatrix([-1]*self.n, self.idr, self.idr, (system.DAE.nx, system.DAE.nx)) \
                        + spmatrix(mul(z, oidr), self.omega_m, self.idr, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy + spmatrix(mul(z, Kv-1/xm), self.idr, self.v, (system.DAE.nx, system.DAE.ny)) \
                        - spmatrix(mul(z, Kv), self.idr, self.vref, (system.DAE.nx, system.DAE.ny))
        system.DAE.Gx = system.DAE.Gx - spmatrix(mul(z, pidr), self.a, self.idr, (system.DAE.ny, system.DAE.nx)) \
                        - spmatrix(mul(z, qidr), self.v, self.idr, (system.DAE.ny, system.DAE.nx))

        # power reference equation
        z = matrix(0, (self.n, 1))
        for i in range(self.n):
            if (omega_m[i] < 1) and (omega_m[i] > 0.5):
                z[i] = 1
        system.DAE.Gx = system.DAE.Gx + spmatrix(2*mul(z, self.Sn)/system.Settings.mva, self.pwa, self.omega_m, (system.DAE.ny, system.DAE.nx))

    def suiji(self, sigmavw, suijih, vw, syst=-1.0, type='1'):

        k = 0
        t = -1.0
        if type == '1':
            while True:
                t1 = time.time()
                try:
                    # tnew = syst.get()
                    tnew = syst.value
                except EOFError:
                    break
                if round(tnew, 2) != round(t, 2):
                    k = k + 1
                    t = tnew
                    # print('风速随机波动模拟模块')

                    if round(round(t, 2) % suijih, 1) == 0.0:  # 是否有更好的表达形式

                        for i in range(self.n):
                            a = random.normalvariate(0, sigmavw)
                            vw[i] = vw[i] + a

                        # j = int(round(t, 2) // suijih)
                        #
                        # for i in range(self.n):
                        #     a = random.normalvariate(0, sigmavw)
                        #     b = list(vw[i])
                        #     b[j+1] = b[j+1] + a
                        #     vw[i] = numpy.array(b)
                    # time.sleep(0.05)

                    # t2 = time.time()
                    # print('风速随机波动模拟耗时：%s' % (t2-t1))
                    # print('风速随机波动模块：第%s次，仿真时间%ss' % (k, t))
        if type == '2':

            t1 = time.time()
            t = syst
            # t1 = time.time()
            if round(round(t, 2) % suijih, 1) == 0.0:  # 是否有更好的表达形式

                # j = int(round(t, 2) // suijih)

                for i in range(self.n):
                    a = random.normalvariate(0, sigmavw)
                    vw[i] = vw[i] + a
                system.DAE.x[system.Wind.vw] = vw
                #     b = list(ssvw[i])
                #     b[j + 1] = b[j + 1] + a
                #     ssvw[i] = numpy.array(b)
                # system.Wind.svw = ssvw
            # print('风速随机波动模块：仿真时间%ss' % t)
            t2 = time.time()
            # print('风速随机波动模拟耗时：%s' % (t2-t1))
