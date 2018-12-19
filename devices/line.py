# _*_ coding:utf-8 _*_
"""

"""
from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, cos, sin, sparse, mul, exp, div
import cvxopt.blas
import math
import random
import numpy as np
import time


class line(base_device):
    def __init__(self):
        base_device.__init__(self)
        self._data.update({'fn': 50, 'f_bus': None, 'to_bus': None, 'l': 0,  'kT': 1, 'r': 2, 'x': 1.1, 'b': 1, 'tap_ratio': 1.1, 'theta': 0.1,
                           'Imax': 0, 'Pmax': 0, 'Smax': 0, 'tf': -1, 'tc': -1})
        self._type = 'Line'
        self._name = 'Line'
        self.n = 0
        self._bus = {'f_bus': ['af', 'vf'], 'to_bus': ['at', 'vt']}
        self._params.extend(['fn', 'kT', 'r', 'x', 'tap_ratio', 'theta', 'Imax', 'Pmax', 'Smax'])
        self.z = ['r', 'x']
        self.properties.update({'gcall': True, 'Gycall': True})

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

        if self.n == 0:
            return
        for i in range(len(self.tap_ratio)):
            if self.tap_ratio[i] == 0:
                self.tap_ratio[i] = 1
        nb = system.Bus.n

        # 处理传输线数据和生成导纳矩阵

        chrg = mul(matrix(self.u), 0.5 * matrix(self.b))
        y = div(matrix(self.u) + 0j, matrix(matrix(self.r) + 1j * matrix(self.x)))
        ts = exp(matrix(self.theta) * math.pi / 180 * 1j)
        ts = mul(matrix(self.tap_ratio), ts)
        ts2 = mul(ts, ts.H.T)

        # ts2real = ts2.real()
        # ts2imag = ts2.imag()
        # for i in range(len(ts2)):
        #     ts2real[i] = round(ts2real[i], 15)
        #     ts2imag[i] = round(ts2imag[i], 15)
        # ts2 = ts2real + 1j*ts2imag

        yft = div(-y, ts.H.T)
        # yftreal = yft.real()
        # yftimag = yft.imag()

        ytf = div(-y, ts)
        # ytfreal = ytf.real()
        # ytfimag = ytf.imag()

        yff = div(y, ts2) + 1j*chrg
        # yffreal = yff.real()
        # yffimag = yff.imag()

        ytt = y + 1j*chrg
        # yttreal = ytt.real()
        # yttimag = ytt.imag()

        # for i in range(len(ytt)):
        #     yftreal[i] = round(yftreal[i], 15)
        #     yftimag[i] = round(yftimag[i], 15)
        #
        #     ytfreal[i] = round(ytfreal[i], 15)
        #     ytfimag[i] = round(ytfimag[i], 15)
        #
        #     yffreal[i] = round(yffreal[i], 15)
        #     yffimag[i] = round(yffimag[i], 15)
        #
        #     yttreal[i] = round(yttreal[i], 15)
        #     yttimag[i] = round(yttimag[i], 15)
        #
        # yft = yftreal + 1j * yftimag
        # ytt = yttreal + 1j * yttimag
        # yff = yffreal + 1j * yffimag
        # ytf = ytfreal + 1j * ytfimag



        self.Y = spmatrix(yft, self.af, self.at, (nb, nb)) +\
            spmatrix(ytf, self.at, self.af, (nb, nb)) +\
            spmatrix(yff, self.af, self.af, (nb, nb)) +\
            spmatrix(ytt, self.at, self.at, (nb, nb))

        for i in range(nb):
            if self.Y[i, i] == 0:
                self.Y[i, i] == -1j*1e-6
        system.DAE.Y = self.Y

        system.DAE.Y_G = self.Y.real()

        system.DAE.Y_B = self.Y.imag()



        # for i in range(len(self.tap_ratio)):
        #     if self.tap_ratio[i] == 0:
        #         self.tap_ratio[i] = 1
        # # print(self.tap_ratio)
        #
        # # chrg = np.mat(self.b) * 0.5
        # # chrg = chrg.T
        # #
        # # print(chrg)
        # chrg = np.array(np.zeros((len(self.b), 1), complex))
        # for i in range(len(self.x)):
        #     chrg[i] = complex(0.0, self.b[i] * 0.5)
        #
        # #print(chrg)
        # #zeros = [0] * len(self.x)
        # #y = matrix(zeros, (len(self.x), 1, complex))
        # y = np.array(np.zeros((len(self.x), 1), complex))
        # for i in range(len(self.x)):
        #     y[i] = 1.0 / complex(self.r[i], self.x[i])
        # #print(y)
        #
        # ts = np.array(np.zeros((len(self.theta), 1), complex))
        # for i in range(len(self.theta)):
        #     ts[i] = complex(self.tap_ratio[i]*cos(self.theta[i]*math.pi/180), self.tap_ratio[i]*sin(self.theta[i]*math.pi/180))
        # #print(ts)
        #
        # ts2 = ts * ts.conj()
        # #print(ts2)
        #
        # y1 = -y / ts.conj()
        # #print(y1)
        #
        # self.Y = spmatrix(y1, self.af, self.at, (system.Bus.n, system.Bus.n)) +\
        #     spmatrix(y1, self.at, self.af, (system.Bus.n, system.Bus.n)) +\
        #     spmatrix(y/ts2 + chrg, self.af, self.af, (system.Bus.n, system.Bus.n)) +\
        #     spmatrix(y + chrg, self.at, self.at, (system.Bus.n, system.Bus.n))
        # system.DAE.Y = self.Y
        #
        # system.DAE.Y_G = self.Y.real()
        #
        # system.DAE.Y_B = self.Y.imag()





    def gcall(self):

        system.DAE.y = matrix(system.DAE.y)
        zeros = [0.0] * system.DAE.ny
        system.DAE.g = zeros
        system.DAE.g = matrix(system.DAE.g)
        Vn = exp(system.DAE.y[system.Bus.a] * 1j)
        Vc = mul(system.DAE.y[system.Bus.v] + 0j, Vn)
        Ic = self.Y * Vc
        S = mul(Vc, Ic.H.T)

        self.p = S.real()

        self.q = S.imag()

        for i in range(system.Bus.n):
             system.DAE.g[i] = self.p[i]
             system.DAE.g[i+system.Bus.n] = self.q[i]




    def Gycall(self):

        system.DAE.y = matrix(system.DAE.y)
        U = exp(system.DAE.y[system.Bus.a]*1j)
        V = mul(system.DAE.y[system.Bus.v] + 0j, U)
        I = self.Y * V
        nb = len(system.Bus.a)
        diagU = spmatrix(U, system.Bus.a, system.Bus.a, (nb, nb), 'z')
        diagV = spmatrix(V, system.Bus.a, system.Bus.a, (nb, nb), 'z')
        diagI = spmatrix(I, system.Bus.a, system.Bus.a, (nb, nb), 'z')
        dS = self.Y*diagU

        dS = diagV*dS.H.T


        dS += diagI.H.T * diagU
        dR = diagI
        dR = dR-self.Y*diagV

        dR = diagV.H.T*dR

        # system.DAE.Gy = matrix(0.0, (system.DAE.ny, system.DAE.ny))
        #
        # system.DAE._list2matrix()
        # system.DAE.Gy = spmatrix(dR.imag().V, dR.imag().I, dR.imag().J, (system.DAE.ny, system.DAE.ny)) \
        #                 + spmatrix(dR.real().V, dR.real().I, dR.real().J+system.Bus.n, (system.DAE.ny, system.DAE.ny)) \
        #                 + spmatrix(dS.real().V, dS.real().I+system.Bus.n, dS.real().J, (system.DAE.ny, system.DAE.ny)) \
        #                 + spmatrix(dS.imag().V, dS.imag().I+system.Bus.n, dS.imag().J+system.Bus.n, (system.DAE.ny, system.DAE.ny))


        Gy = sparse([[dR.imag(),dR.real()],[dS.real(),dS.imag()]])
        system.DAE.Gy = spmatrix(Gy.V, Gy.I, Gy.J, (system.DAE.ny, system.DAE.ny))

        system.DAE.Gy = matrix(system.DAE.Gy)
    def base(self):

        # 找到变压器（kT!=0）
        idx = []
        kT = []
        for i in range(self.n):
            if self.kT[i] != 0:
                idx.append(i)
                kT.append(self.kT[i])
        VL1 = self.Vn

        # 变压器长度设为0

        V1 = []
        V2 = []
        KT = []
        for i in idx:
            self.l[i] = 0
        for i in range(self.n):

            V1.append(system.Bus.Vb[self.af[i]])
            V2.append(system.Bus.Vb[self.at[i]])

        for i in idx:
            KT.append(V1[i]/V2[i])

        for i in range(len(KT)):
            KT[i] = round(KT[i], 15)


        kT = matrix(kT)
        KT = matrix(KT)
        idx1 = []
        if len(idx) != 0:
            corr = div(abs(kT-KT), KT)
            for i in range(len(corr)):
                if corr[i] > 0.1:
                    idx1.append(i)
            # 显示变压器调压比误差超过连个节点的10%
            for i in range(len(idx1)):
                print('显示变压器调压比误差超过两个节点电压比的10%')
        # 调整变压器调压比如果不符合电压基值

            for i in range(len(idx)):
                if self.tap_ratio[idx[i]] == 0.0:
                    self.tap_ratio[idx[i]] = 1.0

            for i in range(len(idx)):
                self.tap_ratio[idx[i]] = self.tap_ratio[idx[i]] * kT[i]/KT[i]

        idx2 = []
        for i in range(self.n):
            if abs(VL1[i]-V1[i])/V1[i] > 0.1:
                idx2.append(i)
        for i in range(len(idx2)):
            print('显示输电线路电压误差超过连接节点的10%')

        V1 = matrix(V1)
        VL1 = matrix(VL1)
        Vb2new = mul(V1, V1)
        Vb2old = mul(VL1, VL1)

        # 化标幺值

        for i in range(self.n):
            self.r[i] = Vb2old[i] * self.r[i] / self.Sn[i] / Vb2new[i] * system.Settings.mva
            self.x[i] = Vb2old[i] * self.x[i] / self.Sn[i] / Vb2new[i] * system.Settings.mva
            self.b[i] = Vb2new[i] * self.b[i] * self.Sn[i] / Vb2old[i] / system.Settings.mva

            self.Imax[i] = self.Imax[i] * self.Sn[i] * V1[i] / VL1[i] / system.Settings.mva
            self.Pmax[i] = self.Pmax[i] * self.Sn[i] / system.Settings.mva
            self.Smax[i] = self.Smax[i] * self.Sn[i] / system.Settings.mva

    def flow(self, y):

        if y is not None:

            system.DAE.y = matrix(y)
        if self.n == 0:
            return

        chrg = mul(matrix(self.u), 0.5 * matrix(self.b))
        y = div(matrix(self.u) + 0j, matrix(matrix(self.r) + 1j * matrix(self.x)))
        ts = exp(matrix(self.theta) * math.pi / 180 * 1j)
        tps = mul(matrix(self.tap_ratio), ts)
        tpj = tps.H.T
        ts2 = mul(tps, tpj)

        Vf = mul(system.DAE.y[self.vf], exp(system.DAE.y[self.af]*1j))
        Vt = mul(system.DAE.y[self.vt], exp(system.DAE.y[self.at] * 1j))

        I1 = div(mul(Vf, y+1j*chrg),ts2)-div(mul(Vt, y), tpj)
        I2 = mul(Vt, y + 1j * chrg) - div(mul(Vf, y), tps)
        MWs = mul(Vf, I1.H.T)
        MWr = mul(Vt, I2.H.T)

        return [MWs, MWr]

    # def guzhang(self, linep):
    #
    #     [MWs, MWr] = self.flow()
    #     pmin = 0.008
    #     pmax = 0.999
    #     ld = 0.95
    #     lu = 1.3
    #     # 系数
    #     alpha = 1.7
    #     pj = []
    #     p = MWs.real()
    #
    #     for i in range(self.n):
    #         if linep[i] == 0:
    #             # print('第%i条线路传输功率为0' % (i+1))
    #             pj.append(-1)
    #         else:
    #             if p[i]/linep[i]/alpha < ld:
    #                 pj.append(pmin)
    #             else:
    #                  if p[i]/linep[i]/alpha > lu:
    #                      pj.append(pmax)
    #                  else:
    #                      pj.append((pmax-pmin)/(lu-ld)*(p[i]/linep[i]/alpha-ld)+pmin)
    #
    #     for i in range(self.n):
    #         a = random.random()
    #         if a < pj[i]:
    #
    #             self.u[i] = 0
    #             print('第%i条线路断开' % (i+1))

    # 设置线路最大传输容量

    def SetSmax(self):


        # 找到输电线路索引
        idx = []

        for i in range(self.n):
            if self.kT[i] == 0:
                idx.append(i)
        # 根据输电线路电压等级设置最大传输容量
        for i in range(len(idx)):
            if self.Smax[idx[i]] == 0.0:
                if self.Vn[idx[i]] > 765.0:
                    self.Smax[idx[i]] = 25
                elif self.Vn[idx[i]] > 500.0:
                    self.Smax[idx[i]] = 15
                elif self.Vn[idx[i]] > 330.0:
                    self.Smax[idx[i]] = 8
                elif self.Vn[idx[i]] > 220.0:
                    self.Smax[idx[i]] = 5
                elif self.Vn[idx[i]] > 110.0:
                    # self.Smax[idx[i]] = 0.5
                    self.Smax[idx[i]] = 5
                elif self.Vn[idx[i]] > 69.0:
                    # self.Smax[idx[i]] = 0.3
                    self.Smax[idx[i]] = 3
                elif self.Vn[idx[i]] > 35.0:
                    self.Smax[idx[i]] = 2.5
                elif self.Vn[idx[i]] > 13.8:
                    self.Smax[idx[i]] = 2
                else:
                    self.Smax[idx[i]] = 1


    def settimes(self, t):

        t1 = time.time()
        # 找到输电线路
        idx = []
        for i in range(self.n):
            if self.kT[i] == 0:
                idx.append(i)
        # 求出当前时刻线路传输容量
        [MWs, MWr] = self.flow(system.DAE.y)
        MWs = abs(MWs)

        # 设置故障概率模型
        pj = []
        for i in range(len(idx)):
            if MWs[idx[i]]/self.Smax[idx[i]] > 0.9:
                pj.append(0.0003)
            elif MWs[idx[i]]/self.Smax[idx[i]] > 0.7:
                pj.append(0.0002)
            elif MWs[idx[i]] / self.Smax[idx[i]] > 0.5:
                pj.append(0.0001)
            elif MWs[idx[i]] / self.Smax[idx[i]] > 0.3:
                pj.append(0.00005)
            else:
                pj.append(0)
        # 设置故障时间
        for i in range(len(idx)):
            a = random.random()
            if a < pj[i]:

                if self.tc[idx[i]] < t:  # 当故障清除时间小于当前仿真时刻时，故障已处理，重新设置故障时间
                    self.tf[idx[i]] = t + 0.01+1e-6
                    self.tc[idx[i]] = t + 0.07
        # print('输电线路随机故障模块：仿真时间%ss' % t)
        t2 = time.time()
        # print('输电线路随机故障模拟耗时：%s' % (t2-t1))

    def setfault(self, tf, tc, suijih, Sline, syst=-1.0):

        k = 0
        t = -1.0
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
                # print('输电线路随机模拟模块')

                if round(round(t, 2) % suijih, 1) == 0.0:  # 是否有更好的表达形式

                    # 找到输电线路
                    idx = []
                    for i in range(self.n):
                        if self.kT[i] == 0:
                            idx.append(i)

                    # 求出当前时刻线路传输容量
                    # t3 = time.time()
                    # [MWs, MWr] = self.flow(y)
                    # MWs = abs(MWs)
                    # t4 = time.time()
                    # print('计算潮流耗时：%s' % (t4 - t3))
                    MWs = Sline

                    # 设置故障概率模型

                    a = matrix(MWs)
                    b = matrix(self.Smax)
                    c = div(a[idx], b[idx])
                    n = len(c)
                    pj = [0.0] * n
                    for i in range(n):
                        if c[i] < 0.3:
                            pj[i] = 0.0
                        elif c[i] < 0.5:
                            pj[i] = 0.00005
                        elif c[i] < 0.7:
                            pj[i] = 0.0001
                        elif c[i] < 0.9:
                            pj[i] = 0.0002
                        else:
                            pj[i] = 0.0003
                    # pj = []
                    # for i in range(len(idx)):
                    #     if MWs[idx[i]] / self.Smax[idx[i]] > 0.9:
                    #         pj.append(0.003)
                    #     elif MWs[idx[i]] / self.Smax[idx[i]] > 0.7:
                    #         pj.append(0.002)
                    #     elif MWs[idx[i]] / self.Smax[idx[i]] > 0.5:
                    #         pj.append(0.001)
                    #     elif MWs[idx[i]] / self.Smax[idx[i]] > 0.3:
                    #         pj.append(0.0005)
                    #     else:
                    #         pj.append(0)


                    # 设置故障时间
                    for i in range(len(idx)):
                        a = random.random()
                        if a < pj[i]:
                            if tc[idx[i]] < t:  # 当故障清除时间小于当前仿真时刻时，故障已处理，重新设置故障时间
                                tf[idx[i]] = t + 0.01 + 1e-6
                                tc[idx[i]] = t + 0.06
                # time.sleep(0.05)
                # t2 = time.time()
                # print('输电线路随机故障模拟耗时：%s' % (t2-t1))
                # print('输电线路随机故障模块：第%s次，仿真时间%ss' % (k, t))

    def gettimes(self, fixed_times):

        for i in range(system.Line.n):
            if self.tf[i] > 0:
                tf = self.tf[i] - 1e-6
                tc = self.tc[i] - 1e-6
                a = [tf, self.tf[i]]
                # print(a)
                b = [tc, self.tc[i]]
                a.extend(b)
                fixed_times.extend(a)

        return fixed_times

    def istime(self, t):

        u = False
        for i in range(system.Line.n):
            if self.tf[i] == t or self.tc[i] == t:
                u = True
        return u

    def intervention(self, t):
        #
        for item in range(system.Line.n):
            if t == self.tf[item]:      #发生故障
                print('线路%i在t = %s s断开' % (item+1, self.tf[item]))

                #启用故障

                system.Line.u[item] = 0

                #存储故障前
                self._ang = matrix(system.DAE.y[system.Bus.a])
                self._vol = matrix(system.DAE.y[system.Bus.n:len(system.DAE.y)])

            if t == self.tc[item]:    #故障清除

                print('线路%i在t = %s s恢复连接' % (item+1, self.tc[item]))

                # 禁用故障
                system.Line.u[item] = 1


                #恢复电压
                if system.Settings.resetangles:
                    system.DAE.y[system.Bus.a] = matrix(self._ang)
                    system.DAE.y[system.Bus.n:len(system.DAE.y)] = matrix(self._vol)

        system.Line.build_y()

    def setstatus(self, idx, u):

        self.u[idx-1] = u
        system.Line.build_y()