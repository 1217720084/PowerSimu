import system
from devices.base_device import base_device
from cvxopt.base import mul, matrix, exp, div, spmatrix
import random

class stoFault():
    def __init__(self):
        self._type = 'stochastic'
        self._name = 'stochastic'
        self.tf = []  # 输电线路故障发生时间
        self.tc = []  # 输电线路故障切除时间


    def setup(self):
        self.tf = [-1 for i in range(system.Line.n)]
        self.tc = [-1 for i in range(system.Line.n)]
    def istime(self, t):

        return

    def settimes(self, linep, t):

        [MWs, MWr] = system.Line.flow()
        pmin = 0.008
        pmax = 0.999
        ld = 0.95
        lu = 1.3
        # 系数
        alpha = 100
        pj = []
        p = MWs.real()
        p = abs(p)

        for i in range(system.Line.n):
            if linep[i] == 0:
                # print('第%i条线路传输功率为0' % i)
                pj.append(-1)
            else:
                if p[i] / linep[i] / alpha < ld:
                    pj.append(pmin)
                else:
                    if p[i] / linep[i] / alpha > lu:
                        pj.append(pmax)
                    else:
                        pj.append((pmax - pmin) / (lu - ld) * (p[i] / linep[i] / alpha - ld) + pmin)

        for i in range(system.Line.n):
            a = random.random()
            if a < pj[i]:
                if self.tc[i] <= t:
                    self.tf[i] = t + 0.1
                    self.tc[i] = t + 0.11

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