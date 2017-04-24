"""

"""
from devices.base_device import base_device
import system


class pq(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update({'Pl': 1, 'bus': None, 'Ql': 1,  'Vmax': 1.1, 'Vmin': 0.95})
        self._type = 'PQ'
        self._name = 'PQ'
        self._bus = {'bus': ['a', 'v']}
        self._algebs = ['V0', 'Va']
        self._params.extend(['Pl', 'Ql',  'Vmax', 'Vmin'])
        self._powers = ['Pl', 'Ql']

    def yinit(self, dae):

        dae.y = system.DAE.y
        dae.g = system.DAE.g
        # dae.Gy = sparse(m*m)
        for key, value in zip(self.v, self.Ql):
            dae.g[key] += value
        for key, value in zip(self.a, self.Pl):
            dae.g[key] += value

    def gcall(self):
        for i in range(system.PQ.n):
            system.DAE.g[self.a[i]] += self.Pl[i]
            system.DAE.g[self.v[i]] += self.Ql[i]
        # 判断是否PQ负荷转为恒阻抗模型

        a = []
        b = []
        for i in range(len(self.n)):
            if system.DAE.y[self.v[i]] < self.Vmin[i]:
                system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] ** system.DAE.y[self.v[i]] / \
                                                                         self.Vmin[i] / self.Vmin[i]
                system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] ** system.DAE.y[self.v[i]] / \
                                                                         self.Vmin[i] / self.Vmin[i]
            if system.DAE.y[self.v[i]] > self.Vmax[i]:
                system.DAE.g[self.a[i]] = system.DAE.g[i] - self.Pl[i] + self.Pl[i] ** system.DAE.y[self.v[i]] / \
                                                                         self.Vmin[i] / self.Vmin[i]
                system.DAE.g[self.v[i]] = system.DAE.g[i] - self.Ql[i] + self.Ql[i] ** system.DAE.y[self.v[i]] / \
                                                                         self.Vmin[i] / self.Vmin[i]









    def Gycall(self,dae):
        i = 0
        while i < dae.n_bus:           #总共n个节点
            if self.a != i:             #i不等于j的情况
                dae.Gy[self.a][i] = - dae.y[v] * dae.y[i] * (system.DAE.Y_G[self.a][i] * sin(dae.y[a] - dae.y[i]) - system.DAE.Y_B[self.a][i] * cos(dae.y[a] - dae.y[i]))     #H[i][j] 实质是某一行的元素，在主程序调用时还需要循环一次
                dae.Gy[self.v][i] = dae.y[v] * dae.y[i] * (system.DAE.Y_G[self.a][i] * cos(dae.y[a] - dae.y[i]) + system.DAE.Y_B[self.a][i] * sin(dae.y[a] - dae.y[i]))        #J[i][j]
                dae.Gy[self.s][i] = - dae.y[v] * dae.y[i] * (system.DAE.Y_G[self.a][i] * cos(dae.y[a] - dae.y[i]) + system.DAE.Y_B[self.a][i] * sin(dae.y[a] - dae.y[i]))      #N[i][j]
                dae.Gy[self.t][i] = - dae.y[v] * dae.y[i] * (system.DAE.Y_G[self.a][i] * sin(dae.y[a] - dae.y[i]) - system.DAE.Y_B[self.a][i] * cos(dae.y[a] - dae.y[i]))      #L[i][j]
            else:
                dae.Gy[self.a][self.a] = dae.y[v] * dae.y[v] * system.DAE.Y_B[self.a][self.a] + self.Ql + self.Qg   #当j不等于i时的表达式，self.Ql + self.Qg表示Qi可能有问题
                dae.Gy[self.v][self.a] = dae.y[v] * dae.y[v]  * system.DAE.Y_G[self.a][self.a] - (self.Pl + self.Qg)
                dae.Gy[self.s][self.a] = - dae.y[v] * dae.y[v] * system.DAE.Y_G[self.a][self.a] - (self.Pl + self.Qg)
                dae.Gy[self.t][self.a] = dae.y[v] * dae.y[v] * system.DAE.Y_B[self.a][self.a] - (self.Ql + self.Qg)
            i += 1


