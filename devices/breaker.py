"""

"""
import system
from devices.base_device import base_device
class breaker(base_device):
    def __init__(self):

        base_device.__init__(self)
        self._data.update(
            {'bus': 0, 'line': 0, 'Sn': 100, 'Vn': 220, 'fn': 50, 'u': 0, 't1': 0, 't2': 0,'u1': 1, 'u2': 1})
        self._type = 'breaker'
        self._name = 'breaker'
        self.n = 0

    def gettimes(self):

        t = []
        if self.n == 0:
            return t

        # u = np.unique(self.tf,self.tc)    #此处省略同时发生故障的情况，待完善
        # print(self.tf[0])

        t1 = self.t1[0] - 0.000001
        t2 = self.t2[0] - 1e-6

        a = [t1, self.t1[0]]
        # print(a)
        b = [t2, self.t2[0]]
        a.extend(b)
        # print(a)
        t = a

        return t
        # t = [u,u]
        # print(t)

    def intervention(self, t):
        #
        action = ['Opening', 'Closing']
        for item in range(self.n):
            if t == self.t1[item]:  # 发生故障
                if self.u[item] == 0:
                    self.u[item] = 1
                else:
                    self.u[item] = 0
                print('%s break at bus %s on line from %s to %s for t = %s s' % (action[self.u[item]], self.bus[item], system.Line.af[self.line[item]-1]+1, system.Line.at[self.line[item]-1]+1, self.t1[item]))

                #启用故障
                system.Line.setstatus(self.line[item], self.u[item])
                # multiagent


            if t == self.t2[item]:  # 发生故障

                if self.u[item] == 0:
                    self.u[item] = 1
                else:
                    self.u[item] = 0
                print(action[self.u] + 'break at bus %s on line from %s to %s for t = %s' % (
                self.bus[item], system.Line.af[self.line[item] - 1], system.Line.at[self.line[item] - 1],
                self.t1[item]))

                # 启用故障
                system.Line.setstatus(self.line[item], self.u[item])

        # 检查网络连通性