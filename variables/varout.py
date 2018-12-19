"""

"""
import system
from cvxopt.base import matrix
class varout():
    def __init__(self):
        self.t = []  # 时域仿真时间
        self.idx = []  # 时域仿真索引
        self.var = []  # 输出变量

    def store(self, t, k):

        self.t.append(t)
        self.idx.append(k-1)
        var = []
        if system.DAE.nx:
            var = list(system.DAE.x)
        if system.DAE.ny:
            var.extend(list(system.DAE.y))

        # 增加存储PQ的Pl、Ql
        var.extend(list(system.PQ.Pl))
        var.extend(list(system.PQ.Ql))

        # 增加存储每个时刻Line的开断时间
        var.extend(system.Line.tf)
        var.extend(system.Line.tc)

        # 增加存储每个时刻vw改变后的值
        # var.extend(list(system.Wind.vwt))

        # 增加存储每个时刻svw的值
        # svw = list(system.Wind.svw)
        # if round(round(t, 2) % 0.01, 1) == 0.0:
        #     j = int(round(t, 2) // 0.01)
        #     for i in range(system.Wind.n):
        #         s = list(svw[i])
        #         var.append(s[j])

        self.var.append([])
        self.var[k-1] = var

        return

    def clear(self):

        self.t = []  # 时域仿真时间
        self.idx = []  # 时域仿真索引
        self.var = []  # 输出变量