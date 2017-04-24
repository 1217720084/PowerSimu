"""

"""
class settings():
    def __init__(self):
        self.pv2pq = 0   # PV是否能转为PQ负荷
        self.iter = 20   # 最后一次潮流的迭代次数，默认20
        self.pv2pqiter = 0  # 经过多少次迭代后PV才能转为PQ
        self.error = 1e-6   # 默认误差
        self.multipvswitch = 0  # 允许多台发电机转为PQ负荷
