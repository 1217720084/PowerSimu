"""

"""
class settings():
    def __init__(self):
        self.pv2pq = 0   # PV是否能转为PQ负荷
        self.iter = 20   # 最后一次潮流的迭代次数，默认20
        self.pv2pqiter = 0  # 经过多少次迭代后PV才能转为PQ
        self.error = 1e-6   # 默认误差
        self.multipvswitch = 0  # 允许多台发电机转为PQ负荷
        self.freq = 50          # 频率
        self.mva = 100
        self.t0 = 0.0  # 初始仿真时间，静态分析时为-1
        self.tf = 10  # 最终仿真时间
        self.deltatmax = 0.125  # 最大时间步长
        self.deltatmin = 0.00039063   # 最小时间步长
        self.chunk = 100  # 输出数组的初始维度
        self.deltat = 6.2178e-17  # 时域仿真的时间步长
        self.fixt = 0  # 设置固定时间步长，默认为0
        self.tstep = 0.05  # 固定的时间步长值
        self.resetangles = 1  # 清除故障后，将节点电压置为故障前的值
        self.method = 2  # 时域仿真才有欧拉法（值为1）或隐式梯形法
        self.dynmit = 20
        self.dyntol = 0.00001
        self.deltadelta = 180
        self.disturbance = False
        self.forcepq = False
        self.pq2z = 0  # 是否把pq转为恒阻抗模型，默认为0


