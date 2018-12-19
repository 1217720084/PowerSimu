"""

"""
from cvxopt.base import matrix

class dae():
    def __init__(self):
        self.g = []          # 代数方程
        self.f = []          # 微分方程
        self.y = []          # 代数变量
        self.x = []          # 状态变量
        self.Y = matrix()          # 导纳矩阵
        self.Y_G = matrix()        # 导纳矩阵实部
        self.Y_B = matrix()        # 导纳矩阵虚部
        self.Gy = []         # 代数方程的雅可比矩阵
        self.Fx = []
        self.Fy = []
        self.Gx = []
        self.nx = 0          # 状态变量个数
        self.ny = 0          # 代数变量个数，默认为0
        self.n_bus = 0        # 母线节点数，默认为0
        self.t = []



        self.factorize =[]  # 如果是True，则对雅可比矩阵进行因式分解
        self._params = ['g', 'x', 'y', 'f', 'Fx', 'Fy', 'Gx', 'Y', 'Y_G', 'Y_B', 'Gy']

        self.factorize = True  # 如果是True，则对雅可比矩阵进行因式分解
        self.tn = []  # 时域仿真的微分方程
        self.Ac = []  # 时域仿真雅可比方程



    def _list2matrix(self):

        for item in self._params:

            self.__dict__[item] = matrix(self.__dict__[item])