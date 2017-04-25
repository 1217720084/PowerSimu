"""

"""
from cvxopt.base import matrix

class dae():
    def __init__(self):
        self.g = []          # 代数方程
        self.y = []          # 代数变量
        self.Y = matrix()          # 导纳矩阵
        self.Y_G = matrix()        # 导纳矩阵实部
        self.Y_B = matrix()        # 导纳矩阵虚部
        self.Gy = []         # 代数方程的雅可比矩阵
        self.nx = 0          # 状态变量个数
        self.ny = 0          # 代数变量个数，默认为0
        self.n_bus = 0        # 母线节点数，默认为0
        self.factorize =[]  # 如果是True，则对雅可比矩阵进行因式分解
        self._params = ['g', 'y', 'Y', 'Y_G', 'Y_B', 'Gy']

    def _list2matrix(self):

        for item in self._params:

            self.__dict__[item] = matrix(self.__dict__[item])