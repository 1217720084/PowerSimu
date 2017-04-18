# _*_ coding:utf-8 _*_
"""

"""
class dae():
    def __init__(self):
        self.g = []          # 代数方程
        self.y = []          # 代数变量
        self.Y = []          # 导纳矩阵
        self.Gy = []         # 代数方程的雅可比矩阵
        self.nx = 0          # 状态变量个数
        self.ny = 0          # 代数变量个数，默认为0
        self.n_bus = 0        # 母线节点数，默认为0
        self.factorize = []   # 如果是True，则对雅可比矩阵进行因式分解