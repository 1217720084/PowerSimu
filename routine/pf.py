#!/usr/bin/env Python
# coding=utf-8
"""

测试读取数据

"""
# import

import system
from routine.readdata import readdata

import math
from cvxopt.base import spmatrix, sparse, matrix
from cvxopt.umfpack import numeric, symbolic, solve, linsolve    #模块cvxopt.umfpack包含四个用于求解稀疏非对称线性方程组的函数

from cvxopt.umfpack import solve as umfsolve
from numpy import multiply
import numpy as np
from time import clock
from numpy.linalg import eigvals
from cvxopt.umfpack import linsolve
import scipy.io as sio
import time
import os


def calcInc():

    global F

    # exec(system.Device.call_pf)
    system.Line.gcall()
    system.PQ.gcall()
    system.Shunt.gcall()
    system.PV.gcall()
    system.SW.gcall()
    system.Line.Gycall()
    system.Shunt.Gycall()
    system.PV.Gycall()
    system.SW.Gycall()
    A = sparse(system.DAE.Gy)
    inc = matrix(system.DAE.g)
    if system.DAE.factorize:
        F = symbolic(A)     #重新排列A矩阵以减少填充并执行LU分解，返回为可以传递的不透明 C object到umfpack.numeric（）。
        system.DAE.factorize = False
    try:
        N = numeric(A, F)
        solve(A, N, inc)
    except:
        print('unexpect')
        F = symbolic(A)
        solve(A, numeric(A, F), inc)

    return inc

    # print(system.DAE.Gy)
    # print(system.DAE.g)
    # system.DAE.g = np.array(system.DAE.g)
    # print(system.DAE.g)
    # y=np.linalg.solve(system.DAE.Gy,system.DAE.g)   #直接调用linalg中的solve求解修正方程
    # return y

def pf():

    start = clock()

    system.Line.gcall()
    system.PQ.gcall()
    system.Shunt.gcall()
    system.PV.gcall()
    system.SW.gcall()

    system.DAE.g = np.array(system.DAE.g)

    iteration = 1
    iter_max = system.Settings.iter
    convergence = True  # 收敛
    tol = system.Settings.error
    cycle = True
    err = []

    # main loop
    while cycle or (max(abs(system.DAE.g)) > tol and iteration <= iter_max):
        inc = calcInc()
        cycle = False
        system.DAE.y -= inc
        err.append(max(abs(system.DAE.g)))
        # print('第%i次迭代最大误差为：<%f>' % (iteration, err[iteration-1]))
        # print(system.DAE.y)
        iteration += 1

        system.DAE.g = np.array(system.DAE.g)

        # stop if the error increases too much
    if iteration > iter_max:
        print('Reached maximum number of iterations')
        convergence = False

    # 结束时间

    finish = clock()

    t = finish - start
    print('潮流计算运行时间：%f' % t)

    # 计算Pl和Ql
    system.DAE.g = matrix(0.0, (system.DAE.ny, 1))
    # print(system.DAE.g)
    system.PQ.gcall()
    system.Shunt.gcall()
    system.DAE._list2matrix()
    system.Bus.Pl = system.DAE.g[system.Bus.a]
    # print(system.Bus.Pl)
    system.Bus.Ql = system.DAE.g[system.Bus.v]
    # print(system.Bus.Ql)

    # 计算Pg 和 Qg
    system.Line.gcall()
    system.PQ.gcall()
    system.Shunt.gcall()
    system.DAE._list2matrix()
    system.Bus.Pg = system.DAE.g[system.Bus.a]
    system.Bus.Qg = system.DAE.g[system.Bus.v]
    # print(system.Bus.Pg)
    # print(system.Bus.Qg)


    # 测试Syn6
    system.Syn6._bus_index()
    system.Syn6._dxy_index()

    # 测试Avr
    system.Avr1._bus_index()
    system.Avr2._bus_index()
    system.Avr1.getbus()
    system.Avr2.getbus()
    system.Avr2._dxy_index()
    system.Avr1._dxy_index()

    # 测试Pss2
    system.Pss2._bus_index()
    system.Pss2._dxy_index()

    # Wind
    system.Wind._dxy_index()

    system.Dfig._bus_index()
    system.Dfig._dxy_index()

    # 重新生成对应维度的x, y, g, f
    system.DAE.x = [1.0] * system.DAE.nx
    system.DAE.f = [0.0] * system.DAE.nx
    system.DAE.y = list(system.DAE.y)
    system.DAE.g = list(system.DAE.g)

    # 重新生成雅可比矩阵
    system.DAE.Gy = matrix(1.0, (system.DAE.ny, system.DAE.ny))
    system.DAE.Fx = matrix(1.0, (system.DAE.nx, system.DAE.nx))
    system.DAE.Fy = matrix(1.0, (system.DAE.nx, system.DAE.ny))
    system.DAE.Gx = matrix(1.0, (system.DAE.ny, system.DAE.nx))

    newy = [0] * (system.DAE.ny - system.Bus.n * 2)
    system.DAE.y.extend(newy)
    system.DAE.g.extend(newy)
    system.DAE.y = matrix(system.DAE.y)

    # 测试Syn6 setx0
    system.Syn6._list2matrix()
    system.Syn6.base(Vb=system.Bus.Vb[system.Syn6.a])
    system.Syn6.setx0()
    # print('Syn6 setx0')
    # print(system.DAE.x)
    # 测试Avr setx0

    system.Avr1._list2matrix()
    system.Avr2._list2matrix()
    # system.Avr1.base(Vb=system.Bus.Vb[system.Avr1.a])
    # system.Avr2.base(Vb=system.Bus.Vb[system.Avr2.a])
    system.Avr1.setx0()
    system.Avr2.setx0()
    # print('Avr setx0')
    # print(system.DAE.y)
    # print(system.Syn6.delta)
    # print('new')
    #
    # print(system.DAE.x)
    # print(system.DAE.x[system.Syn6.delta])

    # 测试Pss2
    system.Pss2._list2matrix()
    system.Pss2.setx0()

    # Wind初始化
    system.Dfig._list2matrix()
    system.Wind._list2matrix()
    system.Dfig.base()
    system.Dfig.setx0()

    system.Wind._list2matrix()
    system.Wind.setx0()

if __name__ == '__main__':

    path = os.path.abspath('..')
    datafile = path + '\\text\d_014_wind.txt'
    readdata(datafile)

    pf()
    print(system.DAE.y)
