# _*_ coding:utf-8 _*_
"""

测试读取数据

"""
# import

import system
from cvxopt.base import spmatrix, sparse,matrix
from cvxopt.umfpack import numeric,symbolic,solve    #模块cvxopt.umfpack包含四个用于求解稀疏非对称线性方程组的函数
from numpy import multiply
import numpy as np
from time import clock
from numpy.linalg import eigvals
from cvxopt.umfpack import linsolve

# 开始时间
start = clock()

system.Bus._init_data()
system.PV._init_data()
system.PQ._init_data()
system.SW._init_data()
system.Shunt._init_data()
system.Line._init_data()
system.Syn6._init_data()
system.Avr1._init_data()
system.Avr2._init_data()

case = open('text_d_036.txt')
for each_line in case:

    data = each_line.split()


    if data[0] == 'Bus':

        bus = data[0] + '_' + str(data[1])
        Vb = float(data[2])
        bus_case = {'bus': bus, 'Vb': Vb}
        system.Bus.add(idx=bus, **bus_case)

    if data[0] == 'PV':

        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        Pg = float(data[4])
        V0 = float(data[5])
        qgmax = float(data[6])
        qgmin = float(data[7])
        Vmax = float(data[8])
        Vmin = float(data[9])
        PV_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'Pg': Pg, 'V0': V0,
                   'qgmax': qgmax, 'qgmin': qgmin, 'Vmax': Vmax, 'Vmin': Vmin}
        system.PV.add(**PV_case)

    if data[0] == 'PQ':

        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        Pl = float(data[4])
        Ql = float(data[5])
        Vmax = float(data[6])
        Vmin = float(data[7])
        PQ_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'Pl': Pl, 'Ql': Ql,
                   'Vmax': Vmax, 'Vmin': Vmin}
        system.PQ.add(**PQ_case)

    if data[0] == 'SW':

        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        V0 = float(data[4])
        Va = float(data[5])
        qgmax = float(data[6])
        qgmin = float(data[7])
        Vmax = float(data[8])
        Vmin = float(data[9])
        SW_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'V0': V0, 'Va': Va,
                   'qgmax': qgmax, 'qgmin': qgmin, 'Vmax': Vmax, 'Vmin': Vmin}
        system.SW.add(**SW_case)

    if data[0] == 'Shunt':

        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        fn = float(data[4])
        g = float(data[5])
        b = float(data[6])
        Shunt_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn, 'g': g, 'b': b}
        system.Shunt.add(**Shunt_case)

    if data[0] == 'Line':
        f_bus = 'Bus_' + str(data[1])
        to_bus = 'Bus_' + str(data[2])
        Sn = float(data[3])
        Vn = float(data[4])
        fn = float(data[5])
        l = float(data[6])
        kT = float(data[7])
        r = float(data[8])
        x = float(data[9])
        b = float(data[10])
        tap_ratio = float(data[11])
        theta = float(data[12])
        Imax = float(data[13])
        Pmax = float(data[14])
        Smax = float(data[15])
        Line_case = {'f_bus': f_bus, 'to_bus': to_bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn,
                     'l': l, 'kT': kT, 'r': r, 'x': x,'b': b, 'tap_ratio': tap_ratio,
                     'theta': theta, 'Imax': Imax, 'Pmax': Pmax, 'Smax': Smax}
        system.Line.add(**Line_case)\

    if data[0] == 'Syn6':
        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        fn = float(data[4])
        m_model = float(data[5])
        xl = float(data[6])
        ra = float(data[7])
        xd = float(data[8])
        xd1 = float(data[9])
        xd2 = float(data[10])
        Td01 = float(data[11])
        Td02 = float(data[12])
        xq = float(data[13])
        xq1 = float(data[14])
        xq2 = float(data[15])
        Tq01 = float(data[16])
        Tq02 = float(data[17])
        M = float(data[18])   # M = 2H
        D = float(data[19])
        Syn6_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn, 'm_model': m_model,
                     'xl': xl, 'ra': ra, 'xd': xd, 'xd1': xd1, 'xd2': xd2, 'Td01': Td01, 'Td02': Td02,
                     'xq': xq, 'xq1': xq1, 'xq2': xq2, 'Tq01': Tq01, 'Tq02': Tq02, 'M': M, 'D': D}
        system.Syn6.add(**Syn6_case)

    if data[0] == 'Avr2':
        bus = 'Bus_' + str(data[1])
        Type = float(data[2])
        vrmax = float(data[3])
        vrmin = float(data[4])
        Ka = float(data[5])
        Ta = float(data[6])
        Kf = float(data[7])
        Tf = float(data[8])
        Ke = float(data[9])
        Te = float(data[10])
        Tr = float(data[11])
        Ae = float(data[12])
        Be = float(data[13])
        Avr2_case = {'bus': bus, 'Type': Type, 'vrmax': vrmax, 'vrmin': vrmin, 'Ka': Ka, 'Ta': Ta, 'Kf': Kf, 'Tf': Tf,
                     'Ke': Ke, 'Te': Te, 'Tr': Tr, 'Ae': Ae, 'Be': Be}
        system.Avr2.add(**Avr2_case)

    if data[0] == 'Avr1':
        bus = 'Bus_' + str(data[1])
        Type = float(data[2])
        vrmax = float(data[3])
        vrmin = float(data[4])
        K0 = float(data[5])
        T1 = float(data[6])
        T2 = float(data[7])
        T3 = float(data[8])
        T4 = float(data[9])
        Te = float(data[10])
        Tr = float(data[11])
        Ae = float(data[12])
        Be = float(data[13])
        Avr1_case = {'bus': bus, 'Type': Type, 'vrmax': vrmax, 'vrmin': vrmin, 'K0': K0, 'T1': T1, 'T2': T2, 'T3': T3,
                     'T4': T4, 'Te': Te, 'Tr': Tr, 'Ae': Ae, 'Be': Be}
        system.Avr1.add(**Avr1_case)

case.close()

# Bus
system.Bus._xy_index()
system.Bus._bus_index()
system.Bus._list2matrix()
system.Bus.yinit(system.DAE)


# PV
system.PV._bus_index()
system.PV._list2matrix()
system.PV.base(Vb=system.Bus.Vb[system.PV.a])
system.PV.yinit(system.DAE)
system.PV._matrix2list()


# PQ
system.PQ._bus_index()
system.PQ._list2matrix()
system.PQ.base(Vb=system.Bus.Vb[system.PQ.a])
system.PQ.yinit(system.DAE)
system.PQ._matrix2list()


# SW
system.SW._bus_index()
system.SW.yinit(system.DAE)


# Shunt
system.Shunt._bus_index()


# Line
system.Line._bus_index()

system.Line.build_y()

#system.DAE.Y = sparse(system.DAE.Y)
#print(system.DAE.Y)
system.Line.gcall()
system.PQ.gcall()
system.Shunt.gcall()
system.PV.gcall()
system.SW.gcall()
print(system.DAE.g)
print(system.DAE.y)


def calcInc():
    global F
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

system.DAE.g = np.array(system.DAE.g)

iteration = 1
iter_max = system.Settings.iter
convergence = True  # 收敛
tol = system.Settings.error
cycle = True
err = []




    #main loop
while cycle or (max(abs(system.DAE.g)) > tol and iteration <= iter_max):

    inc = calcInc()
    cycle = False
    system.DAE.y -= inc
    err.append(max(abs(system.DAE.g)))
    print('第%i次迭代最大误差为：<%f>' % (iteration, err[iteration-1]))
    print(system.DAE.y)
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
# print('Bus.Pg、Qg')
# print(system.Bus.Pg)
# print(system.Bus.Qg)


# 测试Syn6
system.Syn6._bus_index()
system.Syn6._xy_index()

# 测试Avr
system.Avr1._bus_index()
system.Avr2._bus_index()
system.Avr1.getbus()
system.Avr2.getbus()
system.Avr1._xy_index()
system.Avr2._xy_index()

# 重新生成对应维度的x, y, g, f
system.DAE.x = [0.0] * system.DAE.nx
system.DAE.f = [0.0] * system.DAE.nx
system.DAE.y = list(system.DAE.y)
system.DAE.g = list(system.DAE.g)

# 重新生成雅可比矩阵
system.DAE.Gy = matrix(1.0, (system.DAE.ny, system.DAE.ny))
system.DAE.Fx = matrix(1.0, (system.DAE.nx, system.DAE.nx))
system.DAE.Fy = matrix(1.0, (system.DAE.nx, system.DAE.ny))
system.DAE.Gx = matrix(1.0, (system.DAE.ny, system.DAE.nx))



newy = [0]*(system.DAE.ny-system.Bus.n * 2)
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


# 重新生成微分代数方程和雅可比矩阵

system.DAE.g = matrix(0.0, (system.DAE.ny, 1))

# call
system.DAE._list2matrix()
system.Line.gcall()
system.PQ.gcall()
system.Shunt.gcall()
system.Syn6.gcall() # 测试到这里
system.Avr1.gcall()
system.Avr2.gcall()
system.PV.gcall()
#system.SW.gcall()
system.Line.Gycall()
system.Shunt.Gycall()
system.Syn6.Gycall()
system.Avr1.Gycall()
system.Avr2.Gycall()
system.PV.Gycall()
system.SW.Gycall()

system.Syn6.fcall()
system.Avr1.fcall()
system.Avr2.fcall()

system.Syn6.Fxcall()
system.Avr1.Fxcall()
system.Avr2.Fxcall()

# 生成状态矩阵

def state_matrix():

    Gyx = matrix(system.DAE.Gx)
    linsolve(sparse(system.DAE.Gy), Gyx)
    return system.DAE.Fx - system.DAE.Fy * Gyx

state_matrix()

# 计算特征值

def eigs():
    As = state_matrix()
    return eigvals(As)
eigen = eigs()

print(matrix(eigs()))