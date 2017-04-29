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


system.Bus._init_data()
system.PV._init_data()
system.PQ._init_data()
system.SW._init_data()
system.Shunt._init_data()
system.Line._init_data()

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
        system.Line.add(**Line_case)

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




# def calcInc():
#     system.Line.gcall()
#     system.PQ.gcall()
#     system.Shunt.gcall()
#     system.PV.gcall()
#     system.SW.gcall()
#     system.Line.Gycall()
#     system.Shunt.Gycall()
#     system.PV.Gycall()
#     system.SW.Gycall()
#
#
#
#     system.DAE.g = np.array(system.DAE.g)
#
#     y=np.linalg.solve(system.DAE.Gy,system.DAE.g)   #直接调用linalg中的solve求解修正方程
#     y=matrix(y)
#     return y


#

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
    A=sparse(system.DAE.Gy)
    inc=matrix(system.DAE.g)
    if system.DAE.factorize:
        F=symbolic(A)     #重新排列A矩阵以减少填充并执行LU分解，返回为可以传递的不透明 C object到umfpack.numeric（）。
        system.DAE.factorize = False
    try:
        N = numeric(A,F)
        solve(A,N,inc)
    except:
        print('unexpec')
        F=symbolic(A)
        solve(A,numeric(A,F),inc)

    return inc

#






    print(system.DAE.Gy)
    print(system.DAE.g)
    system.DAE.g = np.array(system.DAE.g)
    print(system.DAE.g)
    y=np.linalg.solve(system.DAE.Gy,system.DAE.g)   #直接调用linalg中的solve求解修正方程
    return y
# system.Line.gcall()
# system.PQ.gcall()
# system.Shunt.gcall()
# system.PV.gcall()
system.DAE.g = np.array(system.DAE.g)

iteration = 1
iter_max = system.Settings.iter
convergence = True  # 收敛
tol = system.Settings.error
cycle = True




    #main loop
while cycle or (max(abs(system.DAE.g)) > tol and iteration <= iter_max):

    inc = calcInc()
    cycle = False
    system.DAE.y-= inc
    iteration += 1

    system.DAE.g = np.array(system.DAE.g)

    # stop if the error increases too much
if iteration > iter_max:
    print ('Reached maximum number of iterations')
    convergence = False
print(iteration-1)
print(system.DAE.y)


