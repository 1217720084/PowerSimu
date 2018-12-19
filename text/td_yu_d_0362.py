# _*_ coding:utf-8 _*_
"""

测试读取数据

"""
# import

import system
import math
from cvxopt.base import spmatrix, sparse,matrix
from cvxopt.umfpack import numeric,symbolic,solve,linsolve    #模块cvxopt.umfpack包含四个用于求解稀疏非对称线性方程组的函数

from cvxopt.umfpack import solve as umfsolve
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
system.Fault._init_data()

case = open('text_d_0361.txt')
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
        z = float(data[8])

        PQ_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'Pl': Pl, 'Ql': Ql,
                   'Vmax': Vmax, 'Vmin': Vmin, 'z': z}
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



    if data[0] == 'Fault':
        bus = 'Bus_' + str(data[1])
        Sn = float(data[2])
        Vn = float(data[3])
        fn = float(data[4])
        tf = float(data[5])
        tc = float(data[6])
        rf = float(data[7])
        xf = float(data[8])

        Fault_case = {'bus': bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn, 'tf': tf, 'tc': tc, 'rf': rf, 'xf': xf}
        system.Fault.add(**Fault_case)

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
# yu0406
system.Line.base()
system.Line.build_y()

#Fault
system.Fault._bus_index()
system.Fault._list2matrix()
system.Fault.base(Vb=system.Bus.Vb[system.Fault.a])
system.Fault.setup()

#system.DAE.Y = sparse(system.DAE.Y)
#print(system.DAE.Y)
system.Line.gcall()
system.PQ.gcall()
system.Shunt.gcall()
system.Fault.gcall()
system.PV.gcall()
system.SW.gcall()



def calcInc():
    global F
    system.Line.gcall()
    system.PQ.gcall()
    system.Shunt.gcall()
    system.PV.gcall()
    system.SW.gcall()
    system.Line.Gycall()
    system.PQ.Gycall()
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


# 重新生成对应维度的x, y, g, f
system.DAE.x = [0.0] * system.DAE.nx
system.DAE.f = [0.0] * system.DAE.nx
system.DAE.y = list(system.DAE.y)
system.DAE.g = list(system.DAE.g)

# 重新生成雅可比矩阵
system.DAE.Gy = matrix(0.0, (system.DAE.ny, system.DAE.ny))
system.DAE.Fx = matrix(0.0, (system.DAE.nx, system.DAE.nx))
system.DAE.Fy = matrix(0.0, (system.DAE.nx, system.DAE.ny))
system.DAE.Gx = matrix(0.0, (system.DAE.ny, system.DAE.nx))



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
# print(system.Syn6.delta)
# print('new')
#
# print(system.DAE.x)
# print(system.DAE.x[system.Syn6.delta])



#system.PQ.Vmin = system.DAE.y[system.PQ.v]



# 输出采用的时域仿真算法
print('time domain simulation ')
print('Trapezoidal integration method')

# 检查设置

iter_max = 20
tol = 0.00001
Dn = 1

if system.DAE.nx:
    Dn = system.DAE.nx

identica = spmatrix(1.0, range(max(Dn, 1)), range(max(Dn, 1)), (max(Dn, 1), max(Dn, 1)))

# 把PQ负荷转化为并联阻抗

system.PQ.pqshunt()

# 建立变量

system.DAE.t = system.Settings.t0

# call
system.DAE._list2matrix()
system.Line.gcall()
system.PQ.gcall()
system.Shunt.gcall()
system.Fault.gcall()
system.Syn6.gcall()
system.Avr1.gcall()
system.Avr2.gcall()
# system.PV.gcall()
# system.SW.gcall()

system.Line.Gycall()
system.PQ.Gycall()
system.Shunt.Gycall()
# 故障
system.Fault.Gycall()

system.Syn6.Gycall()
system.Avr1.Gycall()
system.Avr2.Gycall()
# system.PV.Gycall()
# system.SW.Gycall()

system.Syn6.fcall()
system.Avr1.fcall()
system.Avr2.fcall()

if system.DAE.nx > 0:
    system.DAE.Fx = spmatrix([], [], [], (system.DAE.nx, system.DAE.nx))
    system.DAE.Fy = spmatrix([], [], [], (system.DAE.nx, system.DAE.ny))
    system.DAE.Gx = spmatrix([], [], [], (system.DAE.ny, system.DAE.nx))

system.DAE.Fx = matrix(0.0, (system.DAE.nx, system.DAE.nx))
system.DAE.Fy = matrix(0.0, (system.DAE.nx, system.DAE.ny))
system.DAE.Gx = matrix(0.0, (system.DAE.ny, system.DAE.nx))

system.Syn6.Fxcall()
system.Avr1.Fxcall()
system.Avr2.Fxcall()

system.DAE.Gy = sparse(system.DAE.Gy)
system.DAE.Fx = sparse(system.DAE.Fx)
system.DAE.Fy = sparse(system.DAE.Fy)
system.DAE.Gx = sparse(system.DAE.Gx)

# print(system.DAE.f)

system.DAE.tn = system.DAE.f

# 初始化

t = system.Settings.t0
k = 1


# 步长选择

def state_matrix():

    Gyx = matrix(system.DAE.Gx)
    linsolve(sparse(system.DAE.Gy), Gyx)
    I = []
    J = []
    for i in range(system.DAE.nx):
        I.append(i)
        J.append(i)
    return system.DAE.Fx - system.DAE.Fy * Gyx - spmatrix(1e-6, I, J)



def firsttime_tstep():

    if system.DAE.nx == 0:
        freq = 1.0
    elif system.DAE.nx == 1:
        AS =state_matrix()
        freq = max(abs(AS))
    else:  #后续修改
        freq = 40
    if freq > system.Settings.freq:
        freq = float(system.Settings.freq)
    if not freq: freq = 40.0

    # 设置最小时间步长
    deltaT = abs(system.Settings.tf - system.Settings.t0)
    Tstep = 1/freq
    system.Settings.deltatmax = min(5 * Tstep, deltaT / 100.0)
    system.Settings.deltat = min(Tstep, deltaT / 100.0)
    system.Settings.deltatmin = min(Tstep / 64, system.Settings.deltatmax / 20)
    if system.Settings.fixt:
        if system.Settings.tstep <= 0:
            print('Fixed time step is negative or zero')
            print('Automatic time step has been set')
            system.Settings.fixt = 1
        elif system.Settings.tstep < system.Settings.deltatmin:
            print('Fixed time step is less than estimated minimum time step')
            system.Settings.deltat = system.Settings.tstep
        else:
            system.Settings.deltat = system.Settings.tstep
    return system.Settings.deltat



h = firsttime_tstep()

inc = matrix(0.0, (system.DAE.nx+system.DAE.ny, 1))

# 故障时间排序

fixed_times = system.Fault.gettimes()
fixed_times = sorted(fixed_times)

# 计算最大转子间相对摇摆（先忽略）
def anglediff():
    diff_max = 0

    delta = system.DAE.x[system.Syn6.delta].real()


    delta_diff = abs(delta - min(delta))
    diff_max = (max(delta_diff) * 180/math.pi) > system.Settings.deltadelta

diff_max = 0  # diff_max = anglediff()

switch = False

def time_step(convergency, iteration,t):
    if convergency == False:
        system.Settings.deltat = system.Settings.deltat * 0.5
        if system.Settings.deltat < system.Settings.deltatmin:
            system.Settings.deltat = 0


    if convergency == True:
        if iteration >= 15:
            system.Settings.deltat = max(system.Settings.deltat * 0.9, system.Settings.deltatmin)
        if iteration <= 10:
            system.Settings.deltat = min(system.Settings.deltat * 1.3, system.Settings.deltatmax)
        if system.Settings.fixt:
            system.Settings.deltat = min(system.Settings.tstep, system.Settings.deltat)

    if system.Fault.istime(t):
            system.Settings.deltat = min(system.Settings.deltat, 0.0025)


    return system.Settings.deltat

system.DAE.factorize = True

system.DAE.y[system.Bus.v[29]] = round(system.DAE.y[system.Bus.v[29]], 6)
vv = [system.DAE.y[system.Bus.v[29]]]
tt = [0.0]
pyiter = []
# 主循环
while t < system.Settings.tf and t+h > t and not diff_max:



    if t+h > system.Settings.tf:
        h = system.Settings.tf - t
    actual_time = t+h

    # 检查不要跳过扰动（后续修改为更简洁形式）
    for item in fixed_times:
        if item > t and item < t+h:
            actual_time = item
            h = actual_time - t
            switch = True
            break
    # print(actual_time)
    system.DAE.t = actual_time
    # 备份变量
    xa = matrix(system.DAE.x)
    ya = matrix(system.DAE.y)
    # 初始化牛拉法循环
    iteration = 1
    fn = matrix(system.DAE.f)  # 微分方程
    inc[0] = 1

    # 执行扰动

    if switch:

        system.Fault.intervention(actual_time)
        # system.Breaker.intervention(actual_time)
        switch = False

    # 牛拉循环
    system.Settings.error = tol + 1  # 强迫进入循环一次

    while system.Settings.error > tol and iteration < iter_max:
        system.Line.gcall()
        system.PQ.gcall()
        system.Shunt.gcall()
        system.Fault.gcall()
        system.Syn6.gcall()
        system.Avr1.gcall()
        system.Avr2.gcall()
        # system.PV.gcall()
        # system.SW.gcall()
        system.Line.Gycall()
        system.PQ.Gycall()
        system.Shunt.Gycall()
        # 故障
        system.Fault.Gycall()

        system.Syn6.Gycall()
        system.Avr1.Gycall()
        system.Avr2.Gycall()
        # system.PV.Gycall()
        # system.SW.Gycall()

        system.Syn6.fcall()
        system.Avr1.fcall()
        system.Avr2.fcall()

        if system.DAE.nx > 0:
            system.DAE.Fx = spmatrix([], [], [], (system.DAE.nx, system.DAE.nx))
            system.DAE.Fy = spmatrix([], [], [], (system.DAE.nx, system.DAE.ny))
            system.DAE.Gx = spmatrix([], [], [], (system.DAE.ny, system.DAE.nx))

        system.DAE.Fx = matrix(0.0, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = matrix(0.0, (system.DAE.nx, system.DAE.ny))
        system.DAE.Gx = matrix(0.0, (system.DAE.ny, system.DAE.nx))

        system.Syn6.Fxcall()
        system.Avr1.Fxcall()
        system.Avr2.Fxcall()

        system.DAE.Gy = sparse(system.DAE.Gy)
        system.DAE.Fx = sparse(system.DAE.Fx)
        system.DAE.Fy = sparse(system.DAE.Fy)
        system.DAE.Gx = sparse(system.DAE.Gx)

        # 计算雅可比矩阵
        if system.Settings.method == 1:  # 采用欧拉法
            system.DAE.Ac = sparse([[identica - h * system.DAE.Fx, system.DAE.Gx], [-h * system.DAE.Fy, system.DAE.Gy]])
            system.DAE.tn = system.DAE.x - xa - h * system.DAE.f
        elif system.Settings.method == 2:  # 采用隐式梯形法

            system.DAE.Ac = sparse(
                [[identica - h * 0.5 * system.DAE.Fx, system.DAE.Gx], [-h * 0.5 * system.DAE.Fy, system.DAE.Gy]])

            system.DAE.tn = system.DAE.x - xa - h * 0.5 * (system.DAE.f + fn)

        # 限幅器
        system.Avr2.windup('td')

        gg = -matrix([system.DAE.tn, system.DAE.g])
        linsolve(system.DAE.Ac, gg)
        inc = gg

        # if system.DAE.factorize:
        #     F = symbolic(system.DAE.Ac)
        #     system.DAE.factorize = False
        # inc = -matrix([system.DAE.tn, system.DAE.g])
        # try:
        #     umfsolve(system.DAE.Ac, numeric(system.DAE.Ac, F), inc)
        # except ArithmeticError:
        #     print('Singular matrix')
        #     iteration = iter_max+1
        # except ValueError:
        #     F = symbolic(system.DAE.Ac)
        #     try:
        #         umfsolve(system.DAE.Ac, numeric(system.DAE.Ac, F), inc)
        #     except ArithmeticError:
        #         print('Singular matrix')
        #         iteration = iter_max + 1





        system.DAE.x = system.DAE.x + inc[:system.DAE.nx]
        system.DAE.y = system.DAE.y + inc[system.DAE.nx:system.DAE.nx + system.DAE.ny]

        iteration = iteration + 1
        system.Settings.error = max(abs(inc))

    # 计算新的时间步长
    pyiter.append(iteration)
    print(iteration)
    if iteration > iter_max:
        h = time_step(False, iteration, t)

        print('减小步长(delta t = %.5f s)' % h)
        system.DAE.x = matrix(xa)
        system.DAE.y = matrix(ya)
        system.DAE.f = matrix(fn)
    else:
        h = time_step(True, iteration, t)
        t = actual_time


    k = k+1
    vv15 = round(system.DAE.y[system.Bus.v[29]].real, 6)

    # print(vv15)
    vv.append(vv15)
    tt.append(system.DAE.t)

tt = matrix(tt)
vv = matrix(vv)
print('tt')
tt = list(tt)
print(tt)
vv = list(vv)
print('vv')
print(vv)










print('')







