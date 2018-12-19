#!/usr/bin/env Python
# coding=utf-8
"""

测试读取数据

"""
# import

import system
from routine.readdata import readdata
from routine.pf import pf

import math
from cvxopt.base import spmatrix, sparse,matrix
from cvxopt.umfpack import numeric,symbolic,solve,linsolve    #模块cvxopt.umfpack包含四个用于求解稀疏非对称线性方程组的函数

from cvxopt.umfpack import solve as umfsolve
from numpy import multiply
import numpy as np
from time import clock
from numpy.linalg import eigvals
from cvxopt.umfpack import linsolve
import scipy.io as sio
import time
import os
import copy
from multiprocessing import Pool,Process,Manager



def time_step(convergency, iteration, t):
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
        AS = state_matrix()
        freq = max(abs(AS))
    else:  # 后续修改
        freq = 40
    if freq > system.Settings.freq:
        freq = float(system.Settings.freq)
    if not freq: freq = 40.0

    # 设置最小时间步长
    deltaT = abs(system.Settings.tf - system.Settings.t0)
    Tstep = 1 / freq
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


def td(htype, Pl, Ql, vw, tf, tc):
    # 输出采用的时域仿真算法
    tstart = time.time()
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
    system.Syn6.gcall()
    system.Avr2.gcall()
    system.Dfig.gcall()
    system.Wind.gcall()
    system.PV.gcall()
    system.SW.gcall()

    system.Line.Gycall()
    system.PQ.Gycall()
    system.Syn6.Gycall()
    system.Avr2.Gycall()
    system.Dfig.Gycall()
    system.Wind.Gycall()
    system.PV.Gycall()
    system.SW.Gycall()

    system.Syn6.fcall()
    system.Avr2.fcall()
    system.Dfig.fcall()
    system.Wind.fcall()

    if system.DAE.nx > 0:
        system.DAE.Fx = spmatrix([], [], [], (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = spmatrix([], [], [], (system.DAE.nx, system.DAE.ny))
        system.DAE.Gx = spmatrix([], [], [], (system.DAE.ny, system.DAE.nx))

    system.DAE.Fx = matrix(0.0, (system.DAE.nx, system.DAE.nx))
    system.DAE.Fy = matrix(0.0, (system.DAE.nx, system.DAE.ny))
    system.DAE.Gx = matrix(0.0, (system.DAE.ny, system.DAE.nx))

    system.Syn6.Fxcall()
    system.Avr2.Fxcall()
    system.Dfig.Fxcall()
    system.Wind.Fxcall()

    system.DAE.Gy = sparse(system.DAE.Gy)
    system.DAE.Fx = sparse(system.DAE.Fx)
    system.DAE.Fy = sparse(system.DAE.Fy)
    system.DAE.Gx = sparse(system.DAE.Gx)

    # print(system.DAE.f)

    system.DAE.tn = system.DAE.f

    # 初始化

    t = system.Settings.t0
    k = 1
    system.Varout.store(t, k)

    # 步长选择
    if htype == 'fixed':
        h = 0.01
    else:
        h = firsttime_tstep()

    inc = matrix(0.0, (system.DAE.nx + system.DAE.ny, 1))

    # 故障时间排序
    fixed_times = []
    fixed_times = system.Fault.gettimes()

    # fixed_times = sorted(fixed_times)

    # 计算最大转子间相对摇摆（先忽略）
    def anglediff():
        diff_max = 0

        delta = system.DAE.x[system.Syn6.delta].real()

        delta_diff = abs(delta - min(delta))
        diff_max = (max(delta_diff) * 180 / math.pi) > system.Settings.deltadelta

    diff_max = 0  # diff_max = anglediff()

    switch = False

    system.DAE.factorize = True

    # print(MWs)
    # 随机负荷波动参数
    mup = copy.deepcopy(system.PQ.Pl)
    muq = copy.deepcopy(system.PQ.Ql)
    sigmap = 0.01
    sigmaq = 0.01
    # 风力发电机随机风速参数
    suijih = 0.01
    sigmavw = 0.01
    ssvw= system.Wind.svw

    # 发电机随机出力
    # system.Dfig.suiji(sigmavw, t, suijih, ssvw, '2')
    # yu 随机负荷波动
    # system.PQ.suiji(mup, muq, sigmap, sigmaq, t, suijih, '2')
    # 发电机机械功率波动
    # system.Syn6.suiji(musynp, sigmasynp, t, suijih)


    # 主循环


    while t < system.Settings.tf and t + h > t and not diff_max:

        t3 = time.time()
        if t + h > system.Settings.tf:
            h = system.Settings.tf - t
        # print(system.DAE.t)
        actual_time = t + h

        #  检查不要跳过扰动（后续修改为更简洁形式）
        if fixed_times is not None:

            for item in fixed_times:

                if item > t and item < t + h:
                    actual_time = item
                    h = actual_time - t
                    switch = True
                    break

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
            system.Line.intervention(actual_time)
            # system.Breaker.intervention(actual_time)
            switch = False
        ts1 = time.time()

        system.Line.tf = tf[k]
        system.Line.tc = tc[k]
        system.Line.gettimes(fixed_times)
        system.PQ.Pl = Pl[k]
        system.PQ.Ql = Ql[k]
        # system.DAE.x[system.Wind.vw] = vw[k]
        if system.Wind.n == 1:
            system.Wind.vwt[0] = vw[k]
        else:
            system.Wind.vwt = vw[k]

        # # 输电线路故障模拟
        # system.Line.settimes(t)
        # system.Line.gettimes(fixed_times)
        # # yu 随机负荷波动
        # system.PQ.suiji(mup, muq, sigmap, sigmaq, suijih, t, type='2')
        # # 风力发电机随机风速波动
        # vw = system.DAE.x[system.Wind.vw]
        # system.Dfig.suiji(sigmavw, suijih, vw, t, '2')
        ts2 = time.time()
        # print('一个步长随机波动模块仿真耗时：%s' % (ts2-ts1))
        # 牛拉循环
        system.Settings.error = tol + 1  # 强迫进入循环一次
        td1 = time.time()
        while system.Settings.error > tol and iteration < iter_max:

            td5 = time.time()
            system.Line.gcall()
            system.PQ.gcall()
            system.Syn6.gcall()
            system.Avr2.gcall()
            system.Dfig.gcall()

            system.Wind.gcall()
            system.PV.gcall()
            system.SW.gcall()

            system.Line.Gycall()
            system.PQ.Gycall()
            system.Syn6.Gycall()
            system.Avr2.Gycall()
            system.Dfig.Gycall()
            system.Wind.Gycall()
            system.PV.Gycall()
            system.SW.Gycall()

            system.Syn6.fcall()
            system.Avr2.fcall()
            system.Dfig.fcall()
            system.Wind.fcall()

            if system.DAE.nx > 0:
                system.DAE.Fx = spmatrix([], [], [], (system.DAE.nx, system.DAE.nx))
                system.DAE.Fy = spmatrix([], [], [], (system.DAE.nx, system.DAE.ny))
                system.DAE.Gx = spmatrix([], [], [], (system.DAE.ny, system.DAE.nx))

            system.DAE.Fx = matrix(0.0, (system.DAE.nx, system.DAE.nx))
            system.DAE.Fy = matrix(0.0, (system.DAE.nx, system.DAE.ny))
            system.DAE.Gx = matrix(0.0, (system.DAE.ny, system.DAE.nx))

            system.Syn6.Fxcall()
            system.Avr2.Fxcall()
            system.Dfig.Fxcall()
            system.Wind.Fxcall()

            system.DAE.Gy = sparse(system.DAE.Gy)
            system.DAE.Fx = sparse(system.DAE.Fx)
            system.DAE.Fy = sparse(system.DAE.Fy)
            system.DAE.Gx = sparse(system.DAE.Gx)
            td6 = time.time()
            # print('生成雅可比矩阵时间:%s' % (td6-td5))

            # 计算雅可比矩阵
            if system.Settings.method == 1:  # 采用欧拉法
                system.DAE.Ac = sparse(
                    [[identica - h * system.DAE.Fx, system.DAE.Gx], [-h * system.DAE.Fy, system.DAE.Gy]])
                system.DAE.tn = system.DAE.x - xa - h * system.DAE.f
            elif system.Settings.method == 2:  # 采用隐式梯形法

                system.DAE.Ac = sparse(
                    [[identica - h * 0.5 * system.DAE.Fx, system.DAE.Gx], [-h * 0.5 * system.DAE.Fy, system.DAE.Gy]])

                system.DAE.tn = system.DAE.x - xa - h * 0.5 * (system.DAE.f + fn)

            # 限幅器
            system.Avr2.windup('td')
            system.Dfig.windup_final()
            # gg = -matrix([system.DAE.tn, system.DAE.g])
            # linsolve(system.DAE.Ac, gg)
            # inc = gg

            td3 = time.time()
            if system.DAE.factorize:
                F = symbolic(system.DAE.Ac)
                system.DAE.factorize = False
            inc = -matrix([matrix(system.DAE.tn), system.DAE.g])
            try:
                umfsolve(system.DAE.Ac, numeric(system.DAE.Ac, F), inc)
            except ArithmeticError:
                print('Singular matrix')
                iteration = iter_max + 1
            except ValueError:
                F = symbolic(system.DAE.Ac)
                try:
                    umfsolve(system.DAE.Ac, numeric(system.DAE.Ac, F), inc)
                except ArithmeticError:
                    print('Singular matrix')
                    iteration = iter_max + 1
            td4 = time.time()
            # print('一次牛拉迭代时间为:%s' % (td4-td3))
            # print('一次牛拉求解总耗时:%s' % (td4 - td5))





            system.DAE.x = system.DAE.x + inc[:system.DAE.nx]
            system.DAE.y = system.DAE.y + inc[system.DAE.nx:system.DAE.nx + system.DAE.ny]

            iteration = iteration + 1
            system.Settings.error = max(abs(inc))
        td2 = time.time()
        # print('一个步长仿真时间:%s' % (td2 - td1))
        # print('一个步长牛拉迭代次数：%s' % iteration)

        # 计算新的时间步长

        if iteration > iter_max:
            h = time_step(False, iteration, t)

            print('减小步长(delta t = %.5f s)' % h)
            system.DAE.x = matrix(xa)
            system.DAE.y = matrix(ya)
            system.DAE.f = matrix(fn)
        else:
            # h = time_step(True, iteration, t)
            # t = actual_time
            h = 0.01
            t = actual_time

        k = k + 1
        system.Varout.store(t, k)
        t4 = time.time()
        # print('一个步长总仿真耗时：%s' % (t4-t3))

    var = system.Varout.var
    t0 = matrix(system.Varout.t)

    var = matrix(var)
    var = var.T
    sio.savemat('C://Users//user//Desktop//平台//余伟洲//仿真结果//RES 时域仿真//restdvar.mat', {'var': var, 't': t0})
    # sio.savemat('C://Users//user//Desktop//restdvar.mat', {'var': var, 't': t0})
    print('time domain simulation completed')
    print('仿真次数%s', k)
    tend = time.time()
    print('机电暂态仿真总时间%s', tend - tstart)

if __name__ == '__main__':

    # 读取数据文件
    path = os.path.abspath('..')
    datafile = path + '\\text\d_014_wind.txt'
    readdata(datafile)

    # 执行潮流计算
    pf()

    data = sio.loadmat('restdvartd3.mat')

    Pl = data['var'][:, 96:107]
    Ql = data['var'][:, 107:118]
    # vw = data['var'][:, 40]
    tf = data['var'][:, 118:138]
    tc = data['var'][:, 138:158]
    vw = data['var'][:, 158]
    td(htype='fixed', Pl=Pl, Ql=Ql, vw=vw, tf=tf, tc=tc)

