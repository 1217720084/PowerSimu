#!/usr/bin/env Python
# encoding = utf8
"""

测试读取数据

"""
# import

import system
from routine.readdata import readdata
from routine.pf import pf
from routine import td3
from routine import td4

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
from multiprocessing import Pool,Process,Manager,Value

if __name__ == '__main__':
    # 定义循环次数
    n = 1
    path1 = 'C:\\Users\\user\Desktop\powercon2018\仿真结果\\14bus\案例2(全过程和只考虑初始的对比)\\agent'
    path2 = 'C:\\Users\\user\Desktop\powercon2018\仿真结果\\14bus\案例2(全过程和只考虑初始的对比)\\conventional'

    # 读取数据文件
    path = os.path.abspath('..')
    datafile = path + '\\text\d_014_wind.txt'
    readdata(datafile)

    # 执行潮流计算
    pf()

    # 记录潮流结果
    xa = copy.deepcopy(system.DAE.x)
    ya = copy.deepcopy(system.DAE.y)
    aPl = copy.deepcopy(system.PQ.Pl)
    aQl = copy.deepcopy(system.PQ.Ql)

    # 设置初始仿真时间
    system.DAE.t = system.Settings.t0

    # 随机负荷波动参数
    mup = copy.deepcopy(system.PQ.Pl)
    muq = copy.deepcopy(system.PQ.Ql)
    sigmap = 0.01
    sigmaq = 0.01
    # 风机风速随机波动参数
    sigmavw = 0.01
    # 输电线路故障模拟
    [MWs, MWr] = system.Line.flow(system.DAE.y)
    MWs = abs(MWs)
    # 随机波动模拟步长
    suijih = 0.01

    # 生成进程通信manager
    manager = Manager()

    # 负荷波动序列
    Pl = manager.list(system.PQ.Pl)
    Ql = manager.list(system.PQ.Ql)

    # 风速波动序列
    vw = manager.list(system.DAE.x[system.Wind.vw])
    # 线路随机故障模拟
    tf = manager.list(system.Line.tf)
    tc = manager.list(system.Line.tc)
    Sline = manager.list(MWs)

    # 对时时间
    sn = 1
    syst = Value('d', 0.0)

    # 设置子进程
    p1 = Process(target=system.PQ.suiji, args=(mup, muq, sigmap, sigmaq, suijih, syst, Pl, Ql))
    p2 = Process(target=system.Dfig.suiji, args=(sigmavw, suijih, vw, syst))
    p3 = Process(target=system.Line.setfault, args=(tf, tc, suijih, Sline, syst))

    # 子进程
    p1.daemon = True  # 当主进程结束时，该子进程也直接终止运行
    p2.daemon = True  # 当主进程结束时，该子进程也直接终止运行
    p3.daemon = True  # 当主进程结束时，该子进程也直接终止运行
    print('子进程开始:%s' % time.time())
    p1.start()
    p2.start()
    p3.start()

    time.sleep(1.0)

    for i in range(n):

        print('循环次数：%s' % str(i+1))
        print('父进程%s' % os.getpid())
        td3.td(htype='fixed', syst=syst, Pl=Pl, Ql=Ql, sn=sn, vw=vw, tf=tf, tc=tc, Sline=Sline)

        # 保存结果
        var = system.Varout.var
        t0 = matrix(system.Varout.t)
        var = matrix(var)
        var = var.T
        sio.savemat(path1 + '\\tdvar' + str(i+1) + '.mat', {'var': var, 't': t0})
        system.Varout.clear()

        # 恢复潮流前数据
        system.DAE.x = xa
        system.DAE.y = ya
        system.PQ.Pl = copy.deepcopy(aPl)
        system.PQ.Ql = copy.deepcopy(aQl)

        #
        aPl = np.array(var[:, 96:107])
        aQl = np.array(var[:, 107:118])
        avw = np.array(var[:, 40])
        atf = np.array(var[:, 118:138])
        atc = np.array(var[:, 138:])
        td4.td(htype='fixed', Pl=aPl, Ql=aQl, vw=avw, tf=atf, tc=atc)

        # 保存结果
        var = system.Varout.var
        t0 = matrix(system.Varout.t)
        var = matrix(var)
        var = var.T
        sio.savemat(path2 + '\\tdvar' + str(i + 1) + '.mat', {'var': var, 't': t0})
        system.Varout.clear()

        # 恢复潮流前数据
        system.DAE.x = xa
        system.DAE.y = ya
        system.PQ.Pl = copy.deepcopy(aPl)
        system.PQ.Ql = copy.deepcopy(aQl)

    # p1.join()
    # p2.join()
    # p3.join()
    print('子进程结束')