#!/usr/bin/env Python
# coding=utf-8
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
import scipy.io as sio
import time
import os
# 开始时间
def readdata(datafile):

    for i in system.Device.models:
        string = 'system.' + i +'._init_data()'
        exec(compile(string, '', 'exec'))




    # system.Bus._init_data()
    # system.PV._init_data()
    # system.PQ._init_data()
    # system.SW._init_data()
    # system.Shunt._init_data()
    # system.Line._init_data()
    # system.Syn6._init_data()
    # system.Avr1._init_data()
    # system.Avr2._init_data()
    # system.Fault._init_data()
    # system.Wind._init_data()
    # system.Dfig._init_data()

    case = open(datafile)
    for each_line in case:

        data = each_line.split()


        if data[0] == 'Bus':

            bus = data[0] + '_' + str(data[1])
            Vb = float(data[2])
            V_0 = float(data[3])
            theta0 = float(data[4])
            bus_case = {'bus': bus, 'Vb': Vb, 'V_0': V_0, 'theta0': theta0}
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
            u = int(data[16])
            Line_case = {'f_bus': f_bus, 'to_bus': to_bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn,
                         'l': l, 'kT': kT, 'r': r, 'x': x, 'b': b, 'tap_ratio': tap_ratio,
                         'theta': theta, 'Imax': Imax, 'Pmax': Pmax, 'Smax': Smax, 'u': u}
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

        if data[0] == 'Wind':
            type = int(data[1])
            vwn = float(data[2])
            rho = float(data[3])
            Tw = float(data[4])
            deltat = float(data[5])
            cw = float(data[6])
            kw = float(data[7])

            Wind_case = {'type': type, 'vwn': vwn, 'rho': rho, 'Tw': Tw, 'deltat': deltat, 'cw': cw, 'kw': kw}
            system.Wind.add(**Wind_case)

        if data[0] == 'Dfig':
            bus = 'Bus_' + str(data[1])
            windnum = int(data[2]) - 1
            Sn = float(data[3])
            Vn = float(data[4])
            fn = float(data[5])
            rs = float(data[6])
            xs = float(data[7])
            rr = float(data[8])
            xr = float(data[9])
            xm = float(data[10])
            Hm = float(data[11])
            Kp = float(data[12])
            Tp = float(data[13])
            Kv = float(data[14])
            Te = float(data[15])
            R = float(data[16])
            np1 = float(data[17])
            nb = float(data[18])
            nGB = float(data[19])
            pmax = float(data[20])
            pmin = float(data[21])
            qmax = float(data[22])
            qmin = float(data[23])
            ng = float(data[24])

            Dfig_case = {'Sn': Sn, 'Vn': Vn, 'fn': fn, 'bus': bus, 'wind': windnum, 'rs': rs, 'xs': xs,
                         'rr': rr, 'xr': xr, 'xm': xm,
                         'Hm': Hm, 'Kp': Kp, 'Tp': Tp, 'Kv': Kv, 'Te': Te, 'R': R, 'np': np1, 'nb': nb,
                         'pmax': pmax, 'pmin': pmin, 'qmax': qmax, 'qmin': qmin, 'nGB': nGB, 'ng': ng}
            system.Dfig.add(**Dfig_case)

        if data[0] == 'Breaker':

            line = int(data[1])
            bus = int(data[2])
            Sn = float(data[3])
            Vn = float(data[4])
            fn = float(data[5])
            u = int(data[6])
            t1 = float(data[7])
            t2 = float(data[8])
            u1 = int(data[9])
            u2 = int(data[10])
            breaker_case = {'line': line, 'bus': bus, 'Sn': Sn, 'Vn': Vn, 'fn': fn, 'u': u, 't1': t1, 't2': t2, 'u1': u1, 'u2': u2}
            system.Breaker.add(**breaker_case)

        if data[0] == 'Pss2':
            avr = int(data[1]) - 1
            m_model = float(data[2])
            ni = float(data[3])
            vsmax = float(data[4])
            vsmin = float(data[5])
            Kw = float(data[6])
            Tw = float(data[7])
            T1 = float(data[8])
            T2 = float(data[9])
            T3 = float(data[10])
            T4 = float(data[11])
            Ka = float(data[12])
            Ta = float(data[13])
            Kp = float(data[14])
            Kv = float(data[15])
            vamax = float(data[16])
            va1max = float(data[17])
            vs1max = float(data[18])
            vs1min = float(data[19])
            ethr = float(data[20])
            wthr = float(data[21])
            s2 = float(data[22])
            u = int(data[23])
            Pss2_case = {'avr': avr, 'm_model': m_model, 'ni': ni, 'vsmax': vsmax, 'vsmin': vsmin, 'Kw': Kw,
                         'Tw': Tw, 'T1': T1, 'T2': T2, 'T3': T3, 'T4': T4, 'Ka': Ka, 'Ta': Ta,
                         'Kp': Kp, 'Kv': Kv, 'vamax': vamax, 'va1max': va1max, 'vs1max': vs1max,
                         'vs1min': vs1min, 'ethr': ethr, 'wthr': wthr, 's2': s2, 'u': u}
            system.Pss2.add(**Pss2_case)


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
    system.Line.base()
    system.Line.build_y()
    # 设置线路最大传输容量
    system.Line.SetSmax()

    # Fault
    system.Fault._bus_index()
    system.Fault._list2matrix()
    system.Fault.base(Vb=system.Bus.Vb[system.Fault.a])
    system.Fault.setup()

    # call函数
    system.Device.setup()

if __name__ == "__main__":

    path = os.path.abspath('..')
    datafile = path + '\\text\d_014_wind.txt'
    readdata(datafile)
