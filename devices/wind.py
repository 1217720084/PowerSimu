#!/usr/bin/env Python
# coding=utf-8

from devices.base_device import base_device
import system
from cvxopt.base import matrix, spmatrix, mul, div
import numpy
import math
import random
from scipy import interpolate
import scipy.io as sio
import os

class wind(base_device):

    def __init__(self):
        base_device.__init__(self)
        self._data.update(
            {'type': 2, 'vwn': 16,  'rho': 1.225, 'Tw': 0.5, 'deltat': 0.01, 'cw': 20, 'kw': 2})
        self._type = 'Wind'
        self._name = 'Wind'
        self.n = 0

        self._data.update({'vwa': 1, 'time': [], 'svw': [], 'vwt': -1.0})
        self._algebs = ['ws']
        self._states = ['vw']
        self._params.extend(['vwn', 'rho', 'deltat', 'cw', 'kw', 'Tw', 'rho', 'vw'])
        #self.n_PV = 0
        #self._voltages = ['V0']
        #self._powers = ['Pg']

        self.properties.update({'gcall': True, 'Gycall': True,
                                'fcall': True, 'Fxcall': True})

    def setx0(self):

        if self.n == 0:
            return
        check = 1
        # 初始化平均风速（设为1）
        self.vwa = system.DAE.x[self.vw]
        system.DAE.y[self.ws] = system.DAE.x[self.vw]

        # 确保每台风机和一个风速连接

        # 初始化
        for i in range(self.n):
            t0 = system.Settings.t0
            tf = system.Settings.tf
            dt = self.deltat
            self.time[i].extend(numpy.linspace(t0, tf, (tf-t0)/dt[i]+1, endpoint=True))
            self.time[i] = list(self.time[i])

            # 测量数据
            if self.type[i] == 1:
                if os.path.exists('C://Users//user//Desktop//平台//余伟洲//仿真结果//windspeed1.mat'):
                    data = sio.loadmat('C://Users//user//Desktop//平台//余伟洲//仿真结果//windspeed1.mat')
                    told = data['time']
                    vold = data['vw']
                    told1 = [0.0]*len(told)
                    vold1 = [0.0]*len(told)
                    for j in range(len(told)):
                        told1[j] = told[j][0]
                        vold1[j] = vold[j][0]

                    f = interpolate.interp1d(told1, vold1)
                    self.svw[i] = f(self.time[i])/self.vwn[i]
                    self.svw[i][0] = self.vwa[i]

                else:
                    n = len(self.time[i])

                    if self.cw[i] <= 0:
                        self.cw[i] = 5

                    if self.kw[i] <= 0:
                        self.kw[i] = 2

                    self.svw[i].extend((-math.log(random.random()) / self.cw[i]) ** (1 / self.kw[i]) for _ in range(n))
                    mean_vw = sum(self.svw[i]) / n
                    self.svw[i][0] = self.vwa[i]
                    for j in range(n - 1):
                        self.svw[i][j + 1] = abs(1 + self.svw[i][j + 1] - mean_vw) * self.vwa[i]



            # Weibull distribution
            if self.type[i] == 2:

                n = len(self.time[i])

                if self.cw[i] <= 0:
                    self.cw[i] = 5

                if self.kw[i] <= 0:
                    self.kw[i] = 2

                self.svw[i].extend((-math.log(random.random())/self.cw[i]) ** (1/self.kw[i]) for _ in range(n))
                mean_vw = sum(self.svw[i])/n
                self.svw[i][0] = self.vwa[i]
                for j in range(n-1):
                    self.svw[i][j+1] = abs(1+self.svw[i][j+1]-mean_vw)*self.vwa[i]
                # self.svw[i][1:] = abs(1+self.svw[i][1:]-mean_vw)*self.vwa[i]


    def gcall(self):
        # system.DAE.lambda = 1
        if self.n == 0:
            return

        Vw = system.Wind.wspeed()
        for i in range(self.n):

            system.DAE.g[self.ws[i]] = Vw[i] - system.DAE.y[self.ws[i]]

    def wspeed(self):

        if self.n == 0:
            return
        Vw = [0.0] * self.n
        t = system.DAE.t
        if t < 0:
            t = system.Settings.t0
        if t == system.Settings.t0:
            Vw = self.vwa
        else:
            for i in range(self.n):
                f = interpolate.interp1d(self.time[i], self.svw[i])
                Vw[i] = f(t)
        return Vw

    def Gycall(self):

        if self.n == 0:
            return
        system.DAE.Gy = system.DAE.Gy - spmatrix(1, self.ws, self.ws, (system.DAE.ny, system.DAE.ny))

    def fcall(self):

        if self.n == 0:
            return

        system.DAE.f[self.vw] = div(system.DAE.y[self.ws]-system.DAE.x[self.vw], self.Tw)

    def Fxcall(self):

        if self.n == 0:
            return

        k = div(1, self.Tw)

        system.DAE.Fx = system.DAE.Fx - spmatrix(k, self.vw, self.vw, (system.DAE.nx, system.DAE.nx))
        system.DAE.Fy = system.DAE.Fy + spmatrix(k, self.vw, self.ws, (system.DAE.nx, system.DAE.ny))
