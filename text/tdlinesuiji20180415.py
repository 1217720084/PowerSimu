# _*_ coding:utf-8 _*_
from text import tdyu20180415
from cvxopt.base import spmatrix, sparse, matrix
import system
import copy
import scipy.io as sio
from time import clock

# 潮流计算
tdyu20180415.pf20180410()
sio.savemat('pf.mat', {'pf': system.DAE.y})
# 设置仿真参数
n = 100  # 仿真次数

tdstart = clock()
# 备份变量
tf = copy.deepcopy(system.Line.tf)
tc = copy.deepcopy(system.Line.tc)
lineu = copy.deepcopy(system.Line.u)
xa = copy.deepcopy(system.DAE.x)
ya = copy.deepcopy(system.DAE.y)
apm0 = copy.deepcopy(system.Syn6.pm0)
aPl = copy.deepcopy(system.PQ.Pl)
aQl = copy.deepcopy(system.PQ.Ql)
muvref01 = copy.deepcopy(system.Avr1.vref0)
muvref02 = copy.deepcopy(system.Avr2.vref0)
#
names = locals()

for i in range(n):

    tdyu20180415.td20180410()
    t0 = matrix(system.Varout.t)
    var = system.Varout.var
    # print(system.Varout.t)
    var = matrix(var)
    var = var.T
    names['time%s' % (i + 1)] = t0
    names['var%s' % (i + 1)] = var
    v0 = var[:, 76]
    system.Varout.clear()

    system.DAE.x = xa
    system.DAE.y = ya

    # system.Line.tf = copy.deepcopy(tf)
    # system.Line.tc = copy.deepcopy(tc)
    system.PQ.Pl = copy.deepcopy(aPl)
    system.PQ.Ql = copy.deepcopy(aQl)
    # system.Avr1.vref0 = copy.deepcopy(muvref01)
    # system.Avr2.vref0 = copy.deepcopy(muvref02)
    # system.Line.u = copy.deepcopy(lineu)
    # system.Syn6.pm0 = copy.deepcopy(apm0)
    # print(system.PQ.Vmin)


    print('第%i 次仿真' % (i+1))

for i in range(n):

    datanew = 'C://Users//user//Desktop//平台//余伟洲//仿真结果//pythonsimulation//'+'tdvar'+str(i+1)+'.mat'

    sio.savemat(datanew, {'var': names['var%s' % (i + 1)], 't': names['time%s' % (i + 1)]})
tdfinish = clock()
tdt = tdfinish - tdstart
print('时域仿真总时间：')
print(tdt)
sio.savemat('C://Users//user//Desktop//平台//余伟洲//仿真结果//pythonsimulation//tdvar.mat', {'var': var, 't': t0})
sio.savemat('C://Users//user//Desktop//平台//余伟洲//仿真结果//pythonsimulation//suiijpm0.mat', {'suijipm0': matrix(system.Syn6.suijipm0)})