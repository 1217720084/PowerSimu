# _*_ coding:utf-8 _*_
from text import yu20180410
from cvxopt.base import spmatrix, sparse, matrix
import system
import copy
import scipy.io as sio

# 潮流计算
yu20180410.pf20180410()
# 设置仿真参数
n = 100 # 仿真次数
t = matrix(0.0, (1100, n))
v = matrix(0.0, (1100, n))
muvref = matrix(0.0, (1100, n))
# 备份变量
xa = copy.deepcopy(system.DAE.x)
ya = copy.deepcopy(system.DAE.y)
aPl = copy.deepcopy(system.PQ.Pl)
aQl = copy.deepcopy(system.PQ.Ql)
muvref01 = copy.deepcopy(system.Avr1.vref0)
muvref02 = copy.deepcopy(system.Avr2.vref0)
#
names = locals()

for i in range(n):

    [muvref0] = yu20180410.td20180410()
    t0 = matrix(system.Varout.t)
    var = system.Varout.var
    # print(system.Varout.t)
    var = matrix(var)
    var = var.T
    names['time%s' % (i + 1)] = t0
    names['var%s' % (i + 1)] = var
    v0 = var[:, 76]
    system.Varout.clear()
    print(list(system.DAE.y))
    system.DAE.x = xa
    system.DAE.y = ya
    print(list(system.DAE.y))
    system.PQ.Pl = copy.deepcopy(aPl)
    system.PQ.Ql = copy.deepcopy(aQl)
    system.Avr1.vref0 = copy.deepcopy(muvref01)
    system.Avr2.vref0 = copy.deepcopy(muvref02)
    # print(system.PQ.Vmin)
    for j in range(len(t0)):
        t[j, i] = t0[j]
        v[j, i] = v0[j]
        muvref[j, i] = muvref0[j]
    print('第%i 次仿真' % (i+1))
sio.savemat('va0.mat', {'t': t, 'v30': v, 'muvref': muvref})
sio.savemat('pytdvar.mat', {'t': t0, 'var': var})
for i in range(n):

    datanew = 'C://Users//user//Desktop//平台//余伟洲//仿真结果//pythonsimulation//'+'tdvar'+str(i+1)+'.mat'

    sio.savemat(datanew, {'var': names['var%s' % (i + 1)], 't': names['time%s' % (i + 1)]})
