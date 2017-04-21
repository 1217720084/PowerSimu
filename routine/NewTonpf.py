import system
import numpy as np    #用于底下的linalg函数

def calcInc():
    for self.a in range(dae.n_bus):
        for self.v in range(dae.n_bus):
            gcall()
    for self.a in range(dae.n_bus):
        for self.v in range(dae.n_bus):
            for self.s in range(dae.n_bus):
                for self.t in range(dae.n_bus):
                    Gycall()                     #J矩阵里位置应该还要调整


     y=np.linalg.solve(dae.Gy,dae.g)   #直接调用linalg中的solve求解修正方程
     return y                           #假定他存在一个类似于dae.y的结构里，不过这样表达应该不对
def powerflow():
    """main power flow routine"""
    # general settings
    for self.a in range(dae.n_bus):
        for self.v in range(dae.n_bus):
            gcall()
    iteration = 1
    iter_max = system.Settings.pf_max_iter
    convergence = True            #收敛
    tol = system.Settings.tol


    # main loop
    while max(abs(dae.g)) > tol and iteration <= iter_max:
         inc = calcInc()
         i = 0
         j = 1
         #一个元素一个元素的更新的
         while i < dae.n_bus:
             while j < dae.n_bus:
                 dae.y[i]-=inc[i]
                 dae.y[j]-=dae.y[j]*inc(j)dae.y       #注意inc里的Δv不是直接表示而是用Δv/v
                 i+=2
                 j+=2

         iteration += 1
         for self.a in range(dae.n_bus):
             for self.v in range(dae.n_bus):
                 gcall()
         # stop if the error increases too much
    if iteration > iter_max:
        print ('Reached maximum number of iterations')
        convergence = False
    else:
        #求pv节点的q
        sumq[i]=0

        for i in range(dae.n_pv):
            for j in range(dae.n_pv):
                    sumq[i]+= dae.y[j+1] *(system.DAE.Y_G[i][j] * cos(dae.y[i] - dae.y[j]) + system.DAE.Y_B[i][j] * sin(dae.y[i]-dae.y[j])
            Qpv[i]=V[i]*sumq[i]

        #求平衡节点有无功
        sumps= 0
        sumqs= 0
        for i in range(bus.self.n):
                sumps = dae.y[i+1] * (system.DAE.Y_G[bus,self.n][i] * cos(dae.y[bus.self.n] - dae.y[i]) + system.DAE.Y_G[bus,self.n][i] * sin(dae.y[bus.self.n] - dae.y[i])
                sumqs = dae.y[i+1] * (system.DAE.Y_G[bus,self.n][i] * sin(dae.y[bus.self.n] - dae.y[i]) - system.DAE.Y_G[bus,self.n][i] * cos(dae.y[bus.self.n])
        Ps[bus.self.n] = dae.y[bus.self.n+1] * sumps
        Qs[bus.self.n] = dae.y[bus.self.n+1] * sumqs
        print('result')