import system
import numpy as np

def calcInc():
     Gcall()
     Gycall()
     y=np.linalg.solve(system.DAE.Gy,system.DAE.g)   #直接调用linalg中的solve求解修正方程
     return -y
def powerflow():
    """main power flow routine"""
    # general settings
    Gcall()
    iteration = 1
    iter_max = system.Settings.pf_max_iter
    convergence = True            #收敛
    tol = system.Settings.tol


    # main loop
    while max(abs(dae.g)) > tol and iteration <= iter_max:
         inc = calcInc()
         system.DAE.y += inc
         iteration += 1
         Gcall()
         # stop if the error increases too much
    if iteration > iter_max:
        print ('Reached maximum number of iterations')
        convergence = False
    else:
        #先求PV节点无功
        # V[i] = system.PV.__dict__['_data']['V0'][i - 1]
        # V[j] = system.PV.__dict__['_data']['V0'][j - 1]
        # G[i][j] =
        # B[i][j] =
        # Va[i] =
        # Va[j] =      #对应pv节点的电压相角
        # P[i]
        # Q[i]        #对应pv节点的P,Q
        sumq[i]=0
        for i in range(pv.self.n):
            for j in range(bus.self.n):
                if i~j:   #j与i直接相连
                    sumq[i] = V[j] * (G[i][j] *sin(Va[i]-Va[j]) - B[i][j]* cos(Va[i]-Va[j]))
            Q[i]=V[i]*sumq[i]

        #求平衡节点有无功
        sumps[bus.self.n] = 0
        sumqs[bus.self.n] = 0
        for j in range(bus.self.n):
            if i~j:       #j与i直接相连
                sumps[bus.self.n] = V[j] * (G[bus.self.n][j] * cos(Va[bus.self.n] - Va[j]) + B[bus.self.n][j] * sin(Va[bus.self.n] - Va[j]))
                sumqs[bus.self.n] = V[j] * (G[bus.self.n][j] *  sin(Va[bus.self.n] - Va[j]) - B[bus.self.n][j] * cos(Va[bus.self.n] - Va[j]))
        Ps[bus.self.n] = V[bus.self.n] * sumps[bus.self.n]
        Qs[bus.self.n] = V[bus.self.n] * sumqs[bus.self.n]
        print('result')