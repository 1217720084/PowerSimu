import system
import numpy as np

def calcInc():
     Gcall（）
     Gycall（）
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
    system.Settings.error = tol + 1

    # main loop
    while system.Settings.error > tol and iteration <= iter max:
         inc = calcInc()
         system.DAE.y += inc
         iteration += 1
         Gcall()
         # stop if the error increases too much
    if iteration > iter_max:
        print ('Reached maximum number of iterations')
        convergence = False
    else:
