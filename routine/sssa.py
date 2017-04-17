import system
from numpy.linalg import eigvals
from cvxopt.umfpack import linsolve


def state_matrix():
    Gyx = matrix(system.DAE.Gx)
    linsolve(system.DAE.Gy, Gyx)
    return system.DAE.Fx - system.DAE.Fy * Gyx

def eigs():
    As = state_matrix()
    return eigvals(As)