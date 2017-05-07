import system
from numpy.linalg import eig
from cvxopt.lapack import gesv
from cvxopt.base import matrix, spmatrix, mul


def compute_eig(As):

    mu, N = eig(matrix(As))
    N = matrix(N)
    n =len(mu)
    idx = range(n)
    W = matrix(spmatrix(1.0, idx, idx, (n, n), N.typecode))
    gesv(N, W)
    pf = mul(abs(W.T), abs(N))
    b = matrix(1.0, (1, n))
    WN = b * pf
    pf = pf.T

    for item in idx:

        mur = mu[item].real
        mui = mu[item].imag
        mu[item] = complex(round(mur, 5), round(mui, 5))
        pf[item, :] /= WN[item]