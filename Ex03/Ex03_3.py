from math import sin, cos, sqrt, exp, pi
import numpy as np
import cmath as cm
import numpy.lib.scimath as sc
import scipy.optimize as opt
import matplotlib.pyplot as plt

hbar = 1.
m_e = 1.
a_0 = 1.
Hartree = 1.
eV = Hartree / 27.211386
Ang = a_0 / 0.529177210


def psi_n(t_n, r_n, k_, x_):
    return t_n * cm.exp(1j * k_ * x_) + r_n * cm.exp(-1j * k_ * x_)


def psiPrime_n(t_n, r_n, k_, x_):
    return 1j * k_ * (t_n * cm.exp(1j * k_ * x_) - r_n * cm.exp(-1j * k_ * x_))


def F(x, k, kw):
    return [psi_n(1., x[0], k, -L / 2.) - psi_n(x[1], x[2], kw, -L / 2.),
            psiPrime_n(1., x[0], k, -L / 2.) - psiPrime_n(x[1], x[2], kw, -L / 2.),
            psi_n(x[1], x[2], kw, L / 2.) - psi_n(x[3], 0, k, L / 2.),
            psiPrime_n(x[1], x[2], kw, L / 2.) - psiPrime_n(x[3], 0, k, L / 2.)]


def reals_as_complex(x):
    cIn = np.transpose([x[:-1:2], x[1::2]])
    cIn = [complex(c[0], c[1]) for c in cIn]
    return cIn


def f2p_as_reals(x, F, k, kw):
    cIn = reals_as_complex(x)
    cOut = F(cIn, k, kw)
    res = []
    for c in cOut:
        res.append(c.real)
        res.append(c.imag)
    return res


def genSolve(e):
    k = sc.sqrt(2. * (e) * m_e / (hbar ** 2))
    kw = sc.sqrt(2. * (e - V0) * m_e / (hbar ** 2))
    res = opt.fsolve(f2p_as_reals, [10., -740., -10., 25, 58.213, -127., 150., 0.], (F, k, kw))
    return reals_as_complex(res)


def checkData(V, Z):
    if len(V) != len(Z) + 1:
        raise AttributeError
    print 'AllClear'


def getK(e, V):
    return sc.sqrt(2. * (e - V) * m_e / (hbar ** 2))


def Fun(x, K, Z):
    f = []
    i = 1
    N = len(Z)
    n = 2 * N
    if N == 1:
        f.append(psi_n(1., x[0], K[0], Z[0]) - psi_n(x[n - 1], 0, K[N], Z[N - 1]))
        f.append(psiPrime_n(1., x[0], K[0], Z[0]) - psiPrime_n(x[n - 1], 0, K[N], Z[N - 1]))
    if (N > 1):
        f.append(psi_n(1., x[0], K[0], Z[0]) - psi_n(x[1], x[2], K[1], Z[0]))
        f.append(psiPrime_n(1., x[0], K[0], Z[0]) - psiPrime_n(x[1], x[2], K[1], Z[0]))
    while (i < N - 1):
        f.append(psi_n(x[i], x[i + 1], K[i], Z[i]) - psi_n(x[i + 2], x[i + 3], K[i + 1], Z[i]))
        f.append(psiPrime_n(x[i], x[i + 1], K[i], Z[i]) - psiPrime_n(x[i + 2], x[i + 3], K[i + 1], Z[i]))
        i += 1
    if (N > 1):
        f.append(psi_n(x[n - 3], x[n - 2], K[N - 1], Z[N - 1]) - psi_n(x[n - 1], 0, K[N], Z[N - 1]))
        f.append(psiPrime_n(x[n - 3], x[n - 2], K[N - 1], Z[N - 1]) - psiPrime_n(x[n - 1], 0, K[N], Z[N - 1]))
    return np.array(f)


def func_as_reals(x, F, K, Z):
    cIn = reals_as_complex(x)
    cOut = F(cIn, K, Z)
    res = []
    for c in cOut:
        res.append(c.real)
        res.append(c.imag)
    return res


def genSolveVec(K, Z):
    try:
        res = opt.fsolve(func_as_reals, np.zeros(len(Z) * 2 * 2), (Fun, K, Z))
    except RuntimeWarning:
        res = np.zeros(len(Z) * 2 * 2)
    return res


def findIndex(x, Z):
    for i in range(len(Z) - 1):
        if Z[i] <= x <= Z[i + 1]:
            return i


def getPsi(x, coeffs, K, Z):
    N = len(Z)
    n = 2 * N
    if x < Z[0]:
        return psi_n(1, coeffs[0], K[0], x)
    if Z[N - 1] <= x:
        return psi_n(coeffs[n - 1], 0, K[N], x)
    index = findIndex(x, Z)
    return psi_n(coeffs[2 * index + 1], coeffs[2 * (index + 1)], K[index + 1], x)

def modSquaredNormalized(psi_old, dx):
    area = np.trapz(psi_old, dx=dx)
    psi_new = psi_old / area
    return psi_new

def modSquared(x, coeffs, K, Z):
    N = len(Z)
    n = 2 * N
    if N == 1:
        if x < Z[0]:
            return abs(psi_n(1, coeffs[0], K[0], x)) ** 2
        if Z[N - 1] <= x:
            return abs(psi_n(coeffs[n - 1], 0, K[N], x)) ** 2
    if N > 1:
        if x < Z[0]:
            return abs(psi_n(1, coeffs[0], K[0], x)) ** 2
        if Z[N - 1] < x:
            return abs(psi_n(coeffs[n - 1], 0, K[N], x)) ** 2
        index = findIndex(x, Z)
        return abs(psi_n(coeffs[2 * index + 1], coeffs[2 * (index + 1)], K[index + 1], x)) ** 2

def potential(V,Z,x):
    N = len(Z)
    if N == 1:
        if x < Z[0]:
            return V[0]
        if Z[N - 1] <= x:
            return V[1]
    if N > 1:
        if x < Z[0]:
            return V[0]
        if Z[N - 1] < x:
            return V[len(V)-1]
        index = findIndex(x, Z)
        return V[index+1]



V0 = 1. * eV  # hartree
L = 10. * Ang  # bohr

E = np.linspace(0., 3 * eV, 1001, True)
k = sc.sqrt(2. * (E - V0) * m_e / (hbar ** 2))


# e = 0.01

Energies = [0.01, 0.02, 0.03, 0.05, 0.091]
for e in Energies:
    V = np.array([-4.*eV, -1.*eV, 1. * eV, 2. * eV, 3.*eV])
    Z = np.array([-L, 0., L, 2.*L])
    checkData(V, Z)

    pot = lambda x: potential(V,Z,x)

    K = getK(e, V)
    print K
    k = sc.sqrt(2. * (e) * m_e / (hbar ** 2))
    kw = sc.sqrt(2. * (e - V0) * m_e / (hbar ** 2))

    coeffsV = genSolveVec(K, Z)
    if np.all(coeffsV==0.) == False:
        coeffsV = reals_as_complex(coeffsV)
    # print coeffsV
        print Fun(coeffsV,K,Z)
        infinity = 50 * Ang
        dx = 2*(infinity)/10000
        X = np.linspace(-infinity, infinity, 10000)
        Y = np.array([modSquared(x, coeffsV, K, Z) for x in X])
        Y = modSquaredNormalized(Y, dx)
        print np.trapz(Y, dx=dx)
        Yv = [pot(x) for x in X]
        print len(Yv)
        print Yv[0], pot(-infinity)
        print len(X)
        # plt.figure(figsize=(13,13))
        if e==Energies[0]:
            plt.plot(X,Yv,label='V(X)', color='black', linewidth=3)
        plt.plot(X, e+Y, label='E=' + str(e)+' hartree')
        plt.axvline(0, label='Left_border', linestyle='--', color='black')
        plt.axvline(L, label='Right_border', linestyle='--', color='black')
        plt.legend(loc='best')
        plt.axhline(e, label='E='+str(e)+' hartree', color='black')
        plt.xlabel('$x,\ bohrs$',fontsize=20)
        plt.ylabel('$|\phi(x)|^{2}$',fontsize=20)
plt.show()
