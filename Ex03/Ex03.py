from math import sin, cos, sqrt, exp
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

hbar = 1.
m_e = 1.
a_0 = 1.
Hartree = 1.
eV = Hartree / 27.211386
Ang = a_0 / 0.529177210

V0 = 4. * eV  # hartree
L = 20. * Ang  # bohr

t = 2 * m_e * V0 / (hbar ** 2)


def psi(x, k, kappa, coef):
    if x < -L/2.:
        return exp(kappa*x)
    if -L/2. <= x <= L/2.:
        return coef*cos(k*x)
    if L/2. < x:
        return exp(-kappa*x)

def matchEven(k):
    return k * sin(k * L / 2) - sqrt(t - k ** 2) * cos(k * L / 2)

def matchOdd(k):
    return k * cos(k * L / 2) + sqrt(t - k ** 2) * sin(k * L / 2)

def solver(space, fun):
    sol = []
    tt = fun(space)
    A = tt[:-1]
    B = tt[1:]
    sgn = np.sign(A*B)
    for i in range(len(sgn)):
        if sgn[i] == -1.:
            ttt = opt.bisect(fun,space[i],space[i+1])
            if ttt!=0:
                sol.append(ttt)
    return np.array(sol)

def kappa(k):
    return sqrt(t-k**2)

# matchOddV = np.vectorize(matchOdd)
matchEvenV = np.vectorize(matchEven)
kappaV = np.vectorize(kappa)

K = np.linspace(-sqrt(t), sqrt(t), 100, True)

kEven = solver(K, matchEvenV)
kappaEven = kappaV(kEven)
coeffEven = np.exp(-kappaEven*L/2)/np.cos(kEven*L/2)

# print coeffEven

psiV = np.vectorize(psi)

# print psiV(L/4, kEven, kappaEven, coeffEven)

#
X = np.linspace(-2*L,2*L,1000,True)
for i in range(len(kappaEven)):
    Y = psiV(X, kEven[i], kappaEven[i], coeffEven[i])
# print Y
    plt.plot(X,Y,label=str(kEven[i]))
plt.axvline(-L/2)
plt.axvline(L/2)
plt.legend(loc='best')
plt.show()








# kOdd = solver(K, matchOddV)
# kappaOdd = kappaV(kOdd)
# coeffOdd = np.exp(-kappaOdd*L/2)/np.sin(kOdd*L/2)



# psiEven = lambda x: psi(x,kEven[0],kappaEven[0],coeffEven[0])

# psiEvenV = np.vectorize(psiEven)
# Y = np.transpose(Y)
# print Y[0]
# print kEven
# print len(kEven)
# print kOdd
# print len(kOdd)
# eK = match_evenV(K)
# oK = match_oddV(K)
# plt.plot(K, eK)
# plt.plot(K, oK)
# plt.show()
# print V0
# print L
