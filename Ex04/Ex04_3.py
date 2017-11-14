import numpy as np
from math import cos, sin, exp, sqrt, pi
import scipy.special as sp
import matplotlib.pyplot as plt

def cmpVV(x1, x2):
    abs = np.abs(x1-x2)
    return np.max(abs)

def Vparametrized(x, epsilon_):
    return epsilon_*x

def kSquaredParametrized(x, E_, V):
    return 2*(E_-V(x))

def AiMinusInf(x):
    x = sqrt(x)
    return cos(2./3.*x**3+pi/4.)/sqrt(pi*x)

def BiMinusInf(x):
    x = sqrt(x)
    return sin(2./3.*x**3-pi/4.)/sqrt(pi*x)

def AiInf(x):
    x = sqrt(x)
    return exp(-2./3.*x**3)/2./sqrt(pi*x)

def BiInf(x):
    x = sqrt(x)
    return exp(2./3.*x**3)/sqrt(pi*x)

def NumerovIterations(gr, assimptotics):
    h = gr[1]-gr[0]
    K = kSq(gr)
    K *= (h**2)/12.
    u = np.zeros(len(gr))
    u[0] = assimptotics(gr[0])
    u[1] = assimptotics(gr[1])
    # u[2:] = (2.*(1.-5*K[1:-1])*u[1:-1]-(1.+K[0:-2])*u[0:-2])/(1.+K[2:])
    for i in range(1,len(gr)-1):
        u[i+1] = (2.*(1.-5*K[i])*u[i]-(1.+K[i-1])*u[i-1])/(1.+K[i+1])
    return u

hbar = 1.
m_e = 1.
a_0 = 1.
e = 1. # Charge of an electron
Hartree = 1.
HartreeBohr = 1. # Electric field in SI
eV = Hartree / 27.211386
Ang = a_0 / 0.529177210
VAng = HartreeBohr / 51.4220652 # Electric field in AU

E = 4. * eV  # hartree
L = 5. * Ang  # bohr
Epslon = 1. * VAng

alpha = (2*Epslon)**(1/3.)
x0 = E/e/Epslon
dzita = lambda x: alpha * (x - x0)

V = lambda x: Vparametrized(x, Epslon)
kSq = lambda x: kSquaredParametrized(x, E, V)

Ai = lambda x: sp.airy(dzita(x))[0]
Bi = lambda x: sp.airy(dzita(x))[2]

a = -30.*Ang
b = x0 + 15.*Ang
delta = b-a
epsMachine = 1.11e-16
N = [100*2**x + 1 for x in range(0, 15)]
print N
H = [delta/(n-1) for n in N]
# EPS = [epsMachine/h for h in H]
errors = []
# print N
f = Ai
for n in N:
    X = np.linspace(a, b, n, True)
    # Y = NumerovIterations(X,f)
    Y = NumerovIterations(X[::-1],f)
    Y = Y[::-1]
    exactY = f(X)
    errors.append(cmpVV(Y,exactY))

plt.loglog(H,errors)
# plt.loglog(H,EPS)
plt.show()
# sp.airy()

# plt.plot(X,Y)
# plt.axvline(x0)
# plt.show()



