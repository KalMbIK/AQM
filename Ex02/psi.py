import scipy.constants as sc
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
# from math import exp, sqrt
from cmath import exp,sqrt

# def densN(x,t,sk,k0,x0):
#     w_s = 2*sc.hbar/sc.m_e*(sk**2)
#     w_0 = sc.hbar/2/sc.m_e * (k0 ** 2)
#     vg = sc.hbar*k0/sc.m_e
#     d = 1+(w_s*t)**2
#     return sqrt(2./sc.pi/d)*sk*exp(-2*(sk*(x-x0-vg*t))**2/d)

def psiN(x,t,sk,k0,x0):
    w_s = 2*sc.hbar/sc.m_e*(sk**2)
    w_0 = sc.hbar/2/sc.m_e * (k0 ** 2)
    vg = sc.hbar*k0/sc.m_e
    return sqrt(2*sk/sqrt(2*sc.pi)/(1+1j*w_s*t))*exp(1j*(k0*x-w_0*t))*exp(-1/(1+1j*w_s*t)*(sk*(x-x0-vg*t))**2)


nm = 10**-9
sX = nm
la = 0.1*nm
x0 = -15*nm
k0 = 2*sc.pi/la
sk = 1./2/sX
x1 = 15*nm
k1 = -2*sc.pi/la
sk1 = 1./2/sX

# p = 2+1j
# p.conjugate()

psi1 = lambda x,t: psiN(x,t,sk,k0,x0)
psi2 = lambda x,t: psiN(x,t,sk1,k1,x1)

def dens(x,t):
    return abs(psi1(x,t))**2

def dens2(x,t):
    return abs(psi1(x,t)-psi2(x,t))**2


T = np.linspace(0., 10**-14, 1000)
X = np.linspace(-70*nm, 70*nm, 100000)
easyDensNV = np.vectorize(dens2)
for time,idx in zip(T,range(len(T))):
    # if time/10 == 0:
    if idx % 25 == 0:
        plt.plot(X, easyDensNV(X, time))
        # print integrate.quad(lambda x: easyDensNV(x, time), -50*nm,50*nm)
# plt.legend(loc='best')
# plt.plot(X, easyDensNV(X, 0.1*10**-15))
plt.show()