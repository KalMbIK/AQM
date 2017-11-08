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

V0 = 1. * eV  # hartree
L = 10. * Ang  # bohr

E = np.linspace(0., 3 * eV, 1001, True)
k = sc.sqrt(2.*(E-V0)*m_e/(hbar**2))

t = 2 * m_e * V0 / (hbar ** 2)

def psi_n(t_n, r_n, k_, x_):
    return t_n*cm.exp(1j * k_ * x_) + r_n*cm.exp(-1j*k_*x_)

def psiPrime_n(t_n, r_n, k_, x_):
    return 1j*k_*(t_n*cm.exp(1j * k_ * x_) - r_n*cm.exp(-1j*k_*x_))

def func_as_reals(x, F):
    cIn = np.transpose([x[:-1:2],x[1::2]])
    cIn = [complex(c[0],c[1]) for c in cIn]
    cOut = F(cIn)
    res = []
    for c in cOut:
        res.append(c.real)
        res.append(c.imag)
    return res

def genSolve(e):
    k = sc.sqrt(2.*(e)*m_e/(hbar**2))
    kw = sc.sqrt(2.*(e-V0)*m_e/(hbar**2))
    # print psiPrime(0,0,0)
    def F(x):
        return [psi_n(1.,x[0],k,-L/2.)-psi_n(x[1],x[2],kw,-L/2.),
                psiPrime_n(1.,x[0],k,-L/2.)-psiPrime_n(x[1],x[2],kw,-L/2.),
                psi_n(x[1],x[2],kw,L/2.)-psi_n(x[3],0,k,L/2.),
                psiPrime_n(x[1],x[2],kw,L/2.)-psiPrime_n(x[3],0,k,L/2.)]
    res = opt.fsolve(func_as_reals,[10.,-740.,-10.,25,58.213,-127.,150.,0.],F)
    return res
e1 = 0.05
e = e2= 0.091
# e = e1
print pi*sqrt(2/e1)/Ang
print pi*sqrt(2/e2)/Ang
k = sc.sqrt(2.*(e)*m_e/(hbar**2))
kw = sc.sqrt(2.*(e-V0)*m_e/(hbar**2))
coeffs = genSolve(e)
r1 = complex(coeffs[0],coeffs[1])
t2 = complex(coeffs[2],coeffs[3])
r2 = complex(coeffs[4],coeffs[5])
t3 = complex(coeffs[6],coeffs[7])

psi1 = lambda x: psi_n(1,r1,k,x)
psi2 = lambda x: psi_n(t2,r2,kw,x)
psi3 = lambda x: psi_n(t3,0,k,x)

def modSquared(x):
    if x < -L/2.:
        return abs(psi1(x))**2
    if -L/2. <= x < L/2.:
        return abs(psi2(x))**2
    if L/2. <= x:
        return abs(psi3(x))**2
alpha = 1.
X = np.linspace(alpha*(-L/2-1*Ang),alpha*(L/2+1*Ang),1000)
Y = [modSquared(x) for x in X]
plt.plot(X,Y, label='E='+str(e/eV)+' eV')
plt.axvline(-L/2, label='Left_border', linestyle='--', color='black')
plt.axvline(L/2, label='Right_border', linestyle='--', color='black')
plt.legend(loc='best')
plt.xlabel('$x,\ bohrs$',fontsize=20)
plt.ylabel('$|\phi(x)|^{2}$',fontsize=20)
plt.show()