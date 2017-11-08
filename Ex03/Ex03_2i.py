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

# def psi_n(t_n, r_n, k_, x_):
#     return t_n*np.exp(1j * k_ * x_) + r_n*np.exp(-1j*k_*x_)

# def psiPrime_n(t_n, r_n, k_, x_):
#     return 1j*k_*(t_n*np.exp(1j * k_ * x_) - r_n*np.exp(-1j*k_*x_))

# psi = lambda t, r, x: psi_n(t,r,k,x)
# psiPrime = lambda t, r, x: psiPrime_n(t,r,k,x)
# psi2 = lambda t, r, x: psi_n(t,r,kw,x)
# psi2Prime = lambda t, r, x: psiPrime_n(t,r,kw,x)

def genSolve(e):
    k = sc.sqrt(2.*(e)*m_e/(hbar**2))
    kw = sc.sqrt(2.*(e-V0)*m_e/(hbar**2))
    # print psiPrime(0,0,0)
    def F(x):
        return [psi_n(1.,x[0],k,-L/2.)-psi_n(x[1],x[2],kw,-L/2.),
                psiPrime_n(1.,x[0],k,-L/2.)-psiPrime_n(x[1],x[2],kw,-L/2.),
                psi_n(x[1],x[2],kw,L/2.)-psi_n(x[3],0,k,L/2.),
                psiPrime_n(x[1],x[2],kw,L/2.)-psiPrime_n(x[3],0,k,L/2.)]
    def func_as_reals(x):
        r1, c1, r2, c2, r3, c3, r4, c4 = x
        a1, a2, a3, a4 = F([complex(r1, c1), complex(r2, c2), complex(r3, c3), complex(r4, c4)])
        return [a1.real, a1.imag, a2.real, a2.imag, a3.real, a3.imag, a4.real, a4.imag]
    res = opt.fsolve(func_as_reals,[10.,-740.,-10.,25,58.213,-127.,150.,0.])
    # print func_as_reals(res)
    return res

# psi1 = lambda x: psi_n(1.,-0.5,sc.sqrt(2.*(1.2*eV-V0)*m_e/(hbar**2)),x)
# X = np.linspace(-2*L,-L/2,1000)
# Y = [psi1(x) for x in X]
# plt.plot(X,Y)
# plt.show()
e = 1.4*eV
genSolve(e)
k = sc.sqrt(2.*(e)*m_e/(hbar**2))
kw = sc.sqrt(2.*(e-V0)*m_e/(hbar**2))
def eT3(e):
    arr = genSolve(e)
    return sqrt(arr[6]**2+arr[7]**2)
Y = [eT3(e) for e in E]
plt.plot(E,Y, label='T=T(e)')
plt.axvline(1.*eV, label='e=V0', linestyle='--', color='black')
plt.axvline(0.05, label='e=0.05 (1st_Max)', linestyle='-', color='red')
plt.axvline(0.091, label='e=0.091 (2nd_Max)', linestyle='-', color='red')
plt.legend(loc='best')
plt.xlabel('$E,\ hartree$',fontsize=20)
plt.ylabel('$|T(E)|^{2}$',fontsize=20)
plt.show()