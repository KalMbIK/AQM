from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *

nm = 10**-9
L = 1.*nm
A = sqrt(256/(63*L))

coeffs = {1:10*A/16, 3: -5*A/16, 5:A/16}

def sinN(n, z):
    return sin(pi*n*z/L)
def E(n):
    return (pi**2)*(hbar**2)*(n**2)/2/m_e/L**2

E1 = E(1)

def density(z,t):
    return (1.0/63.0/L)*(100.0*(sin(pi*z/L))**2+25.0*(sin(3.0*pi*z/L))**2+(sin(5.0*pi*z/L))**2-100.0*sin(3.0*pi*z/L)*sin(pi*z/L)*cos(8.0*E1*t/hbar)+20.0*sin(5.0*pi*z/L)*sin(pi*z/L)*cos(24.0*E1*t/hbar)-10.0*sin(3*pi*z/L)*sin(5*pi*z/L)*cos(16.0*E1*t/hbar))
# def mySol(z, t):
#     res = 0.
#     for k in coeffs.keys():
#         res+=coeffs.get(k)*sinN(k,z)*cm.exp(-1j/hbar*E(k)*t)
#     return res
#
# def mySolC(z,t):
#     return np.conjugate(mySol(z,t))
#
# def density(z,t):
#     return mySol(z,t)*mySolC(z,t)
#
densityV=np.vectorize(density)



Z = np.linspace(0,L,100)
T = np.linspace(0, nm*0.000014, 1000)
z0 = 0.5*nm
plt.plot(T,densityV(z0,T))
# for t in T:
    # plt.plot(Z,densityV(Z,t))
plt.show()