from sympy import *
from scipy.constants import *
import numpy as np
import matplotlib.pyplot as plt

z = symbols('z')
t = symbols('t')

nm = 10**-9
L = 1*nm


psi_init = sin(pi*z/L)**5
dens_init = psi_init**2
# print dens_init
res = integrate(dens_init, (z, 0., L))
A = sqrt(simplify(res))

def E (n):
    return (pi**2)*(hbar**2)*(n**2)/2/m_e/L**2

sol_coeffs = {1:10*A/16, 3: -5*A/16, 5:A/16}
# test_coeffs = {1:10, 3: -5, 5:1}

sol = 0.

for key in sol_coeffs.keys():
    sol += sol_coeffs.get(key)*sin(pi*key*z/L)*exp(-I/hbar*E(key)*t)

density = simplify(sol*conjugate(sol))

print density

solZT = lambdify((z,t), density)
sV = np.vectorize(solZT)

Z = np.linspace(0,L,100,endpoint=True)
T = np.linspace(0,nm*0.00014,10,endpoint=True)

for tt in T:
    plt.plot(Z,sV(Z,tt))
plt.show()


# for key in test_coeffs.keys():
#     sol += test_coeffs.get(key)*sin(key*z)


# print integrate(sol**2, (z, 0., L))
# print sol