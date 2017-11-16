import scipy.special as sp
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def normalize(psi_old, dx):
    area = np.trapz(psi_old**2, dx=dx)
    sq = np.sqrt(area)
    psi_new = psi_old / sq
    return psi_new, sq

Ai = lambda x: sp.airy(x)[0]
Bi = lambda x: sp.airy(x)[2]
# alpha = 0.198
# alpha2 = 0.249
# dz1 = -2.338
# dz2 = -4.088
# dz3 = -5.521
# X = np.linspace(0, 200, 10000)
# Y = np.array([Ai(alpha*x+dz2) for x in X])
# Y, sq = normalize(Y,X[1]-X[0])
# print sq
# print np.trapz(Y**2, dx=(X[1]-X[0]))
# print sq*2
# plt.plot(X,Y)
# plt.show()

Ang = 1. / 0.529177210
VAng = 1 / 51.4220652 # Electric field in AU
L = 12.*Ang
epsilon = (10**-4)*(500 * VAng)
alpha = (2*epsilon)**(1./3)
arg = lambda z, dz: alpha*(z+L/2.)+dz
ci = lambda dz: -Bi(dz)/Ai(dz)
# print epsilon/alpha*2.58
# print alpha*L
# print (12.52/alpha-L/2.)*epsilon
# alpha = 0.198
# alpha2 = 0.249
dz1 = -2.58
dz2 = -6.36
dz3 = -12.518
dzx = dz1
X = np.linspace(-L/2., L/2., 10000)
Y = np.array([ci(dzx)*Ai(arg(x,dzx))+Bi(arg(x,dzx)) for x in X])
Y1 = np.array([Ai(arg(x,dzx))+1./ci(dzx)*Bi(arg(x,dzx)) for x in X])
Y, sq = normalize(Y,X[1]-X[0])
Y1, sq2 = normalize(Y1,X[1]-X[0])
print sq
print ci(dzx)*sq
print alpha
print arg(0,dzx)
print np.trapz(Y**2, dx=(X[1]-X[0]))
print "state1="+str(pi*pi*(3**2)/2/L**2)
# print sq*2
plt.plot(X,Y)
plt.plot(X,Y1)
plt.show()