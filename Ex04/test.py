import scipy.special as sp
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

def normalize(psi_old, dx):
    area = np.trapz(psi_old**2, dx=dx)
    sq = np.sqrt(area)
    psi_new = psi_old / sq
    return psi_new, sq

Ai = lambda x: sp.airy(x)[0]
alpha = 0.198
alpha2 = 0.249
dz1 = -2.338
dz2 = -4.088
dz3 = -5.521
X = np.linspace(0, 200, 10000)
Y = np.array([Ai(alpha*x+dz2) for x in X])
Y, sq = normalize(Y,X[1]-X[0])
print sq
print np.trapz(Y**2, dx=(X[1]-X[0]))
print sq*2
plt.plot(X,Y)
plt.show()