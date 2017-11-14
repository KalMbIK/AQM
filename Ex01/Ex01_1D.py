import numpy as np
import matplotlib.pyplot as plt
import math as m

def eigState(x, L, n):
    return m.sqrt(2./L)*m.sin(m.pi*n/L*x)

eigStateV = np.vectorize(eigState)

def normalize(psi_old, dx):
    area = np.trapz(psi_old**2, dx=dx)
    sq = np.sqrt(area)
    psi_new = psi_old / sq
    return psi_new, sq

nm = 10**-9
# L=1.*nm
Ang = 1. / 0.529177210
L = 12.*Ang
h=0.01*nm
# x = np.arange(0, L+h, h)
x = np.linspace(0, L, 10000)
N = 3
states = range(1,N)

Y = [eigStateV(x, L, n) for n in states]
for y,n in zip(Y,states):
    # psi, sq = normalize(y, dx=(x[1] - x[0]))
    plt.plot(x, y, label='State='+str(n))
plt.legend(loc='lower left')
plt.xlim(0,L)
plt.grid()
# plt.xlabel('m')
plt._show()