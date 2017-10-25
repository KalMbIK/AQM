import numpy as np
import matplotlib.pyplot as plt
import math as m

def eigState(x, L, n):
    return m.sqrt(2./L)*m.sin(m.pi*n/L*x)

eigStateV = np.vectorize(eigState)

nm = 10**-9
L=1.*nm
h=0.01*nm
x = np.arange(0, L+h, h)
N = 3
states = range(1,N)

Y = [eigStateV(x, L, n) for n in states]
for y,n in zip(Y,states):
    plt.plot(x, y, label='State='+str(n))
plt.legend(loc='lower left')
plt.xlim(0,L)
plt.grid()
# plt.xlabel('m')
plt._show()