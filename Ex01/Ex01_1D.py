import numpy as np
import matplotlib.pyplot as plt
import math as m

def eigState(x, L, n):
    return m.sqrt(2./L)*m.sin(m.pi*n/L*x)

eigStateV = np.vectorize(eigState)

L=1.
x = np.arange(0, L, 0.01)
N = 4

Y = [eigStateV(x, L, n) for n in range(1,N)]
for y,n in zip(Y,range(1,N)):
    plt.plot(x, y, label='State='+str(n))
plt.legend(loc='best')
plt._show()