import numpy as np
import matplotlib.pyplot as plt

L = 1.
# endpoint=True means that the last point is also included
x = np.linspace(0, L, 100, endpoint=True)

A = np.sqrt(2.0/L)
psi = A*np.sin(np.pi*x/L)
psi2 = A*np.sin(2*np.pi*x/L)

plt.plot(x,psi, label='State1')
plt.plot(x,psi2, label='State2')
plt.legend(loc='best')
plt.show()
