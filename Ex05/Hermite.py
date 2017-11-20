from sympy import Symbol, sympify, lambdify
import matplotlib.pyplot as plt
import numpy as np

x = Symbol('x')

H = []

H.append(1.)
H.append(2 * x)

def hermite(H_n, H_1, n):
    return sympify((2*x * H_n - 2*n*H_1)).simplify().expand()

for i in range(2, 10):
    H.append(hermite(H[i - 1], H[i - 2], i - 1))

Hf = [lambdify(x, h) for h in H]
X = np.linspace(-2.,2.,100)
Y = [f(X) for f in Hf]
plt.axhline(1.0, label=str(0), color='black')

for i in range(1,8):
    plt.plot(X,Y[i],label=str(i))
plt.legend(loc='best')
plt.show()