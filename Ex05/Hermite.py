from sympy import Symbol, sympify, lambdify
import matplotlib.pyplot as plt
import numpy as np

x = Symbol('x')

H = []

H.append(1.)
H.append(2 * x)

def hermite(H_n, H_1, n):
    return sympify((2*x * H_n - 2*n*H_1)).simplify().expand()
# Last polynomial I store is the 20th
for i in range(2, 21):
    H.append(hermite(H[i - 1], H[i - 2], i - 1))

Hf = [lambdify(x, h) for h in H]
X = np.linspace(-3.,3.,100)
# You can chose what polynomials to print and plot as you want.
desirablePolynomials = [4, 8, 2]
Y = [Hf[i](X) if i >= 1 else np.ones_like(X) for i in desirablePolynomials]

for y,i in zip(Y,desirablePolynomials):
    plt.plot(X,y,label=str(i))
    print 'H[' + str(i) + ']=' + str(H[i])
plt.legend(loc='best')
plt.show()