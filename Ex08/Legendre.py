from sympy import Symbol, sympify, lambdify, integrate
import matplotlib.pyplot as plt
import numpy as np

x = Symbol('x')

P = []

P.append(1.)
P.append(x)

def legendre(P_n, P_1, n):
    return sympify(((2 * n + 1) * x * P_n - n * P_1) / (n + 1)).simplify().expand()

# Last polynomial I store is the 20th
for i in range(2, 21):
    P.append(legendre(P[i - 1], P[i - 2], i - 1))

# for i in range(2, 21):
#     H.append(hermite(H[i - 1], H[i - 2], i - 1))

Pf = [lambdify(x, h) for h in P]
X = np.linspace(-1.,1.,100)
# You can chose what polynomials to print and plot as you want.

desirablePolynomials = [1, 2, 3, 4, 5]
GramMatrix = []
for n, idx in zip(desirablePolynomials,range(len(desirablePolynomials))):
    GramMatrix.append([])
    for m in desirablePolynomials:
        GramMatrix[idx].append(integrate(P[n] * P[m], (x, -1., 1.)))

Y = [Pf[i](X) if i >= 1 else np.ones_like(X) for i in desirablePolynomials]

for y,i in zip(Y,desirablePolynomials):
    plt.plot(X,y,label=str(i))
    print 'P[' + str(i) + ']=' + str(P[i])
print GramMatrix
plt.legend(loc='best')
plt.show()