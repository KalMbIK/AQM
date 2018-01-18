import matplotlib.pyplot as plt
import numpy as np
from math import exp, sqrt


R = []
R.append(lambda r: 2.*exp(-r))
R.append([])
R[1].append(lambda r: 1./sqrt(2)*(1-1/2.*r)*exp(-r/2))
R[1].append(lambda r: 1./sqrt(24)*r*exp(-r/2))
R.append([])
R[2].append(lambda r: 2./sqrt(27)*(1-2./3.*r+2./27.*r*r)*exp(-r/3))
R[2].append(lambda r: 8./27/sqrt(6)*r*(1-1./6.*r)*exp(-r/3))
R[2].append(lambda r: 4./81/sqrt(30)*r*r*exp(-r/3))

u = []
u.append(lambda r: r*R[0](r))
u.append([])
u.append([])
u[1].append(lambda r: r*R[1][0](r))
u[1].append(lambda r: r*R[1][1](r))
u[2].append(lambda r: r*R[2][0](r))
u[2].append(lambda r: r*R[2][1](r))
u[2].append(lambda r: r*R[2][2](r))

X = np.linspace(0, 25, 1000)

def plotter(X, funs, type):
    Y1 = [funs[0](x) for x in X]
    plt.plot(X,Y1,label=type+'10')

    for i in [1,2]:
        for Rfun, idx in zip(funs[i],range(0,len(funs[i]))):
            Y1 = [Rfun(x) for x in X]
            plt.plot(X,Y1,label=type+str(i+1)+str(idx))

plotter(X,u,'u')

plt.legend(loc='best')
plt.show()