from math import sin, cos, sqrt, exp
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

hbar = 1.
m_e = 1.
a_0 = 1.
Hartree = 1.
eV = Hartree / 27.211386
Ang = a_0 / 0.529177210

V0 = 4. * eV  # hartree
L = 5. * Ang  # bohr

t = 2 * m_e * V0 / (hbar ** 2)

class AbstractFiniteWell(object):
    def __init__(self, V0_, L_, Nk):
        self.V0 = V0_
        self.L = L_
        self.t = 2 * m_e * V0_ / (hbar ** 2)
        # self.space = np.linspace(-sqrt(self.t), sqrt(self.t), Nk, True)
        self.space = np.linspace(0., sqrt(self.t), Nk, True)
        self.matcher = np.vectorize(self.match)
        self.k = self.solve()
        self.kappa = np.sqrt(self.t - self.k ** 2)
        self.c = self.coeff()

    def getPsiFunctions(self):
        psiV = np.vectorize(self.psi)
        fArray = lambda x: psiV(x,self.k,self.kappa,self.c)
        # fArray = lambda x: psiV(x,self.k[self.allowed:],self.kappa[self.allowed:],self.c[self.allowed:])
        return fArray

    def getEnergies(self):
        eArray = []
        for i  in range(len(self.k)):
        # for i  in range(self.allowed,len(self.k)):
            eArray.append(self.energy(self.k[i]))
        return eArray

    def energy(self,k):
        return (k**2)/2-self.V0

    def psi(self, x, k, kappa, coef):
        raise NotImplementedError

    def match(self, k):
        raise NotImplementedError

    def coeff(self):
        raise NotImplementedError

    def coeff2(self):
        raise NotImplementedError

    def solve(self):
        sol = []
        tt = self.matcher(self.space)
        A = tt[:-1]
        B = tt[1:]
        sgn = np.sign(A * B)
        for i in range(len(sgn)):
            if sgn[i] == -1.:
                ttt = opt.bisect(self.matcher, self.space[i], self.space[i + 1])
                if ttt != 0:
                    sol.append(ttt)
        return np.array(sol)

class OddFiniteWell(AbstractFiniteWell):
    def __init__(self, V0_, L_, Nk):
        AbstractFiniteWell.__init__(self, V0_, L_, Nk)
        # self.allowed = 0

    def psi(self, x, k, kappa, coef):
        if x < -L/2.:
            return -exp(kappa*x)
        if -L/2. <= x <= L/2.:
            return coef*sin(k*x)
        if L/2. < x:
            return exp(-kappa*x)

    def match(self, k):
        return k * cos(k * L / 2) + sqrt(t - k ** 2) * sin(k * L / 2)

    def coeff(self):
        return np.exp(-self.kappa * L / 2) / np.sin(self.k * L / 2)

    def coeff2(self):
        return -self.kappa/self.k*np.exp(-self.kappa * L / 2) / np.cos(self.k * L / 2)

class EvenFiniteWell(AbstractFiniteWell):
    def __init__(self, V0_, L_, Nk):
        AbstractFiniteWell.__init__(self, V0_, L_, Nk)
        # self.allowed = 0
        # for k in self.k:
        #     if k < 0:
        #         self.allowed+=1

    def psi(self, x, k, kappa, coef):
        if x < -L/2.:
            return exp(kappa*x)
        if -L/2. <= x <= L/2.:
            return coef*cos(k*x)
        if L/2. < x:
            return exp(-kappa*x)

    def match(self, k):
        return k * sin(k * L / 2) - sqrt(t - k ** 2) * cos(k * L / 2)

    def coeff(self):
        return np.exp(-self.kappa * L / 2) / np.cos(self.k * L / 2)

    def coeff2(self):
        return self.kappa/self.k*np.exp(-self.kappa * L / 2) / np.sin(self.k * L / 2)

def normalize(psi_old, dx):
    area = np.trapz(psi_old**2, dx=dx)
    psi_new = psi_old / np.sqrt(area)
    return psi_new

def getNormalizedFunctions(psis,X):
    Y = [psis(x) for x in X]
    Y = np.transpose(Y)
    Z = []
    for y in Y:
        Z.append(normalize(y,(X[len(X)-1]-X[0])/1000))
    return Z

def fromEtoK(e):
    return sqrt(2.*(e+V0))

V0 = 4. * eV  # hartree
L = 20. * Ang  # bohr
X = np.linspace(-2*L,2*L,1001,True)
even = EvenFiniteWell(V0,L,101)
odd = OddFiniteWell(V0,L,101)

psisOdd = odd.getPsiFunctions()
energiesOdd = odd.getEnergies()
psisEven = even.getPsiFunctions()
energiesEven = even.getEnergies()

Z = getNormalizedFunctions(psisOdd, X)
for y,e in zip(Z, energiesOdd):
    plt.plot(X,y, label=str(e/eV))

Z = getNormalizedFunctions(psisEven, X)
for y,e in zip(Z, energiesEven):
    plt.plot(X,y, label=str(e/eV))
    # plt.plot(X,y, label=str(fromEtoK(e)))
plt.axvline(-L/2, label='Left_border', linestyle='--', color='black')
plt.axvline(L/2, label='Right_border', linestyle='--', color='black')
plt.legend(loc='best')
plt.show()