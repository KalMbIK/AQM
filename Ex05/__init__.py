import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Set up a grid

N = 1000
x, dx = np.linspace(-2.,2.,N,endpoint=True,retstep=True)


# Set up a potential

V = 50.*x**2

# Set up a figure

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)
plt.xlim(x[0],x[-1])
plt.ylim(-5.,60.)
lines = plt.plot(x, V, '-k',[],[],'--m',[],'-b',linewidth=2)

EnergyAxes = plt.axes([0.1,0.07,0.8,0.03])
EnergySlider = Slider(EnergyAxes, '$E$', 0, 50, valinit=0, valfmt='%0.1f')

def normalize(psi_old, dx):
    area = np.trapz(psi_old[0:N/2]**2, dx=dx)
    psi_new = psi_old / np.sqrt(area)
    return psi_new

def EnergyChangeHandler(val):
    E = round(EnergySlider.val,1)
    kk = 2*(E-V)

    psi = np.zeros(N)
    psi[0] = 1.0e-8
    psi[1] = 1.1e-8
    # Numerov integration
    for i in range(1,N-1):
        psi[i+1] = (2.*(1.-5/12.*dx*dx*kk[i])*psi[i]-(1.+dx*dx/12.*kk[i-1])*psi[i-1])/(1.+dx*dx/12.*kk[i+1])

    psi = normalize(psi,dx)

    # Update plots
    lines[1].set_data([x[0],x[-1]],[E,E])
    lines[2].set_data(x,5*psi+E)
    fig.canvas.draw

EnergySlider.on_changed(EnergyChangeHandler)
EnergyChangeHandler(0)

plt.show()