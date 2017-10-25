from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from numpy import exp,arange,vectorize
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from math import pi, sin, sqrt

def eigenState(x, y, Lx, Ly, nx, ny):
    return sqrt(2*2/Lx/Ly)*sin(pi*nx/Lx*x)*sin(pi*ny/Ly*y)

eigenStateV = vectorize(eigenState)

# the function that I'm going to plot

nm = 10**-9
Lx = 1.*nm
Ly = 2.*nm
h = nm*10**-2
x = arange(0.0, Lx+h, h)
y = arange(0.0, Ly+h, h)
X, Y = meshgrid(x, y)  # grid of point
Z = eigenStateV(X, Y, Lx, Ly, 1, 1)  # evaluation of the function on the grid

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, 
                      cmap=cm.rainbow,linewidth=0, antialiased=False)

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# ax.set_xlabel('m')
# ax.set_ylabel('m')
ax.set_xlim(0, Lx)
ax.set_ylim(0, Ly)

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()