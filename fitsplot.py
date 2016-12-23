#!/usr/bin/python

# Make a surface plot of the model we got from R

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig=plt.figure()
ax = fig.gca(projection='3d')

X = np.linspace(0.01,0.1,50)
Y = np.linspace(0.3,0.9,50)

Z = 0.5490 + ( -0.6414 * X ) + ( 0.1296 * Y ) + ( 2.0293 * X * X ) + ( -0.1975 * Y * Y ) + ( 0.4279 * X * Y )

X, Y = np.meshgrid(X, Y)

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap = cm.RdPu, vmax=max(Z), vmin=min(Z), linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

plt.show()
