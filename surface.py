#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, shutil, re
import collections
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata
from scipy import ndimage
from scipy.ndimage.filters import maximum_filter, minimum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.signal import argrelmax, argrelmin

path = '/data/mhoffman/NUQ/fine_grid/tmps'
os.chdir(path)

ofs = []
for file in os.listdir(os.getcwd()):
    if file.endswith('.out'):
        nf = file
        ofs.append(nf)

data = []
for file in ofs:
    with open(file) as f:
        dat = []
        for line in f:
            line = float(line)
            dat.append(line)
    data.append(dat)

nof = len(ofs)
nums = range(nof)
keys = dict(zip(ofs,nums))

for key, value in sorted(keys.iteritems(), key=lambda (k,v): (v,k)):
    print "%s: %s" % (key, value)


Xax = 'Block'
Yax = 'Reims'
Cols = 'StarM'

kX = keys[Xax+'.out']
kY = keys[Yax+'.out']
kC = keys[Cols+'.out']
X = data[kX]
Y = data[kY]
C = data[kC]
nX = np.asarray(X)
nY = np.asarray(Y)
nC = np.asarray(C)

xi = np.linspace(min(nX),max(nX))
yi = np.linspace(min(nY),max(nY))

Q = griddata(nX,nY,nC,xi,yi,interp='linear')
XX, YY = np.meshgrid(xi,yi)










#ix, iy = np.where(ndimage.maximum_filter(ZZ, size=(10,10), mode='constant') == ZZ)

#neighborhood = generate_binary_structure(2,2)
#LMX = maximum_filter(ZZ, footprint=neighborhood)==ZZ
#bkgd = (ZZ==0)
#er_bkgd = binary_erosion(bkgd, structure=neighborhood, border_value=1) 
#maxs = LMX - er_bkgd

fig, ax = plt.subplots()
ax.scatter(X,Y,s=100,c=C,cmap=plt.cm.gray,vmin=min(C),vmax=max(C))

print nX
print xi


plt.show()
