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

# This is my data as np arrays

nX = np.asarray(X)
nY = np.asarray(Y)
nC = np.asarray(C)

def REORDER(nX,nY,nC):
    # Reordering everything based on X

    Xord = np.argsort(nX)
    XX = nX[Xord]
    XY = nY[Xord]
    XZ = nC[Xord]
    
    # Find the places where the values of
    # X are different and then split things
    # up accordingly 
    
    Q = np.where(np.diff(XX[:]))[0]+1
    FirstX = XX[0:Q[0]]
    FirstY = XY[0:Q[0]]
    FirstZ = XZ[0:Q[0]]

    Xdiv = []
    Ydiv = []
    Zdiv = []
    Xdiv.append(FirstX)
    Ydiv.append(FirstY)
    Zdiv.append(FirstZ)
    
    for i in range(len(Q)-1):
        a = Q[i]+1
        b = Q[i+1]
        NewX = XX[a:b]
        NewY = XY[a:b]
        NewZ = XZ[a:b]
        Xdiv.append(NewX)
        Ydiv.append(NewY)
        Zdiv.append(NewZ)

    LastX = XX[Q[len(Q)-1]:len(XX)-1]    
    LastY = XY[Q[len(Q)-1]:len(XX)-1]
    LastZ = XZ[Q[len(Q)-1]:len(XX)-1]
    Xdiv.append(LastX)
    Ydiv.append(LastY)
    Zdiv.append(LastZ)

    SortedX = []
    SortedY = []
    SortedZ = []

    # Now find order according to Y data
    # For each segment of X data 

    for i in range(len(Ydiv)):
            
        YOrd = np.argsort(Ydiv[i])
        USX = Xdiv[i]
        USY = Ydiv[i]
        USZ = Zdiv[i]
        SorX = USX[YOrd]
        SorY = USY[YOrd]
        SorZ = USZ[YOrd]
        SortedX.append(SorX)
        SortedY.append(SorY)
        SortedZ.append(SorZ)

    # Now Z (stellar mass) should be organized into
    # a list of arrays in the correct order
    # ascending X order, ascending Y order

    result = np.asarray(SortedZ)

    return result

NZ = REORDER(nX,nY,nC)
TNZ = NZ.transpose()

locs = []
for i in range(len(NZ)):
    row = i
    xcol = argrelmax(NZ[i],order=1)
    ncol = argrelmin(NZ[i],order=1)
    txcol = argrelmax(TNZ[i],order=1)
    tncol = argrelmin(TNZ[i],order=1)
    if (len(xcol[0]) == 0):
        pass
    else:
       for j in xcol[0]:
           r = row
           n = int(j)
           a = (r,n)
           locs.append(a)
    if (len(ncol[0])== 0):
        pass
    else:
        for j in ncol[0]:
            r = row
            n = int(j)
            a = (r,n)
            locs.append(a)
    if (len(tncol[0])==0):
        pass
    else:
        for j in tncol[0]:
            r = row
            n = int(j)
            a = (n, r)
            locs.append(a)
    if (len(txcol[0])==0):
        pass
    else:
        for j in txcol[0]:
            r = row
            n = int(j)
            a = (n,r)
            locs.append(a)
            
print locs
#print( argrelmax(TNZ,order=1,axis=1) )

    
xi = np.linspace(min(nX),max(nX))
yi = np.linspace(min(nY),max(nY))

Q = griddata(nX,nY,nC,xi,yi,interp='linear')
#XX, YY = np.meshgrid(xi,yi)

#ix, iy = np.where(ndimage.maximum_filter(ZZ, size=(10,10), mode='constant') == ZZ)

#neighborhood = generate_binary_structure(2,2)
#LMX = maximum_filter(ZZ, footprint=neighborhood)==ZZ
#bkgd = (ZZ==0)
#er_bkgd = binary_erosion(bkgd, structure=neighborhood, border_value=1) 
#maxs = LMX - er_bkgd

#fig, ax = plt.subplots()
#ax.scatter(X,Y,s=100,c=C,cmap=plt.cm.gray,vmin=min(C),vmax=max(C))

plt.show()
