#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, shutil, re
import collections
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata

path = '/data/mhoffman/NUQ/1M_grid/data/'
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

CX = []
CY = []
CZ = [] 

for i in SortedX:
    a = i.tolist()
    CX.append(a)
for i in SortedY:
    a = i.tolist()
    CY.append(a)
for i in SortedZ:
    a = i.tolist()
    CZ.append(a)

NX = np.asarray(CX)
NY = np.asarray(CY)
NZ = np.asarray(CZ)

# Now we want to break everything into smaller arrays
# size 6, as defined by the formula in the QRSM paper

def Break(M,size):
    result = []

    num = M.shape[0]/size
    nx = int(num)
    num2 = M.shape[1]/size
    ny = int(num2)
    
    Q = np.array_split(M,nx,axis=1)
    
    for i in Q:
        QN = np.array_split(i,ny)
        result.append(QN)

    return result

splitX = Break(NX,4)
splitY = Break(NY,4)
splitZ = Break(NZ,4)

#print splitX
#print splitY
print splitZ


#xi = np.linspace(min(nX),max(nX))
#yi = np.linspace(min(nY),max(nY))

#Q = griddata(nX,nY,nC,xi,yi,interp='linear')

#fig, ax = plt.subplots()
#ax.scatter(SortedX,SortedY,s=100,c=SortedZ,cmap=plt.cm.gray,vmin=min(C),vmax=max(C))

