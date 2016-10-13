#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, shutil, re
import collections
import datetime
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata

#topdir = raw_input('What directory should we go to? ')

topdirs = ['/data/mhoffman/NUQ/1M_grid','/data/mhoffman/NUQ/1M_grid_highres']

# Read in the .out files, obtained from readdata.py
# can change extension to match whatever textfile your
# data is stored in
dd = []

for i in topdirs:
    data = []
    
    os.chdir(i)
    h = os.getcwd()
    os.chdir(h+'/data/')
    th = os.getcwd()
    
    ofs = []
    for file in os.listdir(th):
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
    #for key, value in sorted(keys.iteritems(), key=lambda (k,v): (v,k)):
    #    print "%s: %s" % (key, value)
            
    kFC = keys['FailM.out']
    kFX = keys['FailR.out']
    kFY = keys['FailB.out']
    FC = data[kFC]
    FX = data[kFX]
    FY = data[kFY]
    nFX = np.asarray(FX)
    nFY = np.asarray(FY)
    nFC = np.asarray(FC)
    
    kX = keys['Reims.out']
    kY = keys['Block.out']
    X = data[kX]
    Y = data[kY]
    nX = np.asarray(X)
    nY = np.asarray(Y)
    plott = 'scat'
    kC = keys['StarM.out']
    C = data[kC]
    nC = np.asarray(C)

    dd.append(nX)
    dd.append(nY)
    dd.append(nC)
    dd.append(nFX)
    dd.append(nFY)
    
# Set plot font to the Computer Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# fig, axes = plt.subplots(nrows=1, ncols=2)
# plot1 = ax.scatter(data[0][0],data[0][1],s=50,c=data[0][2],cmap=plt.cm.gnuplot,vmin=(max(C)),vmax=(min(C)),linewidth=1)
# plot2 = ax.scatter(data[1][0],data[1][1],s=50, c='g', vmin=min(C), vmax=max(C), edgecolor='limegreen', linewidth=2, label='Not converged')
# plot3 = ax.scatter(data[2][0],data[2][1],s=50,c=data[2][2],cmap=plt.cm.gnuplot,vmin=(max(C)),vmax=(min(C)),linewidth=1)

cmax = max(max(dd[2]),max(dd[7]))
cmin = min(min(dd[2]),min(dd[7]))

fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
p1 = ax1.scatter(dd[0],dd[1],c=dd[2],cmap=plt.cm.gnuplot,vmin=cmin,vmax=cmax,linewidth=1)
p2 = ax1.scatter(dd[3],dd[4],c='g',edgecolor='limegreen',linewidth=2,label='Not converged')
p3 = ax2.scatter(dd[5],dd[6],c=dd[7],cmap=plt.cm.gnuplot,vmin=cmin,vmax=cmax,linewidth=1)
p4 = ax2.scatter(dd[8],dd[9],c='g',edgecolor='limegreen',linewidth=2)

fig.subplots_adjust(hspace=0)

cax = fig.add_axes([0.9,0.1,0.03,0.8])
cbar = fig.colorbar(p1,cax=cax)

handles = [p2]
labels = [h.get_label() for h in handles]

#box = ax2.get_position()
#ax2.set_position([box.x0, box.y0, box.width * 0.7, box.height])
legend = fig.legend(handles=handles, labels = labels, fancybox=True)

ax1.set_xlim([0.28,0.92])
ax1.set_ylim([-0.01,0.12])

ax2.set_xlabel('$\eta_{R}$', fontsize='large')
#fig.set_ylabel('&\eta_{B}$')
fig.text(0.04,0.5,'$\eta_{B}$', va='center', rotation='vertical',fontsize='large')

fig.savefig('twoplot.eps',format='eps',dpi=1000)

plt.show()

