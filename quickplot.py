#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, shutil, re
import collections
import datetime
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.mlab import griddata

#wheretogo = raw_input('What directory should we go to? ')
wheretogo = '/data/mhoffman/NUQ/fine_grid/'
os.chdir(wheretogo+'/tmps/')
here = os.getcwd()

# Read in the .out files, aka the text files produced from the earlier plotting code
ofs = []
for file in os.listdir(here):
    if file.endswith('.out'):
        nf = file
        ofs.append(nf)
        
# append data to a big list of data  
data = []
for file in ofs:
    with open(file) as f:
        dat = []
        for line in f:
            line = float(line)
            dat.append(line)
    data.append(dat)                                                                                            

# key which tells you which index corresponds to which dataset 
nof = len(ofs)
nums = range(nof)
keys = dict(zip(ofs,nums))
for key, value in sorted(keys.iteritems(), key=lambda (k,v): (v,k)):
    print "%s: %s" % (key, value)

#Xax = raw_input('Which quantity would you like to be your xaxis? ')
#Yax = raw_input('Which quantity would you like to be your yaxis? ')
#Q = raw_input('Is this a scatter plot (s) or line plot (l) or surface plot (su)? [s/l/su] ') 
#Titl = raw_input('Fancy titles? [y/n]')
#Q2 = raw_input('Would you like to see failed points [y/n]? ')
Xax = 'Block'
Yax = 'Reims'
Q = 'su'
Titl = 'y'
Q2 = 'n'

if Q2 == 'y':
    #FXax = raw_input('Failed x axis? ')
    #FYax = raw_input('Failed y axis? ')
    #FCol = raw_input('Failed color? ')
    FXax = 'failB'
    FYax = 'failR'
    FCol = 'failM'
    kFC = keys[FCol+'.out']
    kFX = keys[FXax+'.out']
    kFY = keys[FYax+'.out']
    FC = data[kFC]
    FX = data[kFX]
    FY = data[kFY]
    nFX = np.asarray(FX)
    nFY = np.asarray(FY)
    nFC = np.asarray(FC)
    
kX = keys[Xax+'.out']
kY = keys[Yax+'.out']
X = data[kX]
Y = data[kY]
nX = np.asarray(X)
nY = np.asarray(Y)

if Q == 's':
    plott = 'scat'
    #Cols = raw_input('Which quantity would you like to be your colormap?')
    Cols = 'StarM'
    kC = keys[Cols+'.out']
    C = data[kC]
    fig, ax = plt.subplots()
    plot = ax.scatter(X,Y,s=50,c=C,cmap=plt.cm.gnuplot,vmin=(max(C)),vmax=(min(C)),linewidth=1)
    if Q2 == 'y':
        plot2 = ax.scatter(FX,FY,s=50, c='g', vmin=min(C), vmax=max(C), edgecolor='limegreen', linewidth=2, label='Not converged')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    cbar = fig.colorbar(plot)
    cbar.set_label(Cols)
    legend = ax.legend(scatterpoints=1, fancybox = 'True', loc='upper left', bbox_to_anchor=(0.8,-0.035), fontsize='small')

if Q == 'l':
    plott = 'line'
    print('In Construction...')
    exit()

if Q == 'su':
    plott = 'surf'
    Cols = raw_input('Which quantity would you like to be your colormap?')
    kC = keys[Cols+'.out']
    C = data[kC]
    nC = np.asarray(C)
    xi = np.linspace(min(nX),max(nX))
    yi = np.linspace(min(nY),max(nY))
    Z = griddata(nX,nY,nC,xi,yi,interp='linear')
    XX, YY = np.meshgrid(xi,yi)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=plt.cm.gnuplot, vmax=np.max(Z), vmin=np.min(Z))
    cbar = fig.colorbar(surf)
    
if Titl == 'y':
    print('In Construction...')
    Xname = raw_input('What should the X axis be called?')
    Yname = raw_input('What should the Y axis be called?')
    Gtitle = raw_input('What should the title be?')
    Colba = raw_input('What is the colorbar label?')
    cbar.set_label(Colba)
    ax.set_title(Gtitle)
    ax.set_xlabel(Xname)
    ax.set_ylabel(Yname)
    xmin = min(X) - 0.01
    xmax = max(X) + 0.01
    ymin = min(Y) - 0.01
    ymax = max(Y) + 0.01
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
if Titl == 'n':
    ax.set_xlabel(Xax)
    ax.set_ylabel(Yax)
    ax.set_title(Yax+' v. '+Xax)
    if Q == 'su':
        ax.set_zlabel(Cols)
    if Q == 's':
        xmin = min(X) - 0.01
        xmax = max(X) + 0.01
        ymin = min(Y) - 0.01
        ymax = max(Y) + 0.01
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)                    
#else:
#    print('Errr... not sure what you want for titles, then...?')
#    exit()

ax.title.set_fontsize(20)
ax.yaxis.label.set_fontsize(18)
ax.xaxis.label.set_fontsize(18)
ax.grid()

now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

plt.savefig(Yax+'v'+Xax+'_'+plott+month+day+'_'+hour+minute+'.png')
#plt.savefig('test.png')
plt.show()

