#!/usr/bin/python

import os, shutil, re
import numpy as np
import matplotlib.pyplot as plt

# Make arrays to store data
gnus    = []
T_res   = []
N_tiles = []

m = '/home/mhoffman/scratch/TR_'

gnu_res = [0.5,0.2,0.1,0.05,0.02]
folders = []

for i in gnu_res:
    a = m+str(i)
    folders.append(a)

#for i in folders:
#        print i

tiles = []
avgs = []

# go in each folder and obtain values of T_r

top = os.getcwd()
for i in folders:
        os.chdir(i)
        print(os.getcwd())
        gr = i.split('_')[-1]
        a = 0
        s = 0
        with open('report.txt','r') as f:    
            for line in f:
                if 'GEOM. MEAN NORM RES' in line:
                    b = line
                    Tr = float(b.split()[6])
                    gnus.append(gr)
                    T_res.append(Tr)
                    a = a+1
                    s = s+Tr
                    print('Okay, done for one tile')
        avg = s/a
        print('There are '+str(a)+' tile(s)')
        print('with an average of '+str(avg)+' geom. norm res')
        tiles.append(a)
        avgs.append(avg)
        os.chdir(top)
#print(gnus)
#print(T_res)

fig = plt.figure()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ax1 = fig.add_subplot(111)
resid = ax1.scatter(gnus, T_res, color='b', alpha=0.7, edgecolor='k', linewidth=1, s=55)
res_av = ax1.scatter(gnu_res, avgs, color='k', alpha=0.9, marker='x', s=90, linewidth=1.3)
ax1.set_ylabel("$l^2$-norm of the Normalized Residuals")

ax2 = ax1.twinx()
tileN = ax2.scatter(gnu_res,tiles, color='violet', marker='^', alpha=0.9, edgecolor='k', linewidth=1, s=90)
ax2.set_ylabel('Number of tiles',rotation=270,labelpad=14)

lims = [
    np.min([ax1.get_xlim(),ax1.get_ylim()]),
    np.max([ax1.get_xlim(),ax1.get_ylim()]),
    ]

ax1.plot(lims, lims, 'k-', alpha=0.75)

leg1 = fig.legend((resid,res_av,tileN),("$l^2$-norm (1 tile)","Avg. $l^2$-norm","Number of Tiles"), scatterpoints=1,loc="upper right", bbox_to_anchor=[0.88,0.88], shadow=True)

ax1.text(0.35, 0.37, "$l^2$-norm = Residual Threshold", ha="center", va="center", rotation=34, alpha=0.7)

ax1.set_xlabel('Residual Threshold')

ax1.set_xlim(0.0,0.6)
ax1.set_ylim(0,0.7)

fig.savefig('tileplot3.eps',format='eps',dpi=1000)

plt.show()
