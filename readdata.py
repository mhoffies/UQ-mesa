#!/usr/bin/python

import os, shutil, re
import mesa as ms
import numpy as np
import matplotlib.pyplot as plt
import datetime

def ReadInls(inlist,value):
        inl = open(inlist,'r')
        for line in inl:
                if str(value)in line:
                        a = line
                        val = a.split()[-1]
                        v = re.sub('d','e',val)
                        v = float(v)
                        return v

# Make arrays to store Final Stellar Mass, Reimers, and Blocker
masses = []
reims = []
block = []
lumi = []
failm = []
failr = []
failb = []

#paths = ['/data/mhoffman/NUQ/1M_grid_highres']
paths = ['/data/mhoffman/NUQ/CD_grid_highres_2']

for i in paths:
    os.chdir(i)

    # Get the final masses from the folders
    #masses = []
    for file in os.listdir(i):
    
        # Check if it's a directory, if it's a file, we'll ignore it
        if os.path.isdir(i+'/'+file) == True:

            if os.path.isdir(i+'/'+file) == i+'/data':
                    print(file+' is the data folder, ignoring...')
            else:
                    os.chdir(i+'/'+file)
                    currd = os.getcwd()
                    print('Now in directory... '+os.getcwd())
                    os.system('rm ./LOGS/history.datasa')
                    
                    s = ms.history_data()
                    lum = s.get('log_L')
                    mass = s.get('star_mass')
                    
                    fl = lum[ len(lum)-1 ]
                    fm = mass[ len(mass)-1 ]
                    
                    rv = ReadInls('inlist_1.0','Reimers_scaling')
                    bv = ReadInls('inlist_1.0','Blocker_scaling')
                    
                    if fl < 0. :
                            masses.append(fm)
                            reims.append(rv)
                            block.append(bv)
                            print('Final mass is: '+str(fm))
                            print('Reimers: '+str(rv))
                            print('Blocker: '+str(bv))
                    else:
                            print('Luminosity too low here, failed point!')
                            failm.append(fm)
                            failb.append(bv)
                            failr.append(rv)
                
                    os.chdir(i)

        else:
            print(file+' is not a folder, ignoring it!')
         
#os.chdir('/data/mhoffman/NUQ/1M_grid_highres/')
os.chdir('/data/mhoffman/NUQ/CD_grid_highres_2')

if not os.path.exists('data'):
    os.mkdir('data')
#os.chdir('/data/mhoffman/NUQ/1M_grid_highres/data')
os.chdir('/data/mhoffman/NUQ/CD_grid_highres_2/data')
np.savetxt('StarM.out',masses,delimiter=',')
np.savetxt('Reims.out',reims,delimiter=',')
np.savetxt('Block.out',block,delimiter=',')
np.savetxt('FailM.out',failm,delimiter=',')
np.savetxt('FailB.out',failb,delimiter=',')
np.savetxt('FailR.out',failr,delimiter=',')

import csv

gooddata = [block,reims,masses]
baddata = [failb,failr,failm]

with open('gooddata.csv','a') as f:
        w = csv.writer(f)
        w.writerows(zip(*gooddata))

with open('baddata.csv','a') as f:
        w = csv.writer(f)
        w.writerows(zip(*baddata))
