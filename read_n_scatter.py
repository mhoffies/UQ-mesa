#!/usr/bin/python

# This code should finally be okay.
# Here's what it does: using the NuGridPy module (website here) this runs through a library of directories
# With the same name format (i.e. bloc-rand_reim-#.#), then runs through their subdirectories (i.e. c##)
# And extracts the final parameter (Mass, Luminosity) using the mesa module and extracts information from
# the inlists (blocker_scaling_factor, reimers_scaling_factor). It then plots Final Mass V. Blocker
# and colormaps is using Reimers. Saves this as date.time.png

import os, shutil, re
import mesa as ms
import numpy as np
import matplotlib.pyplot as plt
import datetime

# Make arrays to store Final Stellar Mass, Reimers, and Blocker
StarM = []
Reims = []
Block = []
Lumi = []

youarehere = os.getcwd()
# Move into directories 
# main_dir = os.getcwd()
# Input top most directory
main_dir = input('Please provide topmost directory (i.e. we want main_dir from main_dir/name/c1 with \'\'')
os.chdir(main_dir)
print('Currently in directory '+os.getcwd())
for file in os.listdir(main_dir):
    matchf = re.match('\Abloc\-rand\_reim\-[0-9]\.[0-9]\Z',file)
    if matchf:
        os.chdir(file)
        file_cab = os.getcwd()
        # print file_cab
        for file in os.listdir(file_cab):
            matchc = re.match('\Ac([0-9]*)\Z',file)
            if matchc:
                os.chdir(file)
                # Get Stellar Mass
                s = ms.history_data()
                mass = s.get('star_mass')
                lm = len(mass)
                lmm = lm - 1
                fin_M = mass[lmm]
                lum = s.get('log_L')
                fin_L = lum[lmm]
                print(os.getcwd())
                print(fin_M)
                print(fin_L)
                if fin_L < 0:
                    Lumi.append(fin_L)
                    StarM.append(fin_M)
                    # Get Reimers and Blocker
                    inl = open('inlist_1.0','r')
                    for line in inl:
                        if 'Reimers_scaling_' in line:
                            a = line
                            lastR = re.sub('      Reimers\_scaling\_factor \= ','',a)
                            Rval = re.sub('d','e',lastR)
                            R = float(Rval)
                            print(R)
                            Reims.append(R)
                        if 'Blocker_scaling_' in line:
                            b = line
                            lastB = re.sub('      Blocker\_scaling\_factor \= ','',b)
                            Bval = re.sub('d','e',lastB)
                            B = float(Bval)
                            print(B)
                            Block.append(B)
                            os.chdir(file_cab)
                else:
                    print('Final luminosity was not low enough, this point has been excluded!')
                    os.chdir(file_cab)
    os.chdir(main_dir)
    
if not os.path.exists('tmps'):
    os.mkdir('tmps')
os.chdir(main_dir+'/tmps')
np.savetxt('StarM.out',StarM,delimiter=',')
np.savetxt('Reims.out',Reims,delimiter=',')
np.savetxt('Block.out',Block,delimiter=',')

# Now let's graph this! We have here Mass v. Blocker, with Reimers colormap but you could totally change it up! 

rmax = round(max(Reims),1)
rmin = round(min(Reims),1)
bmax = round(max(Block),1)
bmin = round(min(Block),1)
mmax = round(max(StarM),1)
mmin = round(min(StarM),1)

cols=plt.cm.gnuplot_r

fig, ax = plt.subplots()
plot = ax.scatter(Block,StarM,s=50,c=Reims,cmap=cols,vmin=rmin,vmax=rmax,edgecolor='none')
cbar = fig.colorbar(plot)

ax.set_title('Final Star Mass v. Blocker Scaling Factor')
ax.set_xlabel('Blocker')
ax.set_ylabel('Final Mass')
ax.set_xlim([bmin,bmax])
ax.set_ylim([mmin,mmax])
ax.grid()

now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

plt.savefig(month+day+'.'+hour+minute+'.png')

plt.show()

os.chdir(youarehere)
