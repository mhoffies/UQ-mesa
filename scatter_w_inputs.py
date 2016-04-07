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
Temp = []

failM = []
failR = []
failB = []
failL = []

youarehere = os.getcwd()
# Move into directories 
# main_dir = os.getcwd()
# Input top most directory
main_dir = input('Please provide topmost directory (i.e. we want main_dir from main_dir/name/c1) with \'\' :')
os.chdir(main_dir)
print('Currently in directory '+os.getcwd())
for file in os.listdir(main_dir):
    matchf = re.match('\Afine\_grid\_([0-9]*)\-([0-9]*)\Z',file)
    #matchf = re.match('\Abloc\-rand\_reim\-[0-9]\.[0-9]\Z',file)
    #matchf = re.match('\Areim\_breadth\_test\Z',file)
    #matchf = re.match('\Ares\_study\_reim\-[0-9]\.[0-9]\Z',file)
    #matchf = re.match('\Areim\-probe\_bloc\-0\.1\Z',file)
    if matchf:
        os.chdir(file)
        file_cab = os.getcwd()
        # print file_cab
        for file in os.listdir(file_cab):
            matchc = re.match('\Ac([0-9]*)\Z',file)
            if matchc:
                os.chdir(file)
                print('In directory:'+os.getcwd())
                # Get Stellar Mass
                s = ms.history_data()
                mass = s.get('star_mass')
                lm = len(mass)
                lmm = lm - 1
                fin_M = mass[lmm]
                lum = s.get('log_L')
                fin_L = lum[lmm]
                tem = s.get('log_Teff')
                fin_T = tem[lmm]
                #print(os.getcwd())
                print('Final mass:'+str(fin_M))
                print('Final lumi:'+str(fin_L))
                #print('Final temp:'+str(fin_T))
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
                    failL.append(fin_L)
                    failM.append(fin_M)
                    inl = open('inlist_1.0','r')
                    for line in inl:
                        if 'Reimers_scaling_' in line:
                            a = line
                            lastR = re.sub('      Reimers\_scaling\_factor \= ','',a)
                            Rval = re.sub('d','e',lastR)
                            R = float(Rval)
                            print(R)
                            failR.append(R)
                        if 'Blocker_scaling_' in line:
                            b = line
                            lastB = re.sub('      Blocker\_scaling\_factor \= ','',b)
                            Bval = re.sub('d','e',lastB)
                            B = float(Bval)
                            print(B)
                            failB.append(B)
                            os.chdir(file_cab)
    os.chdir(main_dir)
    
if not os.path.exists('tmps'):
    os.mkdir('tmps')
os.chdir(main_dir+'/tmps')
np.savetxt('StarM.out',StarM,delimiter=',')
np.savetxt('Reims.out',Reims,delimiter=',')
np.savetxt('Block.out',Block,delimiter=',')
np.savetxt('failM.out',failM,delimiter=',')
np.savetxt('failB.out',failB,delimiter=',')
np.savetxt('failR.out',failR,delimiter=',')
np.savetxt('failL.out',failL,delimiter=',')

# Now let's graph this! We have here Mass v. Blocker, with Reimers colormap but you could totally change it up! 
# Ask User Questions about what to Plot instead!

# color = input('What colormap would you like (not rainbow!)')

Xax = input('What would you like your x-axis to be (StarM,Lumi,Block,Reims,Temp)?')
xlab = raw_input('X-axis name:')
Yax = input('What would you like your y-axis to be (StarM,Lumi,Block,Reims,Temp)?')
ylab = raw_input('Y-axis name:')
Colr = input('What should the color mapping be to (StarM,Lumi,Block,Reims,Temp)?')
clab = raw_input('Colormap label name:')

xmx = round(max(Xax),2)+0.01
xmn = round(min(Xax),2)-0.01
ymx = round(max(Yax),2)+0.01
ymn = round(min(Yax),2)-0.01
cmx = round(max(Colr),2)
cmn = round(min(Colr),2)

cols=plt.cm.gnuplot_r
fig, ax = plt.subplots()
plot = ax.scatter(Xax,Yax,s=50,c=Colr,cmap=cols,vmin=cmn,vmax=cmx,edgecolor='none')
cbar = fig.colorbar(plot)

ax.set_title(ylab+' v. '+xlab)
ax.title.set_fontsize(20)
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax.xaxis.label.set_fontsize(18)
ax.set_xlim([xmn,xmx])
ax.set_ylim([ymn,ymx])
ax.yaxis.label.set_fontsize(18)
ax.grid()

cbar.set_label(clab)

now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

plt.savefig(month+day+'.'+hour+minute+'.png')

plt.show()

os.chdir(youarehere)
