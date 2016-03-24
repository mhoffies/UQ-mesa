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

youarehere = os.getcwd()
# Move into directories 
# main_dir = os.getcwd()
# Input top most directory
main_dir = input('Please provide topmost directory (i.e. we want main_dir from main_dir/name/c1) with \'\' :')
os.chdir(main_dir)
print('Currently in directory '+os.getcwd())
for file in os.listdir(main_dir):
    matchf = re.match('\Abloc\-rand\_reim\-[0-9]\.[0-9]\Z',file)
    #matchf = re.match('\Areim\_breadth\_test\Z',file)
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
                    os.chdir(file_cab)
    os.chdir(main_dir)
    
if not os.path.exists('tmps'):
    os.mkdir('tmps')
os.chdir(main_dir+'/tmps')
np.savetxt('StarM.out',StarM,delimiter=',')
np.savetxt('Reims.out',Reims,delimiter=',')
np.savetxt('Block.out',Block,delimiter=',')

# Now let's graph this! We have here Mass v. Blocker, with Reimers colormap but you could totally change it up! 
# Ask User Questions about what to Plot instead!

# color = input('What colormap would you like (not rainbow!)')

# Xax = input('What would you like your x-axis to be (StarM,Lumi,Block,Reims,Temp)?')
# xlab = raw_input('X-axis name:')
# Yax = input('What would you like your y-axis to be (StarM,Lumi,Block,Reims,Temp)?')
# ylab = raw_input('Y-axis name:')
# Colr = input('What should the color mapping be to (StarM,Lumi,Block,Reims,Temp)?')
# clab = raw_input('Colormap label name:')

# xmx = round(max(Xax),2)+0.01
# xmn = round(min(Xax),2)-0.01
# ymx = round(max(Yax),2)+0.01
# ymn = round(min(Yax),2)-0.01
# cmx = round(max(Colr),2)
# cmn = round(min(Colr),2)

#cols=plt.cm.gnuplot_r
#fig, ax = plt.subplots()
#plot = ax.scatter(Xax,Yax,s=50,c=Colr,cmap=cols,vmin=cmn,vmax=cmx,edgecolor='none')
#cbar = fig.colorbar(plot)
# Make numpy arrays

R = np.asarray(Reims)
B = np.asarray(Block)
M = np.asarray(StarM)

R_Order = np.argsort(R)
RS = R[R_Order]
BS = B[R_Order]
MS = M[R_Order]

Q = np.where(np.diff(RS[:]))[0]+1
FirstM = MS[0:Q[0]]
FirstB = BS[0:Q[0]]

Mdiv = []
Bdiv = []
Mdiv.append(FirstM)
Bdiv.append(FirstB)

for i in range(len(Q)-1):
    a = Q[i]+1
    b = Q[i+1]
    NewM = MS[a:b]
    NewB = BS[a:b]
    Mdiv.append(NewM)
    Bdiv.append(NewB)

LastM = MS[Q[len(Q)-1]:len(RS)-1]
LastB = BS[Q[len(Q)-1]:len(RS)-1]

Mdiv.append(LastM)
Bdiv.append(LastB)

SortedM = []
SortedB = []

for i in range(len(Mdiv)):
    B_order = np.argsort(Bdiv[i])
    UnSorM = Mdiv[i]
    UnSorB = Bdiv[i]
    SorM = UnSorM[B_order]
    SorB = UnSorB[B_order]
    SortedM.append(SorM)
    SortedB.append(SorB)

plt.plot(SortedB[0],SortedM[0],color='violet',label='0.3')
plt.plot(SortedB[1],SortedM[1],color='r',label='0.4')
plt.plot(SortedB[2],SortedM[2],color='orange',label='0.5')
plt.plot(SortedB[3],SortedM[3],color='y',label='0.6')
plt.plot(SortedB[4],SortedM[4],color='g',label='0.7')
plt.plot(SortedB[5],SortedM[5],color='b',label='0.8')
plt.plot(SortedB[6],SortedM[6],color='indigo',label='0.9')
plt.plot(SortedB[7],SortedM[7],color='gray',label='1.0')
plt.plot(SortedB[8],SortedM[8],color='k',label='1.1')


# ax.set_title(ylab+' v. '+xlab)
# ax.title.set_fontsize(20)
# ax.set_xlabel(xlab)
# ax.set_ylabel(ylab)
# ax.xaxis.label.set_fontsize(18)
# ax.set_xlim([xmn,xmx])
# ax.set_ylim([ymn,ymx])
# ax.yaxis.label.set_fontsize(18)
# ax.grid()

# cbar.set_label(clab)

now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

plt.savefig(month+day+'.'+hour+minute+'.png')

plt.show()

os.chdir(youarehere)
