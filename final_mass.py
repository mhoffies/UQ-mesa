#!/usr/bin/python

# Here's what it does: using the NuGridPy module, this runs through a library of directories
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
DMC = []

youarehere = os.getcwd()

main_dir = input('Please provide topmost directory (i.e. we want main_dir from main_dir/name/c1) with \'\' :')
os.chdir(main_dir)
print('Currently in directory '+os.getcwd())
for file in os.listdir(main_dir):
    matchf = re.match('\Amdc\Z',file)
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
                        if 'mesh_delta_coeff' in line:
                            c = line
                            dmc1 = re.sub('      mesh\_delta\_coeff \= ','',c)
                            dmc2 = re.sub('d','e',dmc1)
                            DelMC = float(dmc2)
                            print(DelMC)
                            DMC.append(DelMC)
                            os.chdir(file_cab)
                else:
                    print('Final luminosity was not low enough, this point has been excluded!')
                    os.chdir(file_cab)
    os.chdir(main_dir)

R = np.asarray(Reims)
B = np.asarray(Block)
M = np.asarray(StarM)
D = np.asarray(DMC)

os.chdir(main_dir+'/mdc/')
if not os.path.exists('tmps'):
    os.mkdir('tmps')
os.chdir(main_dir+'/mdc/tmps')
print(os.getcwd())
print(StarM)
print(Block)
np.savetxt('StarM.out',StarM,delimiter=',')
np.savetxt('Reims.out',Reims,delimiter=',')
np.savetxt('Block.out',Block,delimiter=',')
np.savetxt('DMC.out',DMC,delimiter=',')

#R = np.asarray(Reims)
#B = np.asarray(Block)
#M = np.asarray(StarM)
#D = np.asarray(DMC)

R_Order = np.argsort(R)
RS = R[R_Order]
BS = B[R_Order]
MS = M[R_Order]
DS = D[R_Order]

Q = np.where(np.diff(RS[:]))[0]+1
FirstM = MS[0:Q[0]]
FirstB = BS[0:Q[0]]
FirstD = DS[0:Q[0]]
FirstR = RS[0:Q[0]]

Mdiv = []
Bdiv = []
Ddiv = []
Rdiv = []

Mdiv.append(FirstM)
Bdiv.append(FirstB)
Ddiv.append(FirstD)
Rdiv.append(FirstR)

for i in range(len(Q)-1):
    a = Q[i]+1
    b = Q[i+1]
    NewM = MS[a:b]
    NewB = BS[a:b]
    NewD = DS[a:b]
    NewR = RS[a:b]
    Mdiv.append(NewM)
    Bdiv.append(NewB)
    Ddiv.append(NewD)
    Rdiv.append(NewR)

LastM = MS[Q[len(Q)-1]:len(RS)-1]
LastB = BS[Q[len(Q)-1]:len(RS)-1]
LastD = DS[Q[len(Q)-1]:len(RS)-1]
LastR = RS[Q[len(Q)-1]:len(RS)-1]

Mdiv.append(LastM)
Bdiv.append(LastB)
Ddiv.append(LastD)
Rdiv.append(LastR)

SortedM = []
SortedB = []
SortedD = []
SortedR = []

for i in range(len(Mdiv)):
    B_order = np.argsort(Bdiv[i])
    UnSorM = Mdiv[i]
    UnSorB = Bdiv[i]
    UnSorD = Ddiv[i]
    UnSorR = Rdiv[i]
    SorM = UnSorM[B_order]
    SorB = UnSorB[B_order]
    SorD = UnSorD[B_order]
    SorR = UnSorR[B_order]
    SortedM.append(SorM)
    SortedB.append(SorB)
    SortedD.append(SorD)
    SortedR.append(SorR)

now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

os.chdir(youarehere)
