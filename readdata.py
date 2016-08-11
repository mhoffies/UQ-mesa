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
                        val = float( a.split()[-1] )
                        return val

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

paths = ''

for i in paths:
    os.chdir(i)

    # Get the final masses from the folders
    #masses = []
    for file in os.listdir(i):
    
        # Check if it's a directory, if it's a file, we'll ignore it
        if os.path.isdir(i+'/'+file) == True:

            os.chdir(i+'/'+file)
            currd = os.getcwd()

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
            else:
                print('In '+currd+' luminosity too low, excluding point!')

                os.chdir(i)

        else:
            print(file+' is not a folder, ignoring it!')




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
