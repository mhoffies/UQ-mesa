
#!/usr/bin/python

# This code is for checking convergence. That's all. 

import os, shutil, re
import mesa as ms
import numpy as np
import matplotlib.pyplot as plt
import datetime
import collections

wheretogo = raw_input('Name of folder path: ')
os.chdir(wheretogo+'/tmps/')
here = os.getcwd()

# Read .out files
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

xR = keys['Reims.out']
xB = keys['Block.out']
xM = keys['StarM.out']
xD = keys['DMC.out']

Reims = data[xR]
Block = data[xB]
StarM = data[xM]
DMC = data[xD]

R = np.asarray(Reims)
B = np.asarray(Block)
M = np.asarray(StarM)
D = np.asarray(DMC)

# Breaks by Reimers
R_Order = np.argsort(R)
RS = R[R_Order]
BS = B[R_Order]
MS = M[R_Order]
DS = D[R_Order]
Q = np.where(np.diff(RS[:]))[0]+1
Mdiv = []
Bdiv = []
Ddiv = []
Rdiv = []

Mdiv.append(MS[0:Q[0]])
Bdiv.append(BS[0:Q[0]])
Ddiv.append(DS[0:Q[0]])
Rdiv.append(RS[0:Q[0]])

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

Mdiv.append(MS[Q[len(Q)-1]:len(RS)-1])
Bdiv.append(BS[Q[len(Q)-1]:len(RS)-1])
Ddiv.append(DS[Q[len(Q)-1]:len(RS)-1])
Rdiv.append(RS[Q[len(Q)-1]:len(RS)-1])
print(len(Ddiv))

SortedM = []
SortedB = []
SortedD = []
SortedR = []

# Breaks by Blocker 
for i in range(len(Mdiv)):
    B_order = np.argsort(Bdiv[i])
    SorM = Mdiv[i][B_order]
    SorB = Bdiv[i][B_order]
    SorD = Ddiv[i][B_order]
    SorR = Rdiv[i][B_order]
    QQ = np.where(np.diff(SorB[:]))[0]+1
    SortedM.append(SorM[0:QQ[0]])
    SortedB.append(SorB[0:QQ[0]])
    SortedD.append(SorD[0:QQ[0]])
    SortedR.append(SorR[0:QQ[0]])
    for i in range(len(QQ)-1):
        a = QQ[i]
        b = QQ[i+1]
        NewM = SorM[a:b]
        NewB = SorB[a:b]
        NewD = SorD[a:b]
        NewR = SorR[a:b]
        SortedM.append(NewM)
        SortedB.append(NewB)
        SortedD.append(NewD)
        SortedR.append(NewR)
    SortedM.append(SorM[QQ[len(QQ)-1]:len(SorM)-1])
    SortedB.append(SorB[QQ[len(QQ)-1]:len(SorM)-1])
    SortedD.append(SorD[QQ[len(QQ)-1]:len(SorM)-1])
    SortedR.append(SorR[QQ[len(QQ)-1]:len(SorM)-1])

# Orders by mdc
FinalM = []
FinalR = [] 
FinalB = []
FinalD = []
for i in range(len(SortedR)):
    D_order = np.argsort(SortedD[i])
    NewM = SortedM[i][D_order]
    NewB = SortedB[i][D_order]
    NewR = SortedR[i][D_order]
    NewD = SortedD[i][D_order]
    FinalM.append(NewM)
    FinalB.append(NewB)
    FinalR.append(NewR)
    FinalD.append(NewD)

Diff = []
Res = []
RGB = []
AGB = []
#Calculate difference in final masses
for i in range(len(FinalM)):
    aDiff = []
    aRes = []
    aRGB = []
    aAGB = []
    for j in range(len(FinalM[i])-1):
        diff = FinalM[i][j] - 0.5
        aDiff.append(FinalM[i][j])
        aRes.append(FinalD[i][j])
        aRGB.append(FinalR[i][j])
        aAGB.append(FinalB[i][j])    
    Diff.append(aDiff)
    Res.append(aRes)
    RGB.append(aRGB)
    AGB.append(aAGB)

# cmap = plt.get_cmap('jet')
# #Colors = [cmap(i) for i in
          
# for i in range(len(Res)):
#     if len(Diff[i])>0:
#         A = Diff[i]
#         print Diff[i]
#         B = Res[i]
#         C = RGB[i][0]+AGB[i][0]
#         label = 'R: '+str(RGB[i][0])+' B: '+str(AGB[i][0])
#         plt.plot(B,A,label=label) 
#     else:
#         print('whoops!')
# plt.show()    
    
now = datetime.datetime.now()
month = str(now.month)
day = str(now.day)
hour = str(now.hour)
minute = str(now.minute)

#os.chdir(youarehere)
