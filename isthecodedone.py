#!/usr/bin/python

# This code runs through completed mesa subdirectors (i.e. c1,c2,...) and reads the 
# out file to determine a) if the code completed b) whether it terminated for one reason or another
# At the moment, it looks for log_L_lower_limit or max_model_number, as these are the two termination
# cases that I'm interested in but you could change these to anything

import os, shutil, re
import numpy as np

listy = []
term_lum = []
term_mod = []
term_tim = []
term_other = []

file_cab = raw_input('Please provide path of top folder using: ')
nf = input('Number of folders? ')
os.chdir(file_cab)
print('Now in directory'+os.getcwd())
for file in os.listdir(file_cab):
    matchc = re.match('\Ac([0-9]*)\Z',file)
    if matchc:
#        print file
        os.chdir(file)
        whichc = os.getcwd()
#        print(whichc)
        for file in os.listdir(whichc):
            logtempl = re.match('\Amesa\_pmswd_([0-9]*)\.log\Z',file)
            if logtempl:
                print(file)
                with open(file,'r') as f:
                    for line in f:
                        # Uncomment this section if you only care
                        # whether the run has finished or not
                        if 'termination code:' in line:
                            #print(line)
                            #name = os.path.split(whichc)[-1]
                            #nme = str(name)
                            #num = nme.replace('c','')
                            #numb = int(num)
                            #listy.append(numb)
                            #print(os.path.split(whichc)[-1]+' has completed its run')
                            #os.chdir(file_cab)
                            if 'termination code: log_L_lower_limit' in line:
                                name = os.path.split(whichc)[-1]
                                nme = str(name)
                                num = nme.replace('c','')
                                numb = int(num)
                                listy.append(numb)
                                os.chdir(file_cab)
                                term_lum.append(numb)
                            if 'termination code: max_model_number' in line:
                                name = os.path.split(whichc)[-1]
                                nme = str(name)
                                num = nme.replace('c','')
                                numb = int(num)
                                listy.append(numb)
                                os.chdir(file_cab)
                                term_mod.append(numb)
                            if 'termination code: min_timestep_limit' in line:
                                name = os.path.split(whichc)[-1]
                                nme = str(name)
                                num = nme.replace('c','')
                                numb = int(num)
                                listy.append(numb)
                                os.chdir(file_cab)
                                term_tim.append(numb)
                            if not 'termination code: log_L_lower_limit' in line:
                                if not 'termination code: max_model_number' in line:
                                    if not 'termination code: min_timestep_limit' in line:
                                        name = os.path.split(whichc)[-1]
                                        nme = str(name)
                                        num = nme.replace('c','')
                                        numb = int(num)
                                        listy.append(numb)
                                        os.chdir(file_cab)
                                        term_other.append(numb)
                        else:
                            os.chdir(file_cab)

mylisty = np.asarray(listy)
q = np.sort(mylisty)
lq = len(q)
missing = []
my_term_l = np.asarray(term_lum)
reached_l = np.sort(my_term_l)
rl = len(reached_l)
my_term_m = np.asarray(term_mod)
max_models = np.sort(my_term_m)
mm = len(max_models)
my_term_t = np.asarray(term_tim)
min_timst = np.sort(my_term_t)
mt = len(min_timst)
my_term_o = np.asarray(term_other)
other_terms = np.sort(my_term_o)
mo = len(other_terms)

for i in range(1,(int(nf)+1)):
    if not i in q:
        missing.append(i)
missingstr = ', '.join(map(str,missing))
lowlumstr = ', '.join(map(str,reached_l))
maxmodstr = ', '.join(map(str,max_models))
mintimstr = ', '.join(map(str,min_timst))
otherstr = ', '.join(map(str,other_terms))

with open('report.out','a') as f:
    f.write(str(lq)+' directories have completed their runs. Those missing:\n')
    f.write(missingstr+'\n')
    f.write(str(rl)+' directories reached 0.1 solar luminosities\n')
    f.write(lowlumstr+'\n')
    f.write(str(mm)+' directories terminated because maximum model number was achieved\n')
    f.write(maxmodstr+'\n')
    f.write(str(mt)+' directories had problems converging -- investigate furthur!\n')
    f.write(mintimstr+'\n')
    f.write(str(mo)+' directors terminated for some other reason - GO LOOK WHY!\n')
    f.write(otherstr+'\n')

f = open('report.out','r')
file_contents = f.read()
print(file_contents)
f.close()
