#!/usr/bin/python

# This code runs through completed mesa subdirectors (i.e. c1,c2,...) and reads the 
# out file to determine a) if the code completed b) whether it terminated for one reason or another
# At the moment, it looks for log_L_lower_limit or max_model_number, as these are the two termination
# cases that I'm interested in but you could change these to anything

import os, shutil, re
import numpy as np

term_lum = []
term_mod = []
term_tim = []
term_other = []
missing = []

# Checks the logfile to see if it finished
def FinCheck(dfile,term):
    found = False
    datalog = open(dfile,'r')
    for line in datalog:
        if term in line:
            found = True
            break

    return found

# Name of top directory and current directory
path = raw_input('Please provide path of top folder using: ')

os.chdir(path)
for file in os.listdir(path):

    # Check if it's a directory, if it's a file, we'll ignore it
    if os.path.isdir(path+'/'+file) == True:

        
        os.chdir(path+'/'+file)
        currd = os.getcwd()
        
        for file in os.listdir(currd):

            # Using Regex, try to match the log file 
            logtempl = re.match('\Arun\_out\_([0-9]*)\.log\Z',file)
            if logtempl:

                folder = os.getcwd()
                cn=logtempl.groups()[0]
                logname = 'run_out_'+cn+'.log'
                
                if FinCheck(logname,'termination code'):
                    #print('TRUE')
                    
                    with open(file,'r') as f:
                        for line in f:
                            if 'termination code:' in line:
                                if 'termination code: log_L_lower_limit' in line:
                                    name = os.path.split(folder)[-1]
                                    nme = str(name)
                                    num = nme.replace('c','')
                                    numb = int(num)
                                    #listy.append(numb)
                                    os.chdir(path)
                                    term_lum.append(numb)

                                if 'termination code: max_model_number' in line:
                                    name = os.path.split(folder)[-1]
                                    nme = str(name)
                                    num = nme.replace('c','')
                                    numb = int(num)
                                    #listy.append(numb)
                                    os.chdir(path)
                                    term_mod.append(numb)

                                if 'termination code: min_timestep_limit' in line:
                                    name = os.path.split(folder)[-1]
                                    nme = str(name)
                                    num = nme.replace('c','')
                                    numb = int(num)
                                    #listy.append(numb)
                                    os.chdir(path)
                                    term_tim.append(numb)

                                if not 'termination code: log_L_lower_limit' in line:
                                    if not 'termination code: max_model_number' in line:
                                        if not 'termination code: min_timestep_limit' in line:
                                            name = os.path.split(folder)[-1]
                                            nme = str(name)
                                            num = nme.replace('c','')
                                            numb = int(num)
                                            os.chdir(path)
                                            term_other.append(numb)
                                            os.chdir(path)
                            else:
                                os.chdir(path)
                else:
                    print('FALSE')
                    name = os.path.split(folder)[-1]
                    nme = str(name)
                    num = nme.replace('c','')
                    numb = int(num)
                    missing.append(numb)
                    os.chdir(path)
    else:
        print(file+' is not a directory, did not look for a .log file')

os.chdir(path)

# my_term_l holds those that terminated "correctly"
# because they reached a lower luminosity limit
my_term_l = np.asarray(term_lum)
reached_l = np.sort(my_term_l)
rl = len(reached_l)

# my_term_m holds those that terminated because of model numbers
# THIS IS A SILLY REASON, unless you had model number really high
# Then you've got some weird shit happening
my_term_m = np.asarray(term_mod)
max_models = np.sort(my_term_m)
mm = len(max_models)

# my_term_t holds those that reached minimum timestep
my_term_t = np.asarray(term_tim)
min_timst = np.sort(my_term_t)
mt = len(min_timst)

# those that terminated for other reasons 
my_term_o = np.asarray(term_other)
other_terms = np.sort(my_term_o)
mo = len(other_terms)

my_miss = np.asarray(missing)
unfinished = np.sort(my_miss)
lq = len(unfinished)

# Make lists that will translate well into text files

missingstr = ', '.join(map(str,unfinished))
lowlumstr = ', '.join(map(str,reached_l))
maxmodstr = ', '.join(map(str,max_models))
mintimstr = ', '.join(map(str,min_timst))
otherstr = ', '.join(map(str,other_terms))

if os.path.exists(path+'/report.out'):
    os.remove(path+'/report.out')
    print('Had to delete an old report!')

# Save everything to a report

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
