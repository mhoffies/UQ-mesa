#!/usr/bin/python
import numpy as np
import random
import os, shutil, re

template = '1M_pre_ms_to_wd_template'
main_list = 'inlist_template'

def ChangeValue(inlist,newinl,param,newval):
    with open(inlist) as f:
        q = open(newinl,'a')
        for line in f:
            l = line
            if param in l:
                q.write('      '+param+' = '+str(newval)+'\n')
            if not param in l:
                q.write(l)
        
youarehere = os.getcwd()
gothere = 'test_suite'
pts = 200

print('Now creating top directory... '+gothere)
print('Number of points in paramters space... '+str(pts))
for i in range(pts):
    name = 'c'+str(i)
    rx = str(np.random.uniform(0.3,0.9))
    bx = str(np.random.uniform(0.0,0.1))
    shutil.copytree(template,gothere+'/'+name)
    shutil.copy(main_list,youarehere+'/'+gothere+'/'+name+'/.')
    os.chdir(youarehere+'/'+gothere+'/'+name)
    ChangeValue(main_list,'CHANGE_R','Reimers_scaling_factor',rx)
    ChangeValue('CHANGE_R','inlist_1.0','Blocker_scaling_factor',bx)
    os.chdir(youarehere+'/'+gothere)
    print('-----------------------------------------\n')
    print('Reimers changing to ...'+str(rx))
    print('Blocker changing to ...'+str(bx))
    print('Folder name: '+name+'\n')
    print('-----------------------------------------')
    os.chdir(youarehere)
