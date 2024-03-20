#! /usr/bin/env python3
import numpy as np
import os, getopt, sys
#------------------------------------------------------------
#Make folders used to contain different cell texts 
#recording the lagtime and cross-correlation coefficients 
#of ray paths crossing this cell
Forwardir = '/home/cjq/ModelValid/Forward/forLove/'

dir_3_7s = Forwardir + '/dir_3_7s'
dir_5_10s = Forwardir + '/dir_5_10s'
dir_8_15s = Forwardir +  '/dir_8_15s'
dir_12_20s = Forwardir  + '/dir_12_20s'

if not os.path.exists(dir_3_7s):
    os.makedirs(dir_3_7s)
if not os.path.exists(dir_5_10s):
    os.makedirs(dir_5_10s)
if not os.path.exists(dir_8_15s):
    os.makedirs(dir_8_15s)
if not os.path.exists(dir_12_20s):
    os.makedirs(dir_12_20s)
#=============================================================
for i in range(1,61):
#for i in range(1,2):
    try:
        cell0= np.loadtxt(Forwardir + 'src' + str(i) + '/' + '/cell_3_7s.txt')
    except:
        continue
    for j in range(len(cell0)):
        with open(dir_3_7s + '/' + str(int(cell0[j,0])) + '_' + str(int(cell0[j,1])) + '.txt', 'a+') as f:
            f.write(str(cell0[j,2]) + ' ' + str(cell0[j,3]) + '\n')

for i in range(1,61):
#for i in range(1,2):
    try:
        cell1= np.loadtxt(Forwardir + 'src' + str(i) + '/' + '/cell_5_10s.txt')
    except:
        continue
    for j in range(len(cell1)):
        with open(dir_5_10s + '/' + str(int(cell1[j,0])) + '_' + str(int(cell1[j,1])) + '.txt', 'a+') as f:
            f.write(str(cell1[j,2]) + ' ' + str(cell1[j,3]) + '\n')

for i in range(1,61):
#for i in range(1,2):
    try:
        cell2= np.loadtxt(Forwardir + 'src' + str(i) + '/' + '/cell_8_15s.txt')
    except:
        continue
    for j in range(len(cell2)):
        with open(dir_8_15s + '/' + str(int(cell2[j,0])) + '_' + str(int(cell2[j,1])) + '.txt', 'a+') as f:
            f.write(str(cell2[j,2]) + ' ' + str(cell2[j,3]) + '\n')

for i in range(1,61):
#for i in range(1,2):
    try:
        cell3= np.loadtxt(Forwardir + 'src' + str(i) + '/' + '/cell_12_20s.txt')
    except:
        continue
    for j in range(len(cell3)):
        with open(dir_12_20s + '/' + str(int(cell3[j,0])) + '_' + str(int(cell3[j,1])) + '.txt', 'a+') as f:
            f.write(str(cell3[j,2]) + ' ' + str(cell3[j,3]) + '\n')

    print("src"+ str(i) + " Done!")
