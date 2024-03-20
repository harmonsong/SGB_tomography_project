#!/usr/bin/env python3
import os,pandas,h5py,itertools
import numpy as np
from geopy.distance import great_circle
from multiprocessing import Manager
import matplotlib.pyplot as plt
from multiprocessing import Pool as ThreadPool

#npzdir = './test_outputCFs'
CCnpz = np.load('/data/wp/freq_CF.npz')
outpath = './CFs/' 
out_file = './CFs/WP_summed2.h5'
#out_file = 'test_H1_summed.h5'
stalistname = '/home/cjq/WP/WithinSta.txt'  #netname, staname, lon, lat

if not os.path.exists(outpath):
    os.makedirs(outpath)
#-------------------------------------------------------------------------------
stalist = []
lon = []
lat = []
with open(stalistname,'r') as f:
    while True:
        tmp = f.readline()
        if tmp:
            stalist.append(tmp.split()[0]+'-'+tmp.split()[1])
            lon.append(float(tmp.split()[2]))
            lat.append(float(tmp.split()[3]))
        else:
            break

R = CCnpz['r']
Stalist = CCnpz['stalist']
StationPairs = CCnpz['StationPairs']
Ncfs = CCnpz['ncfs']
CC_name = []
Dis_great=[]
CC_data = []
sta0_lon_lat0=[]
sta1_lon_lat0=[]
jj = 0
for i in range(len(R)):
    sta0 = Stalist[StationPairs[i,0]]
    sta1 = Stalist[StationPairs[i,1]]
    if (sta0 in stalist and sta1 in stalist):
        cc_name = sta0 + '_' + sta1
        CC_name.append(cc_name)
        id_sta0 = stalist.index(sta0)
        id_sta1 = stalist.index(sta1)
        sta0_lon_lat0.append([float(lon[id_sta0]),float(lat[id_sta0])])
        sta1_lon_lat0.append([float(lon[id_sta1]),float(lat[id_sta1])])
        Dis_great.append(R[i])
        CC_data.append(Ncfs[i])
        jj = jj+1
        print(str(jj)+' of '+str(len(R)))

                    

r0 = np.array(Dis_great)
ncfs0 = np.array(CC_data)
#-------------------------------------------------------------------------------#
cc_prof = pandas.DataFrame({'0cc_name': CC_name,
                            '1sta0': sta0_lon_lat0,
                            '2sta1': sta1_lon_lat0,
                            '3dis_great': r0,
                            })

freq = pandas.DataFrame({'f': CCnpz["f"]})
f = h5py.File(out_file,'a')
f.create_dataset('all', data=ncfs0)
f.close()
cc_prof.to_hdf(out_file, key='cc_prof',mode = 'r+')
freq.to_hdf(out_file, key='freq',mode = 'r+')
#---------------------------------------------------------------------------#
