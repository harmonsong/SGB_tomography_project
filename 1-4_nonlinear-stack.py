# %%
import numpy as np
import obspy
from obspy import UTCDateTime
from ccfj import CC
from ccfj import GetStationPairs
from concurrent.futures import ThreadPoolExecutor
import os
import time
from geopy.distance import great_circle
from scipy.signal import hilbert
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import h5py
import yaml

# %%
import sys
sys.path.append(r'../')
from toollib_nonlinearstack import nonlinear_stacklib

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)

#%%
name_CC = 'CC/CC_40_onebit/'
name_stack = "gather_FPS.h5"
dir_CC = os.path.join(dir_CC_workspace,name_CC)

# %%
# config file
filename = os.path.join(dir_CC+'info_CC.npy')
info_CC = np.load(filename, allow_pickle='TRUE').item()      # setting dictionary

# %%
#d_start = info_CC['d_start']
#d_end = info_CC['d_end']
d_start = 128
d_end = 150
d_len = info_CC['d_len']
#d_end = d_start + d_len + 1
y_start = info_CC['y_start']
y_end = info_CC['y_end']
stalistname = info_CC['stalistname']
nf = info_CC['nf']
Fs = info_CC['Fs']
fftlen = info_CC['fftlen']
fstride = info_CC['fstride']
f = info_CC['f']

# %%
stainfo = pd.read_excel(stalistname)
nsta = len(stainfo.iloc[:,0])
StationPairs = GetStationPairs(nsta)
nPairs = int(len(StationPairs)/2)
stalist = stainfo['Station'].tolist()
lat = stainfo['latitude'].tolist() 
lon = stainfo['longitude'].tolist()

# %%
r = np.zeros(nPairs)
for i in range(len(r)):
    r[i] = great_circle((lat[StationPairs[i*2]],lon[StationPairs[i*2]]),(lat[StationPairs[i*2+1]],lon[StationPairs[i*2+1]])).km

# %% [markdown]
# ### Readccfs

# %%
# read ccfs
def Pairs(sta):
    p = []
    nsta = len(sta)
    for ii in range(nsta):
        for jj in range(ii+1,nsta):
            p.append([sta[ii],sta[jj]])
    return p

def cal_indx(pair,nsta):
    pair0 = min(pair)
    pair1 = max(pair)
    indx = int(pair0*(2*nsta-pair0-1)/2+pair1-pair0-1)
    return indx

# %%
def readccf(dir_CC):

    filename = os.path.join(dir_CC+'info_CC.npy')
    info_CC = np.load(filename, allow_pickle='TRUE').item()      # setting dictionary
    d_len = info_CC['d_len']
    d_start = info_CC['d_start']
    d_end = info_CC['d_end']  
    y_start = info_CC['y_start']
    y_end = info_CC['y_end']
    nf = info_CC['nf']
    stalistname = info_CC['stalistname'] + '.xlsx'
    stainfo = pd.read_excel(stalistname)
    nsta = len(stainfo.iloc[:,0])
    StationPairs = GetStationPairs(nsta)
    nPairs = int(len(StationPairs)/2)


    id1s = []
    ncfs = np.zeros([d_len,nPairs,nf],dtype = np.complex64)  
    for y in range(y_start,y_end):
        for d in range(d_start,d_end):
            start0 = time.time()
            year = str(y)
            day = "%03d"%d
            outname = os.path.join(dir_CC,year+'-'+day+'.npz')
            if os.path.exists(outname):
                data = []
                data = np.load(outname)
                #nsta0 = len(data["stalist"])
                indx = [stalist.index(i) for i in data["stalist"] ]
                pairs = Pairs(indx)
                id1 = [cal_indx(pair,nsta) for pair in pairs]


                ncfs[d-d_start,id1,:] = data["ncfs"]
                id1s.append(id1)
                print(year+'-'+day+'   '+str(time.time()-start0)+' s')
            
            idd = []
            for j in range(nPairs):
                idd1 = []
                for i in range(d_len):
                    idd1.append(i)
                idd.append(idd1)
    return ncfs, id1s, idd

# %%
ncfs, id1s, idd = readccf(dir_CC)

# %% [markdown]
# ### Stack

# %%
def freq(nr1,nf,ncfs,idd,amp):
    std = np.zeros([nr1,nf],dtype = complex)
    mean = np.zeros([nr1,nf],dtype = complex)

    for i in range(nr1):
        if len(idd[i]) > 0:
            mean[i,:] = np.mean(ncfs[idd[i],i,:],axis = 0)
            std[i,:] = np.std(ncfs[idd[i],i,:],axis = 0)
    s1 = np.zeros_like(mean)
    for j in range(nr1):
        for k in range(nf):
            s1[j,k] = np.mean(ncfs[np.abs(ncfs[:,j,k]-mean[j,k]) <= amp*std[j,k],j,k])
    return s1

# %%
amp = 1
s0 = freq(nPairs,nf,ncfs,idd,amp)
s1 = np.nan_to_num(s0)

# %%

h5file = h5py.File(dir_CC+name_stack,'w')
h5file.create_dataset('ncfs',data=s1)
h5file.create_dataset('r',data=r)
h5file.create_dataset('f',data=f)
h5file.create_dataset('StationPairs',data=StationPairs)
h5file.close()


