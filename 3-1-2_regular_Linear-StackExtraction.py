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
import folium
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import yaml

# %%
import sys
sys.path.append(r'../tools_F-J/')
from toollib_standard import maplib
from toollib_standard import filelib
from toollib_standard import stacklib

# %%
with open('a-project.yml', 'r', encoding='utf-8') as f:
    proj = yaml.load(f.read(), Loader=yaml.FullLoader)
name_project = proj['name']
#name_project = 'project/output_FJSJ_16-01/'               # Harmon server
name_project

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_project_workspace = dir_config['dir_project_workspace']
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)
print('dir_project_workspace: ', dir_project_workspace)
dir_project = dir_project_workspace + name_project
print('dir_project: ', dir_project)

#%% 
filename = dir_project+'Basic_info.yml'
with open(filename, 'r', encoding='utf-8') as f:
    info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
filename_bi = dir_project+'Basic_info.npy'
info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary

# %%
dir_stack = dir_project+info_basic['dir_stack']
dir_CC = dir_CC_workspace+info_basic['name_CC']
stalistname_all = info_basic['stalistname_all']

# %% [markdown]
# ### Define how many subworks to be done

# %%
key_subworks = ['01-01']
key_subworks = info_basic['key_subworks']
key_subworks

# %% [markdown]
# ### Read All data

# %%
stainfo = pd.read_excel(stalistname_all)
nsta_all = len(stainfo.iloc[:,0])
StationPairs_all = GetStationPairs(nsta_all)
nPairs_all = int(len(StationPairs_all)/2)
stalist_all = stainfo['Station'].tolist()
lat_all = stainfo['latitude'].tolist() 
lon_all = stainfo['longitude'].tolist()

# %%
ncffile = h5py.File(dir_CC + 'gather_all.h5','r')
ncfs = ncffile['ncfs'][:]
f = ncffile['f'][:]
count_all = ncffile['count'][:]
ncffile.close()

# %% [markdown]
# ### Collect Stack files

# %%
info_basic_bi['r_max'] = {}

# %%
for key_subwork in key_subworks:
    print("Collecting ",key_subwork,' ...')
    time0 = time.time()

    nf = info_basic['nf']
    #f = info_basic['f']
    dir_stack= dir_project + info_basic['dir_stack']
    
    nsta = info_basic['nstaS'][key_subwork]
    nPairs = int(nsta*(nsta-1)/2)
    stalistname = dir_project+info_basic['stalistname']

    stainfo = pd.read_excel(stalistname,key_subwork)
    nsta = len(stainfo.iloc[:,0])
    StationPairs = GetStationPairs(nsta)
    nPairs = int(len(StationPairs)/2)
    stalist = stainfo['Station'].tolist()
    lat = stainfo['latitude'].tolist() 
    lon = stainfo['longitude'].tolist()

    ncfs_sum_linear = np.zeros((nPairs,nf),dtype=np.complex64)
    r = np.zeros(nPairs)
    count= np.zeros(nPairs)
    StationPairs = GetStationPairs(nsta)
    names = []
    for i in range(nPairs):
        sta1 = StationPairs[2*i]
        sta2 = StationPairs[2*i+1]
        idx1 = np.min( [int(stalist_all.index(stalist[sta1])),int(stalist_all.index(stalist[sta2]))] )
        idx2 = np.max( [int(stalist_all.index(stalist[sta1])),int(stalist_all.index(stalist[sta2]))] )
        #idx1 = int(stalist_all.index(stalist[sta1]))
        #idx2 = int(stalist_all.index(stalist[sta2]))
        #idx1 = int(stainfo[stainfo[key_pd]==stalist[sta1]].index.values[0])
        #idx2 = int(stainfo[stainfo[key_pd]==stalist[sta2]].index.values[0])
        
        m = 0
        for j in range(nsta_all-idx1,nsta_all):
            m += j
        num = m +idx2 - idx1 -1
        
        ncfs_sum_linear[i,:] = np.nan_to_num(ncfs[num,:])
        count[i] = count_all[num]
        #r[i] = r0[num]
        #if count_all[num] > 0:
        #    ncfs_sum_linear[i,:] = ncfs[num,:]/count_all[num]
        #    count[i] = count_all[num]

        r[i] = great_circle((lat[sta1],lon[sta1]),(lat[sta2],lon[sta2])).km
        #names.append([stalist_all[StationPairs_all[2*num]],stalist_all[StationPairs_all[2*num+1]]])

    outname = key_subwork+'_gather_linear.h5'
    if os.path.exists(dir_stack+outname):
        os.remove(dir_stack+outname)
    ncffile = h5py.File(dir_stack+outname,'w')
    
    ncffile.create_dataset('ncfs',data=ncfs_sum_linear)
    ncffile.create_dataset('r',data=r)
    ncffile.create_dataset('count',data=count)
    ncffile.create_dataset('f',data=f)
    ncffile.create_dataset('StationPairs',data=StationPairs)
    #print(ncffile.keys())
    ncffile.close()
    #np.savez(dir_stack+key_subwork+"_summed-linear.npz",ncfs= ncfs_sum_linear,r = r,stalist=stalist,StationPairs=StationPairs)
    print("Done in ",time.time()-time0," s.")

    info_basic_bi['r_max'][key_subwork] = max(r)

# %%
np.save(filename_bi,info_basic_bi)