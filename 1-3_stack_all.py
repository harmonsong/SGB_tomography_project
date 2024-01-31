# %%
import numpy as np
from ccfj import GetStationPairs
import os
import time
from geopy.distance import great_circle
import pandas as pd
import h5py
import yaml

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)

#%%
name_CC = 'CC/CC_40_prewhiten-onebit/'
name_stack = 'gather_all.h5'
dir_CC = os.path.join(dir_CC_workspace,name_CC)

# %%
# config file
filename = os.path.join(dir_CC+'info_CC.npy')
info_CC = np.load(filename, allow_pickle='TRUE').item()      # setting dictionary

# %%
d_start = info_CC['d_start']
d_end = info_CC['d_end']
d_start = 128
d_end = 158
print('d_start: ', d_start)
print('d_end: ', d_end)

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

# %% [markdown]
# ### Read stations info

# %%
stainfo = pd.read_excel(stalistname)
nsta = len(stainfo.iloc[:,0])
StationPairs = GetStationPairs(nsta)
nPairs = int(len(StationPairs)/2)
stalist = stainfo['Station'].tolist()
lat = stainfo['latitude'].tolist() 
lon = stainfo['longitude'].tolist()

# %% [markdown]
# ### calculate iterstational distance

# %%
r = np.zeros(nPairs)
for i in range(len(r)):
    r[i] = great_circle((lat[StationPairs[i*2]],lon[StationPairs[i*2]]),(lat[StationPairs[i*2+1]],lon[StationPairs[i*2+1]])).km

# %% [markdown]
# ### Stack

# %%
ncfs = np.zeros([nPairs,nf],dtype=np.complex64)
count = np.zeros(nPairs)

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

            ncfs[id1,:] = ncfs[id1,:]+data["ncfs"]
            count[id1] = count[id1]+1
            print(year+'-'+day+'   '+str(time.time()-start0)+' s')

# %%
ncfs1 = np.zeros(np.shape(ncfs),dtype=np.complex64)
for i in range(nPairs):
    if count[i]>0:
        ncfs1[i,:] = ncfs[i,:]/count[i]

# %% [markdown]
# ### Save

# %%
if os.path.exists(dir_CC+name_stack):
    os.remove(dir_CC+name_stack)
h5file = h5py.File(dir_CC+name_stack,'w')
h5file.create_dataset('ncfs',data=ncfs1)
h5file.create_dataset('count',data=count)
h5file.create_dataset('r',data=r)
h5file.create_dataset('f',data=f)
h5file.create_dataset('StationPairs',data=StationPairs)
h5file.close()
#np.savez(dir_CC+"summed.npz",ncfs= ncfs1,r = r)

# %%

filename = os.path.join(dir_CC+'info_CC.npy')
info_CC = np.load(filename, allow_pickle='TRUE').item()      # setting dictionary

# %%
info_CC


