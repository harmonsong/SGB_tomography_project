import os,pandas,h5py,itertools
import numpy as np
from geopy.distance import great_circle
from multiprocessing import Manager
import matplotlib.pyplot as plt
from multiprocessing import Pool as ThreadPool
import yaml
import pandas as pd

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_SAC_workspace = dir_config['dir_SAC_workspace']
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_SAC_workspace: ', dir_SAC_workspace)
print('dir_CC_workspace: ', dir_CC_workspace)

#%% 
name_CC = 'CC/CC_100_prewhiten/'  
dir_CC = os.path.join(dir_CC_workspace,name_CC)
filename = dir_CC+'info_CC.npy'
info_CC = np.load(filename, allow_pickle='TRUE').item()

#%% 
outpath = dir_CC + 'CFs_modelvalidate/'
out_file = outpath + 'gather_all_modelvalidate.h5'
if not os.path.exists(outpath):
    os.mkdir(outpath)

#%%
stalist = []
lon = []
lat = []
stalistname_all = info_CC['stalistname']
stainfo = pd.read_excel(stalistname_all)
stalist = stainfo['Station'].tolist() 
lat =  stainfo['latitude'].tolist() 
lon =  stainfo['longitude'].tolist() 

#%% 
original_stack = h5py.File(dir_CC + 'gather_all.h5','r')
R = original_stack['r'][:]
Stalist = stalist
Ncfs = original_stack['ncfs'][:]
StationPairs = original_stack['StationPairs'][:]

#%% 
CC_name = []
Dis_great=[]
CC_data = []
sta0_lon_lat0=[]
sta1_lon_lat0=[]
jj = 0
for i in range(len(R)):
    sta0 = Stalist[StationPairs[i*2]]
    sta1 = Stalist[StationPairs[i*2+1]]
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

freq = pandas.DataFrame({'f': original_stack["f"]})
if os.path.exists(out_file):
    os.remove(out_file)
f = h5py.File(out_file,'a')
f.create_dataset('all', data=ncfs0)
f.close()
cc_prof.to_hdf(out_file, key='cc_prof',mode = 'r+')
freq.to_hdf(out_file, key='freq',mode = 'r+')