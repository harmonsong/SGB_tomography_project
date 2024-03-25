# %%
import matplotlib.pyplot as plt
import numpy as np
import ccfj
from ccfj import CC
from ccfj import GetStationPairs
from geopy.distance import great_circle
import folium
import h5py
import time
import os
from concurrent.futures import ThreadPoolExecutor
import yaml
import pandas as pd

# %%
import sys
sys.path.append(r'../tools_F-J/')
from toollib_standard import maplib
from toollib_standard import mathlib
from toollib_standard import filelib
from toollib_standard import plotlib

# %%
flag_save = 0; # 1: save stack, 0: not save stack

#%%
with open('a-project_repar.yml', 'r', encoding='utf-8') as f:
    proj = yaml.load(f.read(), Loader=yaml.FullLoader)
name_project = proj['name']
#name_project = 'project_repartrition/repartrition_01-03/'               # Harmon server
name_project

#%%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_project_workspace = dir_config['dir_project_workspace']
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)
print('dir_project_workspace: ', dir_project_workspace)
dir_project = os.path.join(dir_project_workspace, name_project)
print('dir_project: ', dir_project)

# %%
filename = dir_project+'Basic_info.yml'
with open(filename, 'r', encoding='utf-8') as f:
    info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
filename_bi = dir_project+'Basic_info.npy'
info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary

#%%
key_subworks = info_basic['key_subworks']
if 'key_subworks_repick' in info_basic.keys():
    key_subworks_repick = info_basic['key_subworks_repick']
else:
    key_subworks_repick = []
#print('key_subworks_repick: ', key_subworks_repick)
#key_subworks = key_subworks[261:300]
#key_subworks = ['1']

# %%
dir_CC = dir_CC_workspace+info_basic['rdir_CC']
dir_partition = dir_project + info_basic['rdir_partition']
stalistname_all = info_basic['stalistname_all']
rdir_ds = 'ds_'+info_basic['tag'] + '/'
info_basic['rdir_ds'] = rdir_ds
dir_ds = dir_project + rdir_ds
if os.path.exists(dir_ds) == False:
    os.makedirs(dir_ds)

#%% 
stainfo = pd.read_excel(stalistname_all)
nsta_all = len(stainfo.iloc[:,0])
StationPairs_all = GetStationPairs(nsta_all)
nPairs_all = int(len(StationPairs_all)/2)
stalist_all = stainfo['Station'].tolist()

#%% 
ncffile = h5py.File(dir_CC + 'gather_all.h5','r')
ncfs = ncffile['ncfs'][:]
f = ncffile['f'][:]
count_all = ncffile['count'][:]
ncffile.close()

#%%
def noise_fj(key_subwork,ncfs_linear,ncfs_remove,r,index):
    global dir_project
    global dir_ds
    global c
    global key_subworks
    global info_basic
    global info_basic_bi

    start0 = time.time()
    if os.path.exists(dir_ds+'ds_'+key_subwork+'.h5'):
        os.remove(dir_ds+'ds_'+key_subwork+'.h5')

    f0 = info_basic_bi['f']
    ds11 = ccfj.fj_noise(np.real(ncfs_linear),r,c,f0,fstride=1,itype=1,func=1)
    ds_linear = np.array(ds11).reshape(1,np.shape(ds11)[0],np.shape(ds11)[1])
    ds11 = ccfj.fj_noise(np.real(ncfs_remove),r,c,f0,fstride=1,itype=1,func=1)
    ds_remove = np.array(ds11).reshape(1,np.shape(ds11)[0],np.shape(ds11)[1])
    
    h5file = h5py.File(dir_ds+'ds_'+key_subwork+'.h5','w')
    h5file.create_dataset('ds_remove',data=ds_remove)
    h5file.create_dataset('ds_linear',data=ds_linear)
    h5file.create_dataset('f',data=f0)
    h5file.create_dataset('c',data=c)
    h5file.create_dataset('index_ncfs',data = index)
    h5file.create_dataset('r',data=r)
    h5file.close()
    print('Finish FJ calculation, time:', time.time()-start0, ' seconds')

#%%
def time_window_filter(t,ncfst0,r,v_min,t0,a):
    ncfst = ncfst0.copy()
    for i in range(len(ncfst)):
        tag = r[i]/v_min
        #print(t0,tag)
        t1 = t[t>-tag-t0][t[t>-tag-t0]< tag+t0]
        start = np.where(t == t1[0])[0][0]
        end = np.where(t == t1[-1])[0][0]
        ncfst[i][start:end+1]= ncfst[i][start:end+1]* np.exp(-a*np.abs((tag-np.abs(t1))))
    return ncfst
#%%
def linear_stack(key_subwork):
    global dir_project
    global key_subworks
    global info_basic
    global flag_save

    start0 = time.time()
    path_partition= dir_partition + str(key_subwork) + '.txt'
    stalist, lat_stations, lon_stations = np.loadtxt(path_partition, dtype='str', unpack=True)
    nsta = len(stalist)
    StationPairs = GetStationPairs(nsta)
    nPairs = int(len(StationPairs)/2)
    nf = info_basic['nf']
    ncfs_sum_linear = np.zeros((nPairs,nf),dtype=np.complex64)
    r = np.zeros(nPairs)
    f0 = info_basic_bi['f']
    count= np.zeros(nPairs)
    StationPairs = GetStationPairs(nsta)
    index = []
    for i in range(nPairs):
        sta1 = StationPairs[2*i]
        sta2 = StationPairs[2*i+1]
        idx1 = np.min( [int(stalist_all.index(stalist[sta1])),int(stalist_all.index(stalist[sta2]))] )
        idx2 = np.max( [int(stalist_all.index(stalist[sta1])),int(stalist_all.index(stalist[sta2]))] )
        
        m = 0
        for j in range(nsta_all-idx1,nsta_all):
            m += j
        num = m +idx2 - idx1 -1
        index.append(num)
        ncfs_sum_linear[i,:] = np.nan_to_num(ncfs[num,:])
        count[i] = count_all[num]
        r[i] = great_circle((lat_stations[sta1],lon_stations[sta1]),(lat_stations[sta2],lon_stations[sta2])).km
    if flag_save:
        tag = name_project[ name_project.rfind('_')+1: -1]
        dir_stack = dir_project + 'stack_'+tag+'/'
        if not os.path.exists(dir_stack):
            os.makedirs(dir_stack)
        info_basic['dir_stack'] = dir_stack
        outname = key_subwork+'_gather_linear.h5'
        if os.path.exists(dir_stack+outname):
            os.remove(dir_stack+outname)
        ncffile = h5py.File(dir_stack+outname,'w')
        ncffile.create_dataset('ncfs', data=ncfs_sum_linear) 
        ncffile.create_dataset('r',data=r)
        ncffile.create_dataset('count',data=count)
        ncffile.create_dataset('f',data=f0)
        ncffile.create_dataset('StationPairs',data=StationPairs)
        ncffile.close()
    print('Finish linear stack, time:', time.time()-start0, ' seconds')
    return ncfs_sum_linear, r, StationPairs, index

#%% remove stack
def remove_stack(key_subwork,ncfs_sum_linear,r,StationPairs):
    global dir_project
    global key_subworks
    global info_basic
    global flag_save
    global v_tag
    global t0
    global a

    start0 = time.time()
    f0 = info_basic_bi['f']
    dt = 1/np.max(f0)
    t = (np.linspace(-len(f0)-1,len(f0)-1,2*(len(f0)-1))+0.5)*dt/2
    ncfst_linear = mathlib.freq_time(ncfs_sum_linear)
    ncfst1 = time_window_filter(t,ncfst_linear,r,v_tag,t0,a)
    ncfs_sum_remove = mathlib.time_freq(ncfst1)
    if flag_save:
        tag = name_project[ name_project.rfind('_')+1: -1]
        dir_stack = dir_project + 'stack_'+tag+'/'
        outname = key_subwork+'_gather_timewindow.h5'
        if os.path.exists(dir_stack+outname):
            os.remove(dir_stack+outname)
        ncffile = h5py.File(dir_stack+outname,'w')
        ncffile.create_dataset('ncfs', data=ncfs_sum_remove) 
        ncffile.create_dataset('r',data=r)
        ncffile.create_dataset('StationPairs',data=StationPairs)
        ncffile.close()
    print('Finish remove stack, time:', time.time()-start0, ' seconds')
    return ncfs_sum_remove

# %%
c_min = 0.200
c_max = 2
c_num = 800
c = np.linspace(c_min,c_max,c_num)
info_basic['fj_c_min'] = c_min
info_basic['fj_c_max'] = c_max
info_basic['fj_c_num'] = c_num
with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

#%% parameter for stacking
v_tag = 2
a = 100
t0 = 0.01
info_basic['remove_v_tag'] = v_tag
info_basic['remove_t0'] = t0
info_basic['remove_a'] = a
with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

# %%
info_basic_bi['r_max'] = {} 
start00 = time.time()
for key_subwork in key_subworks:
    print("Collecting ",key_subwork,' ...')
    # linear stack
    ncfs_sum_linear, r, StationPairs,index = linear_stack(key_subwork)
    info_basic_bi['r_max'][key_subwork] = np.max(r)
    ncfs_sum_remove = remove_stack(key_subwork,ncfs_sum_linear,r,StationPairs)
    # FJ
    if os.path.exists(dir_ds+'ds_'+key_subwork+'.h5'):
        print(key_subwork+' exists')
        continue
    noise_fj(key_subwork,ncfs_sum_linear,ncfs_sum_remove,r,index)
    print('Finish *** '+ key_subwork +   '***; time since start:', time.time()-start00, ' seconds. Proceeded '+str(key_subworks.index(key_subwork)+1)+'/'+str(len(key_subworks))+' subworks.')
    

# %%
np.save(filename_bi,info_basic_bi)