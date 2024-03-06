# %%
import matplotlib.pyplot as plt
import numpy as np
import ccfj
import h5py
import time
import os
from concurrent.futures import ThreadPoolExecutor
import yaml

# %%
import sys
sys.path.append(r'../tools_F-J/')
from toollib_standard import maplib
from toollib_standard import mathlib
from toollib_standard import filelib
from toollib_standard import stacklib
from toollib_standard import plotlib

# %%
nThreads = 1

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
dir_image = dir_project+info_basic['dir_image']
dir_stack = dir_project+info_basic['dir_stack']
dir_ds = dir_project+info_basic['dir_ds']

#%%
def noise_fj(key_subwork):
    global dir_project
    global c
    global key_subworks
    global key_subworks_repick
    global dir_stack
    global dir_ds

    
    start0 = time.time()
    filename = dir_project+'Basic_info.yml'
    with open(filename, 'r', encoding='utf-8') as f:
        info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
    filename_bi = dir_project+'Basic_info.npy'
    info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary
    
    f = info_basic_bi['f']


    filename = dir_ds + 'ds_'+str(key_subwork)+'.h5'
    #print(filename)
    if key_subwork in key_subworks_repick:
        if os.path.exists(filename):
            os.remove(filename)
            print('Remove '+filename)
    if os.path.exists(filename):
        print(str(key_subwork)+' exists')
        return

    h5file = h5py.File(filename,'w')

    #print("F-J scan for "+filename+" stack  "+key_subwork)
    """
    data = np.load(dir_stack+key_subwork+'_summed-'+filename+'.npz')
    ncfs = data['ncfs']
    r = data['r']*1e3
    """
    # linear stack
    outname = str(key_subwork)+'_gather_linear.h5'
    ncffile = h5py.File(dir_stack + outname,'r')
    ncfs = ncffile['ncfs'][:]
    r = ncffile['r'][:]
    ncffile.close()
    
    # timewindow filtered stack
    outname = str(key_subwork)+'_gather_timewindow.h5'
    ncffile = h5py.File(dir_stack + outname,'r')
    ncfs_remove = ncffile['ncfs'][:]
    ncffile.close()
    


    #ds00 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=0,func=0)
    #ds01 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=1,func=0)
    #ds10 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=0,func=1)
    
    ds11 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=1,func=1)
    #ds = [ds10,ds11]
    #ds = np.array(ds)  
    ds = np.array(ds11).reshape(1,np.shape(ds11)[0],np.shape(ds11)[1])
    h5file.create_dataset('ds_linear',data=ds)

    
    ds11 = ccfj.fj_noise(np.real(ncfs_remove),r,c,f,fstride=1,itype=1,func=1)
    ds = np.array(ds11).reshape(1,np.shape(ds11)[0],np.shape(ds11)[1])
    h5file.create_dataset('ds_remove',data=ds)
    

    h5file.create_dataset('f',data=f)
    h5file.create_dataset('c',data=c)
    h5file.close()
    print('Finish '+ str(key_subwork) +   ' time:', time.time()-start0, ' seconds. Proceeded '+str(key_subworks.index(key_subwork)+1)+'/'+str(len(key_subworks))+' subworks.')

# %%
c_min = 0.200
c_max = 2
c_num = 800
c = np.linspace(c_min,c_max,c_num)
info_basic['c_min'] = c_min
info_basic['c_max'] = c_max
info_basic['c_num'] = c_num
with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

# %%
"""
pool = ThreadPoolExecutor(max_workers = nThreads)
for key_subwork in key_subworks:
    pool.submit(noise_fj,key_subwork)
pool.shutdown()
"""

# %%
flag_plot = 0
for key_subwork in key_subworks:
    #key_subwork = str(key_subwork)
    #print(1)
    flag_plot += 1
    #print(key_subwork)
    
    #start0 = time.time()
    noise_fj(key_subwork)

# %%
"""
if 'key_subworks_repick' in info_basic.keys():
    del info_basic['key_subworks_repick']
    with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
        yaml.dump(data=info_basic, stream=f, allow_unicode=True)
"""