# %%
import numpy as np
import h5py
import ccfj
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
import yaml
import os

# %%
import sys
sys.path.append(r'../')
from toollib_standard import mathlib
from toollib_standard import plotlib

# %%
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

# %%
name_project_based = info_basic['name_project_based']
dir_project_based = os.path.join(dir_project_workspace, name_project_based)
print('dir_project_based: ', dir_project_based)

# %%
filename = dir_project_based+'Basic_info.yml'
with open(filename, 'r', encoding='utf-8') as f:
    info_basic_based = yaml.load(f.read(), Loader=yaml.FullLoader)
filename_bi = dir_project+'Basic_info.npy'
info_basic_bi_based = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary

# %%
dir_image = dir_project+info_basic['dir_image']
dir_stack = dir_project+info_basic['dir_stack']

# %%
t = info_basic_bi['t']
f = info_basic_bi['f']
dt = 1/np.max(f)
t = (np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))+0.5)*dt/2

# %%
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

# %%
def cut(key_subwork):
    global flag_save 
    global info_basic
    global info_basic_bi
    global key_subworks
    global v_tag
    global t0
    global a
    global dir_project
    global dir_stack

    #t = info_basic_bi['t']
    f = info_basic_bi['f']
    dt = 1/np.max(f)
    t = np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))*dt/2
    t = (np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))+0.5)*dt/2

    inname = str(key_subwork)+'_gather_linear.h5'
    outname = str(key_subwork)+'_gather_timewindow.h5'
    ncffile_in = h5py.File(dir_stack + inname,'r')
    if os.path.exists(dir_stack+outname):
        os.remove(dir_stack+outname)
    ncffile = h5py.File(dir_stack+outname,'w')

    ncfs_sum_linear = ncffile_in['ncfs'][:]
    r = ncffile_in['r'][:]

    ncfst_linear = mathlib.freq_time(ncfs_sum_linear)
    ncfst1 = time_window_filter(t,ncfst_linear,r,v_tag,t0,a)
    ncfs_sum_remove = mathlib.time_freq(ncfst1)

    ncffile.create_dataset('ncfs', data=ncfs_sum_remove) 
    ncffile.create_dataset('r',data=r)
    ncffile.close()  
    ncffile_in.close()  

    print(str(key_subwork)+' is done! ' + str(key_subworks.index(key_subwork)+1) + '/' + str(len(key_subworks)))

# %%
def remove_zorocor(v_tag,t0,a,dir_project,key_subwork_sample = [],nThreads = 1):
    filename = dir_project+'Basic_info.yml'
    with open(filename, 'r', encoding='utf-8') as f:
        info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
    filename_bi = dir_project+'Basic_info.npy'
    info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary
    key_subworks = info_basic['key_subworks']

    if key_subwork_sample != []:
        key_subworks = key_subwork_sample

    if nThreads == 1:
        for key_subwork in key_subworks:
            cut(key_subwork)
    else:
        pool = ThreadPoolExecutor(max_workers = nThreads)
        for key_subwork in key_subworks:
            pool.submit(cut,key_subwork)
        pool.shutdown()

# %%
v_tag = info_basic['v_tag']
a = info_basic['a']
t0 = info_basic['t0']
# %%
with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

# %%
nThreads = 40
remove_zorocor(v_tag,t0,a,dir_project,nThreads=nThreads)
