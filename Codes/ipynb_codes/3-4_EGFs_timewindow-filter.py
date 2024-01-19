# %%
import numpy as np
import h5py
import ccfj
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor
import yaml
import time
import os

# %%
import sys
sys.path.append(r'../')
from toollib_standard import mathlib
from toollib_standard import plotlib

# %%
with open('a-project.yml', 'r', encoding='utf-8') as f:
    proj = yaml.load(f.read(), Loader=yaml.FullLoader)
proj_name = proj['name']
proj_name = '/shdisk/rem2/Harmon/F-J/San/project/output_FJSJ_17-03/'
proj_name

# %%
filename = proj_name+'Basic_info.yml'
with open(filename, 'r', encoding='utf-8') as f:
    info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
filename_bi = proj_name+'Basic_info.npy'
info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary

# %%
dir_stack = info_basic['dir_stack']
t = info_basic_bi['t']
f = info_basic_bi['f']
dt = 1/np.max(f)
t = (np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))+0.5)*dt/2

# %%
key_subworks = info_basic['key_subworks']
#key_subworks.remove('23-08')
#key_subworks = ['28-02']
"""
key_subworks = []
for key in info_basic['key_subworks']:
    num = int(key[3:])
    if num == 2:
        key_subworks.append(key)

key_subworks
"""

# %% [markdown]
# ### Gaussian time window

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
    global r

    time0 = time.time()
    dir_stack = info_basic['dir_stack']
    #t = info_basic_bi['t']
    f = info_basic_bi['f']
    dt = 1/np.max(f)
    t = np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))*dt/2
    t = (np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))+0.5)*dt/2

    inname = key_subwork+'_gather_linear.h5'
    outname = key_subwork+'_gather_timewindow.h5'
    ncffile_in = h5py.File(dir_stack + inname,'r')
    if os.path.exists(dir_stack+outname):
        os.remove(dir_stack+outname)
    ncffile = h5py.File(dir_stack+outname,'a')

    ncfs_sum_linear = ncffile_in['ncfs'][:]
    r = ncffile_in['r'][:]

    ncfst_linear = mathlib.freq_time(ncfs_sum_linear)
    ncfst1 = time_window_filter(t,ncfst_linear,r,v_tag,t0,a)
    ncfs_sum_remove = mathlib.time_freq(ncfst1)

    ncffile.create_dataset('ncfs', data=ncfs_sum_remove) 
    ncffile.create_dataset('r',data=r)
    ncffile.close()  
    ncffile_in.close()  
    print(key_subwork+' is done! ' + str(key_subworks.index(key_subwork)+1) + '/' + str(len(key_subworks))+ ' time: '+str(time.time()-time0))

# %%
def remove_zorocor(v_tag,t0,a,proj_name,key_subwork_sample = [],nThreads = 1):
    filename = proj_name+'Basic_info.yml'
    with open(filename, 'r', encoding='utf-8') as f:
        info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
    filename_bi = proj_name+'Basic_info.npy'
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
v_tag = 2
a = 100
t0 = 0.01
#tag = 0.00055
#a = 0.5
key_subwork_sample = [key_subworks[0]]
remove_zorocor(v_tag,t0,a,proj_name,key_subwork_sample)

# %%
info_basic['v_tag'] = v_tag
info_basic['t0'] = t0
with open(proj_name+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

# %%
key_subwork = key_subwork_sample[0]
outname = key_subwork+'_gather_linear.h5'
ncffile = h5py.File(dir_stack + outname,'r')
outname = key_subwork+'_gather_timewindow.h5'
ncffile_2 = h5py.File(dir_stack + outname,'r')
ncfs_sum_linear_remove = ncffile_2['ncfs'][:]
ncfst_linear_remove = mathlib.freq_time(ncfs_sum_linear_remove)
ncfs_sum_linear = ncffile['ncfs'][:]
ncfst_linear = mathlib.freq_time(ncfs_sum_linear)
r = ncffile['r'][:]
ncffile.close()
ncffile_2.close()

# %%
fig,ax = plt.subplots(1,2,figsize = (20,10))

start = info_basic['start']
interval = info_basic['interval']
flag_time = info_basic['flag_time']
xlim_T = [-1,1]
title1 = 'Linear'
ax[0] = plotlib.plot_ncfst(ax[0],t,ncfst_linear[start::interval],r[start::interval],title1,flag_time,xlim_T,0)
#ax[0].vlines(tag,0,0.2,linestyles='dashed',colors='g')
#ax[0].vlines(-tag,0,0.2,linestyles='dashed',colors='g')
ax[1] = plotlib.plot_ncfst(ax[1],t,ncfst_linear_remove[start::interval],r[start::interval],title1,flag_time,xlim_T,0)
#ax[1].vlines(tag,0,0.2,linestyles='dashed',colors='g')
#ax[1].vlines(-tag,0,0.2,linestyles='dashed',colors='g')
for i in range(0,len(r),interval):
    ax[0].plot(t0+r[i]/v_tag,r[i],'r.')
    ax[0].plot(-t0-r[i]/v_tag,r[i],'r.')
    ax[1].plot(t0+r[i]/v_tag,r[i],'r.')
    ax[1].plot(-t0-r[i]/v_tag,r[i],'r.')
plt.tight_layout()
plt.show()

# %%
nThreads = 1
remove_zorocor(v_tag,t0,a,proj_name,nThreads=nThreads)

# %%



