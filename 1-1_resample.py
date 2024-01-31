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
import shutil
import pandas as pd
import yaml

#%%
nThreads = 1

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_SAC_workspace = dir_config['dir_SAC_workspace']
print('dir_SAC_workspace: ', dir_SAC_workspace)
# %%
name_SAC = '2014/'
name_SAC_new =  'data_resample_2/'
dir_SAC = os.path.join(dir_SAC_workspace,name_SAC)
dir_resample = os.path.join(dir_SAC_workspace,name_SAC_new)
if not os.path.exists(dir_resample):
    os.mkdir(dir_resample)

# %% [markdown]
# ### Parameters

# %%
stalistname = 'stations_info_all.xlsx'

y_start = 2014
y_end = 2015
#d_start = 132
d_start = 128
#d_end = 159
d_end = 159
fmax = 100                    # 降采样频率

# %%
def Checkdata(dirname):
    global day
    global dir_SAC
    filename = os.path.join(dir_SAC,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')
    if os.path.exists(filename):
        return True
    return False

# Set the rules for reading data
def Resample(i):
    global day
    global dir_SAC
    global dir_resample
    global namelist
    global fmax


    if i%200 == 0 and i != 0:
        print(i)

    dirname = namelist[i]
    filename = os.path.join(dir_SAC,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')
    st = obspy.read(filename)

    # 对每个通道的数据进行重采样
    st_resampled = st.copy().decimate(factor=5, strict_length=False)
    
    # 保存到dir_resample中
    st_resampled.write(os.path.join(dir_resample,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00'),format='SAC')

# %%
stainfo = pd.read_excel(stalistname)
nsta = len(stainfo.iloc[:,0])
StationPairs = GetStationPairs(nsta)
nPairs = int(len(StationPairs)/2)
stalist = stainfo['Station'].tolist()

# %%
days = []
years = set()
flag_d = 0
start = time.time()
for y in range(y_start,y_end):
    for d in range(d_start,d_end): 

        if not os.path.exists(dir_resample+str(d)):
            os.mkdir(dir_resample+str(d))

        year = str(y)
        day = "%03d"%d
        namelist = []
        for dirname in stalist:
            if Checkdata(dirname):
                namelist.append(dirname)
        nsta = len(namelist)
        # at least two stations needs
        if nsta >1:
            flag_d += 1
            years.add(int(year))
            days.append(d)
            print("year"+year+" day"+day +" nsta "+ str(nsta))
            # using multiThreads to read data
            pool = ThreadPoolExecutor(max_workers = nThreads)
            start0 = time.time()
            for i in range(nsta):
                pool.submit(Resample,i)
            pool.shutdown()
# %%



