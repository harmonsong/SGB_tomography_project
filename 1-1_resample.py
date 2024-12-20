# %%
import sys, getopt
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
nThreads = 20
fmax = 250                    # 降采样频率

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_SAC_workspace = dir_config['dir_SAC_workspace']
dir_SAC_workspace = '/share/data/California/SJFZ_nodal/'
print('dir_SAC_workspace: ', dir_SAC_workspace)
# %%
name_SAC = '2014/'
name_SAC_new =  'resample_'+str(fmax)+'Hz/'
dir_SAC = os.path.join(dir_SAC_workspace,name_SAC)
dir_resample = 'data_resample/'+name_SAC_new
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
d_len = 31
#d_end = 159
Fs = 500
factor = int(Fs/fmax)

# %%
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:l:n:", ["start=", "len=", "nThreads="])
except getopt.GetoptError:
    print('1-2-2_correlation_day.py -s <day start> -l <day len>')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-s", "--start"):
        d_start = int(arg)
    elif opt in ("-l", "--len"):
        d_len = int(arg)
    elif opt in ("-n", "--nThreads"):
        nThreads = int(arg)

d_end  = d_start + d_len
print('d_start: ', d_start)
print('d_end: ', d_end)

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
    global factor

    dirname = namelist[i]
    outname = os.path.join(dir_resample,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')


    if i%200 == 0 and i != 0:
        print(i)
    
    filename = os.path.join(dir_SAC,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')
    st = obspy.read(filename)

    # 对每个通道的数据进行重采样
    st_resampled = st[0].copy().decimate(factor=factor, strict_length=False)
    
    # 保存到dir_resample中
    st_resampled.write(outname,format='SAC')

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
        print(namelist)
        # at least two stations needs
        if nsta >1:
            flag_d += 1
            years.add(int(year))
            days.append(d)
            print("year"+year+" day"+day +" nsta "+ str(nsta))
            # using multiThreads to read data
            if nThreads > 1:
                pool = ThreadPoolExecutor(max_workers = nThreads)
                start0 = time.time()
                for i in range(nsta):
                    pool.submit(Resample,i)
                pool.shutdown()
            else:
                for i in range(nsta):
                    Resample(i)
# %%



