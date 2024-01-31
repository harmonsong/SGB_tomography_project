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
nThreads = 240                # 线程数

# %%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_SAC_workspace = dir_config['dir_SAC_workspace']
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_SAC_workspace: ', dir_SAC_workspace)
print('dir_CC_workspace: ', dir_CC_workspace)

# %%
name_SAC = 'data_resample/'
#name_CC = 'CC/CC_40_onebit/'
name_CC = 'CC/CC_40_prewhiten-onebit/'
dir_SAC = os.path.join(dir_SAC_workspace,name_SAC)
dir_CC = os.path.join(dir_CC_workspace,name_CC)
if not os.path.exists(dir_CC):
    os.mkdir(dir_CC)
info = {}

# %% [markdown]
# ### Parameters

# %%
stalistname = 'stations_info_all.xlsx'
y_start = 2014
y_end = 2015


d_start = 128
d_end = 159
d_start1 = d_start
d_end1 = d_end
#d_start1 = 128
#d_end1 = 129


d_len = 32
#Fs = 500
Fs = 100

fmax = 40                    # 降采样频率
fftlen = Fs*60*5            # 用于做户相关的时间窗长度
nf = fmax*10                 # 输出的户相关频点数
fstride = fmax*fftlen/nf/Fs
f = np.arange(0, nf)*Fs/fftlen*fstride
dt = 1/np.max(f)
t = (np.linspace(-len(f)-1,len(f)-1,2*(len(f)-1))+0.5)*dt/2

overlaprate = 0.9

flag_onebit = 1             # 是否进行onebit
flag_prewhiten = 1          # 是否进行prewhiten

segday = 1                  # 切割的时间段，未实现，只能为1
npts = int(60*60*24*Fs/segday)  # total time (day)



# %%
info['dir_SAC'] = name_SAC
info['stalistname'] = stalistname
info['Fs'] = Fs
info['fmax'] = fmax
info['fftlen'] = fftlen
info['nf'] = nf
info['fstride'] = fstride
info['f'] = f
info['t'] = t
info['overlaprate'] = overlaprate
info['nThreads'] = nThreads
info['segday'] = segday
info['npts'] = npts
info['flag_onebit'] = flag_onebit
info['flag_prewhiten'] = flag_prewhiten
info['y_start'] = y_start
info['y_end'] = y_end
info['d_start'] = d_start
info['d_end'] = d_end
info['d_len'] = d_len

# %% [markdown]
# ### Correlation day by day

# %%
def Checkdata(dirname):
    global day
    global dir_SAC
    filename = os.path.join(dir_SAC,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')
    if os.path.exists(filename):
        return True
    return False


# Set the rules for reading data
# this will determined by the data files' distribution
def Readdata(i):
    global day
    global dir_SAC
    global npts
    global namelist
    global Fs
    global data
    global startend


    if i%200 == 0 and i != 0:
        print(i)

    dirname = namelist[i]
    filename = os.path.join(dir_SAC,str(day),dirname+'.EHZ.2014.'+str(day)+'.00.00.00')
    st = obspy.read(filename)
    st[0].detrend(type='constant')

    st[0].detrend('linear')
    st[0].detrend('demean')

    st[0].detrend("spline", order=3, dspline=400)
    if st[0].stats.npts >= npts:
        data[npts*i:npts*(i+1)] = st[0].data[0:npts]
        startend[i*2] = 0
        startend[i*2+1] = npts
    else:
        t0 = UTCDateTime(st[0].stats.starttime)
        idx0 = int((t0.second+t0.minute*60+t0.hour*60*60)*Fs)
        t1 = UTCDateTime(st[0].stats.endtime)
        idx1 = int((t1.second+t1.minute*60+t1.hour*60*60)*Fs)
        #print(i,idx0,idx1)
        data[(npts*i+idx0):(npts*i+idx1)] = st[0].data[0:(idx1-idx0)]
        startend[i*2] = idx0
        startend[i*2+1] = idx1

# %%
stainfo = pd.read_excel(stalistname)
nsta = len(stainfo.iloc[:,0])
StationPairs = GetStationPairs(nsta)
nPairs = int(len(StationPairs)/2)
stalist = stainfo['Station'].tolist()

# %%
filename = os.path.join(dir_CC+'info_CC.npy')
np.save(filename,info)

# %%
days = []
years = set()
flag_d = 0
start = time.time()
for y in range(y_start,y_end):
    for d in range(d_start1,d_end1): 
        for this_seg in range(segday):
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
                data = np.zeros([npts*nsta],dtype = np.float32)
                startend = np.zeros([nsta*2],dtype = np.int32)
                # using multiThreads to read data
                pool = ThreadPoolExecutor(max_workers = nThreads)
                start0 = time.time()
                for i in range(nsta):
                    pool.submit(Readdata,i)
                pool.shutdown()
                ## Generate a list for Station Paris 
                StationPairs = GetStationPairs(nsta)
                nPairs = int(len(StationPairs)/2)
                ## Crosscorrelation
                print('finish reading data, start crosscorrelation', time.time()-start0,' seconds')
                start0 = time.time()
                ncfs = CC(npts, nsta, nf, fftlen, StationPairs, startend, data,
                          fstride=fstride, overlaprate=overlaprate, nThreads=nThreads,ifonebit=flag_onebit,ifspecwhittenning=flag_prewhiten)
                print('finish crosscorrelation, start saving data', time.time()-start0,' seconds')
                outname = os.path.join(dir_CC,year+'-'+day+'.npz')
                if os.path.exists(outname):
                    os.remove(outname)
                np.savez(outname,ncfs=ncfs,StationPairs=StationPairs,stalist=namelist)
                #print('time:',time.time()-start0,' seconds')
                data = []
                ncfs = []
