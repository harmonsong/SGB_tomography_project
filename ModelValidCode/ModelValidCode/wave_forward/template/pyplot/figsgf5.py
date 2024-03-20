#!/usr/bin/env python3
# -*- conding: utf-8 -*-
"""
First, synthesize the noise cross-correlation functions (ccfs) by convoluatating the source time function and the Green's function obtained by 3-D FD wave simulation code; Then, compare observed ccfs with synthetic ccfs within several period windows.

Created on Fri Jun 24 10:48 2022
@author: Juqing Chen, USTC
"""

import sys, getopt, os.path
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import math,h5py, pandas
from scipy import integrate
from scipy.fftpack import fft, hilbert
from scipy.signal import butter, lfilter
from geopy.distance import great_circle
#========================================================================================

def getSNR(cfsf,r,f,SNRmin):
    NoiseWinLength = 300 # in seconds
    SigWinC = np.array([2,4]) #signal window for velocity, in km/s
    SigWinT = np.zeros((len(r),2))
    SNR = np.zeros(np.shape(r))
    for i in range(len(r)):
        SigWinT[i,:] = r[i]/SigWinC

    dt = 1/np.max(f)
    t = (np.linspace(-len(f)/2,len(f)/2-1,len(f))-0.5)*dt
    cfst = np.zeros(np.shape(cfsf))
    for i in range(len(cfsf)):
        cfst[i,:] = np.real(np.fft.fftshift(np.fft.ifft(cfsf[i,:])))
        cfst[i,:] = cfst[i,:] /np.max(cfst[i,:])

    id0 = int(len(t)/2)
    cfst12 = cfst[:,id0::]
    cfst21 = cfst[:,id0::-1]
    for j in range(len(r)):
        c12h= hilbert(cfst12[j,:])
        c21h= hilbert(cfst21[j,:])
        groupCF12 = np.sqrt(cfst12[j,:]**2 + c12h**2)
        groupCF21 = np.sqrt(cfst21[j,:]**2 + c21h**2)
        nnSigWinL = int(np.ceil(SigWinT[j,1]/dt))
        nnSigWinR = int(np.floor(SigWinT[j,0]/dt))
        if nnSigWinL >= nnSigWinR:
            continue
        nnNoiseWinL = nnSigWinR
        nnNoiseWinR = nnSigWinR+int(round(NoiseWinLength/dt))

        #print("nnSigWinL: ",nnSigWinL)
        #print("nnSigWinR: ",nnSigWinR)
        SigAmp12 = np.max(groupCF12[nnSigWinL:nnSigWinR])
        SigAmpAve12 = np.mean(groupCF12[nnSigWinL:nnSigWinR])
        NoiseAmpAve12 = np.mean(groupCF12[nnNoiseWinL:nnNoiseWinR])
        if NoiseAmpAve12 == 0:
            if SigAmp12 > 0:
                SNR12 = 100
            else:
                SNR12 = 0
        else:
            SNR12 = SigAmp12/NoiseAmpAve12

        SigAmp21 = np.max(groupCF21[nnSigWinL:nnSigWinR])
        SigAmpAve21 = np.mean(groupCF21[nnSigWinL:nnSigWinR])
        NoiseAmpAve21 = np.mean(groupCF21[nnNoiseWinL:nnNoiseWinR])
        if NoiseAmpAve21 == 0:
            if SigAmp21 > 0:
                SNR21 = 100
            else:
                SNR21 = 0
        else:
            SNR21 = SigAmp21/NoiseAmpAve21

        SNR[j] = np.nan_to_num(min(SNR12,SNR21))
        
    idx = []
    for k in range(len(r)):
        if SNR[k] < np.mean(SNR)*SNRmin:
            idx.append(k)

    return idx

def ricker(t,fc,t0):
 if (t<=0.0):
    v=0.0
 f0 = np.sqrt(math.pi)/2.0
 u = (t-t0)*2.0*math.pi*fc
 v=(u**2.0/4.0-0.5)*np.exp(-u**2.0/4.0)*f0
 return v

def butter_bandpass(lowcut, highcut, fs, order=2):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandpass')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def print_help_info():
  print('\nUsage: fig4seismo2d.py [-h] | [-i inputpath]' +
    ' [-o outputpath] [-c fileFDconf] [-s isrc] [-d stride]')
  print('  [inputpath  ]: directory of the coordinate data files.')
  print('  [outputpath ]: directory of the seismogram data files.')
  print('  [fileFDconf   ]: name of the configure file for FD.')
  print('  [isrc   ]: code number of the virtrual source')
  print('  [stride   ]: number of stations to skip for plotting.')

#==============================================================================

putinrou = ['/home/cjq/ModelValid/Forward/CCF/CFs/SCB_3600_0.9_summed.h5']
stalist = '/home/cjq/ModelValid/Forward/SCB_locat_v2.txt'
srclist = '/home/cjq/ModelValid/Forward/virsrc.txt'
inputpath = '/home/cjq/ModelValid/Forward/template/model_Chen/input/'
outputpath = 'output/'
fileFDconf = 'SeisFD3D.conf'
stride = 1
isline = False
LenFD = 3
isrc = 1

#====================================================================================   
try:
  opts, args = getopt.getopt(sys.argv[1:], 'hl:i:o:c:s:d:', ['help',
    'inputpath=', 'outputpath=', 'fileFDconf=',
    'isrc', 'stride='])
except getopt.GetoptError:
  print_help_info()
  sys.exit(2)

for opt, arg in opts:
  if opt in ('-h', '--help'):
    print_help_info()
    sys.exit()
  elif opt in ('-i', '--inputpath'):
    inputpath = arg
  elif opt in ('-o', '--outputpath'):
    outputpath  = arg
  elif opt in ('-c', '--fileFDconf'):
    fileFDconf = arg
  elif opt in ('-s', '--isrc'):
    isrc = int(arg)
  elif opt in ('-d', '--stride'):
    stride = int(arg)

print('\nConfigure:')
print('  inputpath:\t', inputpath)
print('  outputpath:\t', outputpath)
print('  fileFDconf:\t', fileFDconf)
print('  isrc:\t', isrc)
print('  stride:\t', stride)
#=============================read virtual sources=================================
net_sta = []
sta_profile = pandas.read_table(srclist, header=None, sep='\\s+')
net = list(sta_profile[0][:])
sta = list(sta_profile[1][:])
for i in range(len(net)):
    net_sta.append(str(net[i] + '-' + sta[i]))
srcsta = net_sta[isrc-1]
print('  srcsta:\t', srcsta)
#=============================read observed Ccfs====================================

CC_file = []
CC_prof = []
for i in range(len(putinrou)):
  CC_prof.append(pandas.read_hdf(putinrou[i], 'cc_prof'))
  with h5py.File(putinrou[i], mode = "r") as file_rou:
      CC_file.append(np.array(file_rou[list(file_rou.keys())[0]]))
CC_array = CC_file[0]
CC_table = CC_prof[0]
f_prof = pandas.read_hdf(putinrou[0], 'freq')
freq = np.array(f_prof[:])[:,0]
#==============selec Ccfs based on srcstat=====================================
part_CC = []
part_dis = []
cc_name = CC_table['0cc_name']
for i in range(len(cc_name)):
  stas = cc_name[i].split('_')
  if srcsta in stas:
     part_CC.append((CC_array[i,:]).copy())
     part_dis.append(CC_table.iloc[i,3]) #km

part_CC = np.array(part_CC)
part_dis = np.array(part_dis).copy()
dele_cc_rows = []
for CC_num in range(part_CC.shape[0]):
  if part_CC[CC_num,:].all() == 0:
      dele_cc_rows.append(CC_num)
  if part_dis[CC_num] == 0 or part_dis[CC_num]>= 800:
      dele_cc_rows.append(CC_num)
dele_cc_rows = list(set(dele_cc_rows))
raw_rou_nonezerocc = np.delete(part_CC, dele_cc_rows, axis = 0)
dis_vin_nonezerocc = np.delete(part_dis, dele_cc_rows, axis = 0)
idel = getSNR(raw_rou_nonezerocc,dis_vin_nonezerocc,freq,SNRmin=0.8)
raw_rou_nonezerocc = np.delete(raw_rou_nonezerocc,idel,axis=0)
dis_vin_nonezerocc = np.delete(dis_vin_nonezerocc,idel,axis=0)
#============obtain time-domain Green's function=================================================
ncfs0 = raw_rou_nonezerocc
robs = dis_vin_nonezerocc #km
dtobs = 1/np.max(freq)
tobs = (np.linspace(-len(freq)/2, len(freq)/2-1, len(freq))-0.5)*dtobs
ncfst = np.zeros(np.shape(ncfs0))
for i in range(len(ncfs0)):
  ncfst[i,:] = np.real(np.fft.fftshift(np.fft.ifft(ncfs0[i,:])))

cfs = (ncfst[:,int(len(freq)/2):1:-1] + ncfst[:,int(len(freq)/2+1)::]) / 2 
cfs_grad = -np.gradient(cfs,axis=1,edge_order=1)/dtobs
#===========convolve with the ricker wavelet=============================================================
tt = tobs[int(len(freq)/2+1)::]
wavelet = np.zeros(len(tt))
waveobs = np.zeros((cfs_grad.shape[0],len(tt)))
for i in range(len(tt)):
   wavelet[i] = ricker(tt[i],0.1,15)
indx_obs = int(15/dtobs)
for j in range(cfs_grad.shape[0]):
   tmp = np.convolve(wavelet,cfs_grad[j,:])
   waveobs[j,:] = tmp[indx_obs:len(tt)+indx_obs]
#===========windowed waveform=================================================
fs_obs = 1.0/(tt[1]-tt[0])
bands = np.array([[3,7],[5,10],[8,15],[12,20]])
aver = np.array([3.12,3.12,3.12,3.30])
lowcut = 1.0 / bands[:,1]
highcut = 1.0 / bands[:,0]
wave0obs = np.zeros(waveobs.shape) #3-7s
wave1obs = np.zeros(waveobs.shape) #5-10s
wave2obs = np.zeros(waveobs.shape) #8-15s
wave3obs = np.zeros(waveobs.shape) #12-20s

for i in range(waveobs.shape[0]):
    wave0obs[i,:] = butter_bandpass_filter(waveobs[i,:], lowcut[0], highcut[0], fs_obs, order=2)
    wave1obs[i,:] = butter_bandpass_filter(waveobs[i,:], lowcut[1], highcut[1], fs_obs, order=2)
    wave2obs[i,:] = butter_bandpass_filter(waveobs[i,:], lowcut[2], highcut[2], fs_obs, order=2)
    wave3obs[i,:] = butter_bandpass_filter(waveobs[i,:], lowcut[3], highcut[3], fs_obs, order=2)

#===================Process synthetic waveforms==========================================================
confFD = open(fileFDconf, 'r')
linesFD = confFD.read().split('\n')
for j in range(len(linesFD)):
  line = linesFD[j]
  if line.find('#') >= 0: line = line[:line.find('#') - 1]
  if line.split(' ')[0] == 'dims':
    dims = [int(v) for v in line.split() if v.isdigit()][0:3]

confFD.close()

staloc = np.loadtxt(stalist, usecols=[2,3])
netsta = []
with open(stalist, 'r') as f:
    while True:
        tmp = f.readline()
        if tmp:
           netsta.append(tmp.split()[0] + '-' + tmp.split()[1])
        else:
            break
Srloc = staloc[netsta.index(srcsta)]
#==============================================================================
nsta = staloc.shape[0]
r0 = np.zeros(nsta)
for k in range(nsta):
  r0[k] =great_circle((Srloc[1], Srloc[0]),(staloc[k,1],staloc[k,0])).km
indx = np.argsort(r0)
rsyn = r0[indx]
#==============================================================================

whitems = ['Vx', 'Vy', 'Vz', 'Txx', 'Tyy', 'Tzz', 'Txy', 'Txz', 'Tyz']

isfirstget = True
isgottime = False
iline = 0
for k in range(dims[2]):
  for j in range(dims[1]):
    for i in range(dims[0]):
      filename = '{}seismo_mpi{:02d}{:02d}{:02d}.nc'.format(outputpath, i, j, k)
      if not os.path.isfile(filename):
        continue
      ssm = nc.Dataset(filename, 'r')
      filename = '{}station_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, i, j, k)
      stn = nc.Dataset(filename, 'r')
      pid = stn.variables['id'][:]#get the station index(two dimension)
      if not isgottime:
        t = ssm.variables['time'][:]
        isgottime = True
        nt = t.size # the acctural points stored in the Seis.nc
      for ip in range(pid.shape[0]):
        if pid[ip, 1] == iline:#equals to the L/P-index(P=0) ////get id from STA.nc ///to confirm line number
          ttvar = np.zeros((1, 1, nt))
          ttvar[0, 0, :] = np.reshape(ssm.variables[whitems[0]][:, ip], [nt]) #get the ip strip in Vx, and flaten into [nt] case
          for iw in range(1,9):
            tvar = np.reshape(ssm.variables[whitems[iw]][:, ip], [1, 1, nt])
            ttvar = np.append(ttvar, tvar, axis = 1) #all stored in falten case and different in layerd direction
          if isfirstget:
            dvar = ttvar[:] #data  dvar[Var][Var][nt]
            ipoint = pid[ip, 0] # relevent index(start from 1)
            isfirstget = False
          else:
            dvar = np.append(dvar, ttvar, axis = 0)
            ipoint = np.append(ipoint, pid[ip, 0])
      stn.close()
      ssm.close()
#ipoint is absolute point location in point series(sequence in line).
if 'ipoint' not in globals():
  raise IOError('Maybe no seismogram data files in the outputpath: ' + outputpath)

# nt point in time series
# npt point in record index series

npt = ipoint.size
pvar = np.zeros((npt, 9, nt)) #RECiver Var Time
for ip in range(npt):
  pvar[ipoint[ip] - 1, :, :] = dvar[ip, :, :] 
for iw in range(9):
  if iw < 3:
    pvar[:, iw, :] = pvar[:, iw, :]/1.0e3
  else:
    pvar[:, iw, :] = pvar[:, iw, :]/1.0e9

print(np.shape(pvar))

# Integrate to obtain Displacement
Vz = pvar[indx, 2, :]
Dz0 = integrate.cumtrapz(Vz, t, axis=1, initial=0)
stpindx = int(15/(t[1]-t[0]))
tc = t[:len(t)-stpindx]
Dz = Dz0[:,stpindx:]


# filter-----------------------------------------------------------
fs_syn = 1.0/(t[1]-t[0])
wave0syn = np.zeros(Dz.shape) #3-7s
wave1syn = np.zeros(Dz.shape) #5-10s
wave2syn = np.zeros(Dz.shape) #8-15s
wave3syn = np.zeros(Dz.shape) #12-20s

for i in range(Dz.shape[0]):
    wave0syn[i,:] = butter_bandpass_filter(Dz[i,:], lowcut[0], highcut[0], fs_syn, order=2)
    wave1syn[i,:] = butter_bandpass_filter(Dz[i,:], lowcut[1], highcut[1], fs_syn, order=2)
    wave2syn[i,:] = butter_bandpass_filter(Dz[i,:], lowcut[2], highcut[2], fs_syn, order=2)
    wave3syn[i,:] = butter_bandpass_filter(Dz[i,:], lowcut[3], highcut[3], fs_syn, order=2)
#==============correct synthesized waveform with observed waveform based on station interval==========================================================
indxcorr = []
for i in range(len(robs)):
  for j in range(len(rsyn)):
    if abs(int(robs[i]*1e3)-int(rsyn[j]*1e3)) == 0:
       indxcorr.append(j)
       break
wave0sync = wave0syn[indxcorr]
wave1sync = wave1syn[indxcorr]
wave2sync = wave2syn[indxcorr]
wave3sync = wave3syn[indxcorr]
#print("wave0syn: ", wave0syn.shape)
#print("wave0sync: ", wave0sync.shape)
#print("wave0obs: ", wave0obs.shape)
#==============================================================================
# -- Dz  --
xlimt = 300 #s
linew = 0.6
fonts = 12
fig = plt.figure(figsize = (8, 6), dpi = 80)
maxlen = min(len(robs), wave0sync.shape[0])
for j in range(0,maxlen,stride):
  if robs[j] < bands[0,1]*aver[0]:
     continue
  t0 = robs[j]/aver[0]-5
  t1 = robs[j]/aver[0]+18
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave0obs[j,:indx0] = 0
  wave0obs[j,indx1:] = 0
  pt,= plt.plot(tt, wave0obs[j,:]/np.max(abs(wave0obs[j,:]))*10+robs[j],'b', antialiased='False',linewidth=linew)
  label = "Observation"
  pt0,= plt.plot(tc, wave0sync[j]/np.max(abs(wave0sync[j]))*10+robs[j], 'r', antialiased='False', linewidth=linew)
  label0 = "This study"
plt.ylim([0, 800])
plt.xlim([0, xlimt])
plt.xlabel('Time (s)')
plt.ylabel('Distance (km)')
plt.legend([pt, pt0],[label, label0], fontsize=fonts, loc=4)
plt.title('Waveform comparison for 3-7s')
plt.savefig('comparison_3_7s.jpg', dpi=400)

fig = plt.figure(figsize = (8, 6), dpi = 80)
for j in range(0,maxlen,stride):
  if robs[j] < bands[1,1]*aver[1]:
     continue
  t0 = robs[j]/aver[1]
  t1 = t0+23
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave1obs[j,:indx0] = 0
  wave1obs[j,indx1:] = 0
  pt,= plt.plot(tt, wave1obs[j,:]/np.max(abs(wave1obs[j,:]))*10+robs[j],'b', antialiased='False',linewidth=linew)
  label = "Observation"
  pt1,= plt.plot(tc, wave1sync[j]/np.max(abs(wave1sync[j]))*10+robs[j], 'r', antialiased='False', linewidth=linew)
  label1 = "This study"
plt.ylim([0, 800])
plt.xlim([0, xlimt])
plt.xlabel('Time (s)')
plt.ylabel('Distance (km)')
plt.legend([pt, pt1],[label, label1], fontsize=fonts, loc=4)
plt.title('Waveform comparison for 5-10s')
plt.savefig('comparison_5_10s.jpg', dpi=400)

fig = plt.figure(figsize = (8, 6), dpi = 80)
for j in range(0,maxlen,stride):
  if robs[j] < bands[2,1]*aver[2]:
     continue
  t0 = robs[j]/aver[2]
  t1 = t0 + 30
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave2obs[j,:indx0] = 0
  wave2obs[j,indx1:] = 0
  pt,= plt.plot(tt, wave2obs[j,:]/np.max(abs(wave2obs[j,:]))*10+robs[j],'b', antialiased='False',linewidth=linew)
  label = "Observation"
  pt2,= plt.plot(tc, wave2sync[j]/np.max(abs(wave2sync[j]))*10+robs[j], 'r', antialiased='False', linewidth=linew)
  label2 = "This study"
plt.ylim([0, 800])
plt.xlim([0, xlimt])
plt.xlabel('Time (s)')
plt.ylabel('Distance (km)')
plt.legend([pt, pt2],[label, label2], fontsize=fonts, loc=4)
plt.title('Waveform comparison for 8-15s')
plt.savefig('comparison_8_15s.jpg', dpi=400)

fig = plt.figure(figsize = (8, 6), dpi = 80)
for j in range(0,maxlen,stride):
  if robs[j] < bands[3,1]*aver[3]:
     continue
  t0 = robs[j]/aver[3]
  t1 = t0 + 45
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave3obs[j,:indx0] = 0
  wave3obs[j,indx1:] = 0
  pt,= plt.plot(tt, wave3obs[j,:]/np.max(abs(wave3obs[j,:]))*10+robs[j],'b', antialiased='False',linewidth=linew)
  label = "Observation"
  pt3,= plt.plot(tc, wave3sync[j]/np.max(abs(wave3sync[j]))*10+robs[j], 'r', antialiased='False', linewidth=linew)
  label3 = "This study"
plt.xlabel('Time (s)')
plt.ylabel('Distance (km)')
plt.ylim([0, 800])
plt.xlim([0, xlimt])
plt.legend([pt, pt3],[label, label3], fontsize=fonts, loc=4)
plt.title('Waveform comparison for 12-20s')
plt.savefig('comparison_12_20s.jpg', dpi=400)

#plt.show()
