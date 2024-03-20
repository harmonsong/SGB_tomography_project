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
import math,h5py
import pandas as pd
from scipy import integrate, interpolate, signal
from scipy.fftpack import fft, hilbert
from scipy.signal import butter, lfilter
from geopy.distance import great_circle
#========================================================================================
def ray_cross_path(stloc, edloc, dx):
    #---------------------------------------
    # stloc: start location (lon, lat)
    # edloc: end location (lon, lat)
    # dx: interval for longitude and latitude 
    #----------------------------------------

    #print("stloc: ", stloc)
    #print("edloc: ", edloc)
    xyn = []
    X0 = 108; Y0 = 20
    idx1 = int(np.ceil((stloc[0]-X0) / dx))
    idx2 = int(np.ceil((edloc[0]-X0) / dx))
    idy1 = int(np.ceil((stloc[1]-Y0) / dx))
    idy2 = int(np.ceil((edloc[1]-Y0) / dx))

    xyn.append([idx1, idy1])
    if (edloc[0] == stloc[0]):
        for idy in range(min(idy1, idy2), max(idy1, idy2)+1):
            #print("idx1: ", idx1)
            #print("idy: ", idy)
            xyn.append([idx1, idy])
        return  np.array(xyn)
    #print("idx1: ", idx1)
    #print("idx2: ", idx2)
    #print("idy1: ", idy1)
    #print("idy2: ", idy2)

    if (edloc[1] <= stloc[1]):
        stepy = -1
        iystart = idy1 - 1
        iyend = idy2 -1
    else:
        stepy = 1
        iystart = idy1
        iyend = idy2

    if (edloc[0] <= stloc[0]):
        stepx = -1
        ixstart = idx1 - 1
        ixend = idx2 - 1
    else:
        stepx = 1
        ixstart = idx1
        ixend = idx2

    #print("iystart: ", iystart)
    #print("iyend: ", iyend)
    #print("stepy: ", stepy)

    k = (edloc[1]-stloc[1]) / (edloc[0]-stloc[0])
    b = edloc[1] - k*edloc[0]

    for iy in range (iystart, iyend, stepy):
        #print("iy: ", iy)
        yh = iy*dx + Y0
        xv = (yh - b) / k
        ixl = int(np.floor((xv-X0) / dx))
        ixr = int(np.ceil((xv-X0) / dx))
        #print("ixr: ", ixr)
        if [ixr, iy] not in xyn:
             xyn.append([ixr, iy])

    for ix in range(ixstart, ixend, stepx):
        xv = ix*dx + X0
        yh = k*xv + b
        iyu = int(np.ceil((yh - Y0) / dx))
        iyd = int(np.floor((yh - Y0) / dx))
        if [ix, iyu] not in xyn:
             xyn.append([ix, iyu])
        if ((yh-(iyd*dx+Y0) > 0.01) and (iyu*dx + Y0-yh>0.01)):
            ii = ix+1
            if [ii, iyu] not in xyn:
                 xyn.append([ii, iyu])

    return np.array(xyn)

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
syn_stalist = '/home/cjq/ModelValid/Forward/SCB_locat_v2.txt'
sele_stalist = '/home/cjq/ModelValid/Forward/SCB_locat_v2.txt'
srclist = '/home/cjq/ModelValid/Forward/virsrc.txt'
inputpath = '/home/cjq/ModelValid/Forward/template/model_Chen/input/'
outputpath = 'output/'
fileFDconf = 'SeisFD3D.conf'
stride = 1
isline = False
LenFD = 3
isrc = 1

save_cell_3_7s = 'cell_3_7s.txt'
save_cell_5_10s = 'cell_5_10s.txt'
save_cell_8_15s = 'cell_8_15s.txt'
save_cell_12_20s = 'cell_12_20s.txt'
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
#=============================Read station location for all stations=================================
staloc = np.loadtxt(syn_stalist, usecols=[2,3])
syn_netsta = []
with open(syn_stalist, 'r') as f:
    while True:
        tmp = f.readline()
        if tmp:
           syn_netsta.append(tmp.split()[0] + '-' + tmp.split()[1])
        else:
            break
#=============================read virtual sources=================================
virsta = []
sta_profile = pd.read_table(srclist, header=None, sep='\\s+')
net = list(sta_profile[0][:])
sta = list(sta_profile[1][:])
srcloc = np.loadtxt(srclist, usecols=(2,3))
for i in range(len(sta)):
    virsta.append(str(net[i] + '-' + sta[i]))
srcsta = virsta[isrc-1]
print('  srcsta:\t', srcsta)
#=============================read observed Ccfs====================================
sele_netsta = []
with open(sele_stalist, 'r') as f:
    while True:
        tmp = f.readline()
        if tmp:
            sele_netsta.append(tmp.split()[0] + '-' + tmp.split()[1])
        else:
            break


CC_file = []
CC_prof = []
for i in range(len(putinrou)):
  CC_prof.append(pd.read_hdf(putinrou[i], 'cc_prof'))
  with h5py.File(putinrou[i], mode = "r") as file_rou:
      CC_file.append(np.array(file_rou[list(file_rou.keys())[0]]))
CC_array = CC_file[0]
CC_table = CC_prof[0]
f_prof = pd.read_hdf(putinrou[0], 'freq')
freq = np.array(f_prof[:])[:,0]
#==============selec Ccfs based on srcstat=====================================
part_CC = []
part_dis = []
sta_loc = []
CC_sta_name = []
cc_name = CC_table['0cc_name']
#print("CC_table: ", CC_table)
for i in range(len(cc_name)):
  cc_pairs_sta2 = cc_name[i].split('_')
  netsta1 = cc_pairs_sta2[0]
  netsta2 = cc_pairs_sta2[1]
  if (srcsta in cc_pairs_sta2) and (netsta1 in sele_netsta) and (netsta2 in sele_netsta):
     part_CC.append((CC_array[i,:]).copy())
     part_dis.append(CC_table.iloc[i,3]) #km
     if srcsta == netsta1:
         CC_sta_name.append(netsta2)
         sta_loc.append(staloc[syn_netsta.index(netsta2)])
     else:
         CC_sta_name.append(netsta1)
         sta_loc.append(staloc[syn_netsta.index(netsta1)])

sta_loc = np.array(sta_loc)
part_CC = np.array(part_CC)
part_dis = np.array(part_dis).copy()
#print("sta_loc: ", sta_loc)
#print("part_CC: ", len(part_CC))
#print("part_dis: ", len(part_dis))
dele_cc_rows = []
for CC_num in range(part_CC.shape[0]):
  if part_CC[CC_num,:].all() == 0:
      dele_cc_rows.append(CC_num)
  if part_dis[CC_num] == 0 or part_dis[CC_num]>= 800:
      dele_cc_rows.append(CC_num)
dele_cc_rows = list(set(dele_cc_rows))
dele_cc_rows.sort()
raw_rou_nonezerocc = np.delete(part_CC, dele_cc_rows, axis = 0)
dis_vin_nonezerocc = np.delete(part_dis, dele_cc_rows, axis = 0)
sta_loc_nonezerocc = np.delete(sta_loc, dele_cc_rows, axis = 0)
for counter, index in enumerate(dele_cc_rows):
    index = index-counter
    CC_sta_name.pop(index)

idel = getSNR(raw_rou_nonezerocc,dis_vin_nonezerocc,freq,SNRmin=0.8)
idel = list(set(idel))
idel.sort()
raw_rou_nonezerocc = np.delete(raw_rou_nonezerocc,idel,axis=0)
dis_vin_nonezerocc = np.delete(dis_vin_nonezerocc,idel,axis=0)
sta_loc_nonezerocc = np.delete(sta_loc_nonezerocc, idel, axis = 0)
for counter, index in enumerate(idel):
    index = index-counter
    CC_sta_name.pop(index)
#print("raw_rou_nonezerocc: ", raw_rou_nonezerocc.shape)
#print("dis_vin_nonezerocc: ", dis_vin_nonezerocc.shape)
#print("sta_loc_nonezerocc: ", sta_loc_nonezerocc.shape)
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

Srloc = staloc[syn_netsta.index(srcsta)]
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
Vz = pvar[:, 2, :]
Dz0 = integrate.cumtrapz(Vz, t, axis=1, initial=0)
stpindx = int(15/(t[1]-t[0]))
tc = t[:len(t)-stpindx]
dtsyn = tc[1] - tc[0]
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
indxcorr_syn = []
indxcorr_obs = []
r_plot = []
for i in range(len(robs)):
  if CC_sta_name[i] in syn_netsta:
      indxcorr_obs.append(i)
      indxcorr_syn.append(syn_netsta.index(CC_sta_name[i]))
      r_plot.append(robs[i])
r_plot = np.array(r_plot)
wave0sync = wave0syn[indxcorr_syn]
wave1sync = wave1syn[indxcorr_syn]
wave2sync = wave2syn[indxcorr_syn]
wave3sync = wave3syn[indxcorr_syn]

wave0obsc = wave0obs[indxcorr_obs]
wave1obsc = wave1obs[indxcorr_obs]
wave2obsc = wave2obs[indxcorr_obs]
wave3obsc = wave3obs[indxcorr_obs]
#==============================================================================
xlimt = 300 #s
lag_thre = 6 #s
ccc_thre = 0.6

for j in range(0,len(r_plot),stride):
  if r_plot[j] < bands[0,1]*aver[0]:
     continue
  xyn = ray_cross_path(srcloc[isrc-1], sta_loc_nonezerocc[j],0.5)
  #print("xyn 3-7s: ", xyn)
  t0 = r_plot[j]/aver[0]-5
  t1 = r_plot[j]/aver[0]+18
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave0obsc[j,:indx0] = 0
  wave0obsc[j,indx1:] = 0
  f = interpolate.interp1d(tc, wave0sync[j,:], kind='linear')
  tnew = tt[:int((xlimt)/dtobs)]
  wave0syn2 = f(tnew)
  indx0syn = int((t0-5)/dtobs)
  wave0syn2[:indx0syn] = 0
  wave0obs2 = wave0obsc[j,:int((xlimt)/dtobs)]
  if wave0obs2.any() == 0 or wave0syn2.any() == 0:
      continue
  wave0obs2 = wave0obs2/np.max(abs(wave0obs2))
  wave0syn2 = wave0syn2/np.max(abs(wave0syn2))
  xcorr = signal.correlate(wave0obs2,wave0syn2, mode='full')
  lags = signal.correlation_lags(len(wave0obs2), len(wave0syn2))
  if np.max(xcorr) == 0:
      continue
  xcorr /= np.max(xcorr)
  idmax = np.argmax(xcorr)
  lagtime = abs(lags[idmax])*dtobs/(r_plot[j]/100) #s/100km
  if(lagtime>lag_thre):
      continue
  indlag = lags[idmax]
  wave0syn_lag = np.zeros(len(wave0syn2))

  if (indlag <= 0):
      wave0syn_lag[0:len(wave0syn2)+indlag] = wave0syn2[0-indlag:]
  else:
      wave0syn_lag[indlag:] = wave0syn2[0:len(wave0syn2)-indlag]
  s1 = pd.Series(wave0obs2)
  s2 = pd.Series(wave0syn_lag)
  ccc = s1.corr(s2, method='pearson')
  if(ccc<ccc_thre) or (math.isnan(ccc)):
      continue
  fo = open(save_cell_3_7s, "a")
  for kk in range(len(xyn)):
      fo.write(str(int(xyn[kk,0])) + '\t')
      fo.write(str(int(xyn[kk,1])) + '\t')
      fo.write(str(lagtime) + '\t')
      fo.write(str(ccc) + '\n')
  fo.close()


for j in range(0,len(r_plot),stride):
  if r_plot[j] < bands[1,1]*aver[1]:
     continue
  xyn = ray_cross_path(srcloc[isrc-1], sta_loc_nonezerocc[j],0.5)
  #print("xyn 5-10s: ", xyn)
  t0 = r_plot[j]/aver[1]
  t1 = t0+23
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave1obsc[j,:indx0] = 0
  wave1obsc[j,indx1:] = 0
  f = interpolate.interp1d(tc, wave1sync[j,:], kind='linear')
  tnew = tt[:int((xlimt)/dtobs)]
  wave1syn2 = f(tnew)
  indx1syn = int((t0-5)/dtobs)
  wave1syn2[:indx1syn] = 0
  wave1obs2 = wave1obsc[j,:int((xlimt)/dtobs)]
  if wave1obs2.any() == 0 or wave1syn2.any() == 0:
      continue
  wave1obs2 = wave1obs2/np.max(abs(wave1obs2))
  wave1syn2 = wave1syn2/np.max(abs(wave1syn2))
  xcorr = signal.correlate(wave1obs2,wave1syn2, mode='full')
  lags = signal.correlation_lags(len(wave1obs2), len(wave1syn2))
  if np.max(xcorr) == 0:
      continue
  xcorr /= np.max(xcorr)
  idmax = np.argmax(xcorr)
  lagtime = abs(lags[idmax])*dtobs/(r_plot[j]/100) #s/100km
  if(lagtime>lag_thre):
      continue
  indlag = lags[idmax]
  wave1syn_lag = np.zeros(len(wave1syn2))

  if (indlag <= 0):
      wave1syn_lag[0:len(wave1syn2)+indlag] = wave1syn2[0-indlag:]
  else:
      wave1syn_lag[indlag:] = wave1syn2[0:len(wave1syn2)-indlag]


  s1 = pd.Series(wave1obs2)
  s2 = pd.Series(wave1syn_lag)
  ccc = s1.corr(s2, method='pearson')
  if(ccc<ccc_thre) or (math.isnan(ccc)):
      continue
  fo = open(save_cell_5_10s, "a")
  for kk in range(len(xyn)):
      fo.write(str(int(xyn[kk,0])) + '\t')
      fo.write(str(int(xyn[kk,1])) + '\t')
      fo.write(str(lagtime) + '\t')
      fo.write(str(ccc) + '\n')
  fo.close()

for j in range(0,len(r_plot),stride):
  if r_plot[j] < bands[2,1]*aver[2]:
     continue
  xyn = ray_cross_path(srcloc[isrc-1], sta_loc_nonezerocc[j],0.5)
  #print("xyn 8-15s: ", xyn)
  t0 = r_plot[j]/aver[2]
  t1 = t0 + 30
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave2obsc[j,:indx0] = 0
  wave2obsc[j,indx1:] = 0
  f = interpolate.interp1d(tc, wave2sync[j,:], kind='linear')
  tnew = tt[:int((xlimt)/dtobs)]
  wave2syn2 = f(tnew)
  indx2syn = int((t0-5)/dtobs)
  wave2syn2[:indx2syn] = 0
  wave2obs2 = wave2obsc[j,:int((xlimt)/dtobs)]
  if wave2obs2.any() == 0 or wave2syn2.any() == 0:
      continue
  wave2obs2 = wave2obs2/np.max(abs(wave2obs2))
  wave2syn2 = wave2syn2/np.max(abs(wave2syn2))
  xcorr = signal.correlate(wave2obs2,wave2syn2, mode='full')
  lags = signal.correlation_lags(len(wave2obs2), len(wave2syn2))
  if np.max(xcorr) == 0:
      continue
  xcorr /= np.max(xcorr)
  idmax = np.argmax(xcorr)
  lagtime = abs(lags[idmax])*dtobs/(r_plot[j]/100) #s/100km
  if(lagtime>lag_thre):
      continue
  indlag = lags[idmax]
  wave2syn_lag = np.zeros(len(wave2syn2))

  if (indlag <= 0):
      wave2syn_lag[0:len(wave2syn2)+indlag] = wave2syn2[0-indlag:]
  else:
      wave2syn_lag[indlag:] = wave2syn2[0:len(wave2syn2)-indlag]
  s1 = pd.Series(wave2obs2)
  s2 = pd.Series(wave2syn_lag)
  ccc = s2.corr(s1, method='pearson')
  if(ccc<ccc_thre) or (math.isnan(ccc)):
      continue
  fo = open(save_cell_8_15s, "a")
  for kk in range(len(xyn)):
      fo.write(str(int(xyn[kk,0])) + '\t')
      fo.write(str(int(xyn[kk,1])) + '\t')
      fo.write(str(lagtime) + '\t')
      fo.write(str(ccc) + '\n')
  fo.close()

for j in range(0,len(r_plot),stride):
  if r_plot[j] < bands[3,1]*aver[3]:
     continue
  xyn = ray_cross_path(srcloc[isrc-1], sta_loc_nonezerocc[j],0.5)
  #print("xyn 12-20s: ", xyn)
  t0 = r_plot[j]/aver[3]
  t1 = t0 + 45
  indx0 = int(t0/dtobs)
  indx1 = int(t1/dtobs)
  wave3obsc[j,:indx0] = 0
  wave3obsc[j,indx1:] = 0
  f = interpolate.interp1d(tc, wave3sync[j,:], kind='linear')
  tnew = tt[:int((xlimt)/dtobs)]
  wave3syn2 = f(tnew)
  indx3syn = int((t0-5)/dtobs)
  wave3syn2[:indx3syn] = 0
  wave3obs2 = wave3obsc[j,:int((xlimt)/dtobs)]
  if wave3obs2.any() == 0 or wave3syn2.any() == 0:
      continue
  wave3obs2 = wave3obs2/np.max(abs(wave3obs2))
  wave3syn2 = wave3syn2/np.max(abs(wave3syn2))
  xcorr = signal.correlate(wave3obs2,wave3syn2, mode='full')
  lags = signal.correlation_lags(len(wave3obs2), len(wave3syn2))
  if np.max(xcorr) == 0:
      continue
  xcorr /= np.max(xcorr)
  idmax = np.argmax(xcorr)
  lagtime = abs(lags[idmax])*dtobs/(r_plot[j]/100) #s/100km
  if(lagtime>lag_thre):
      continue
  indlag = lags[idmax]
  wave3syn_lag = np.zeros(len(wave3syn2))

  if (indlag <= 0):
      wave3syn_lag[0:len(wave3syn2)+indlag] = wave3syn2[0-indlag:]
  else:
      wave3syn_lag[indlag:] = wave3syn2[0:len(wave3syn2)-indlag]
  s1 = pd.Series(wave3obs2)
  s2 = pd.Series(wave3syn_lag)
  ccc = s2.corr(s1, method='pearson')
  if(ccc<ccc_thre) or (math.isnan(ccc)):
      continue
  fo = open(save_cell_12_20s, "a")
  for kk in range(len(xyn)):
      fo.write(str(int(xyn[kk,0])) + '\t')
      fo.write(str(int(xyn[kk,1])) + '\t')
      fo.write(str(lagtime) + '\t')
      fo.write(str(ccc) + '\n')
  fo.close()
