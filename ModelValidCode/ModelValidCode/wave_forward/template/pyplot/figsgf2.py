#!/usr/bin/env python3
# -*- conding: utf-8 -*-
"""
Synthesize the noise cross-correlation function by convoluatating the source time function and the Green's function synthesized by 3-D FD wave simulation code.

Created on Fri Jun 24 10:48 2022
@author: Juqing Chen, USTC
"""

import sys, getopt, os.path
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import math
from scipy import integrate
from scipy.fftpack import fft


def print_help_info():
  print('\nUsage: fig4seismo2d.py [-h] | [-i inputpath]' +
    ' [-o outputpath] [-S fileSrconf] [-c fileFDconf] [-d stride] [-s scalefactor]')
  print('  [inputpath  ]: directory of the coordinate data files.')
  print('  [outputpath ]: directory of the seismogram data files.')
  print('  [fileSrconf   ]: name of the configure file for source.')
  print('  [fileFDconf   ]: name of the configure file for FD.')
  print('  [stride   ]: number of stations to skip for plotting.')
  print('  [scalefactor]: scale factor using to adjust the trace interval.\n')

#==============================================================================

inputpath = 'input/'
outputpath = 'output/'
fileSrconf = 'SeisSource.conf'
fileFDconf = 'SeisFD3D.conf'
scalefactor = 1.0
stride = 1
isline = False
LenFD = 3

#==============================================================================

try:
  opts, args = getopt.getopt(sys.argv[1:], 'hl:i:o:S:c:d:s:', ['help',
    'inputpath=', 'outputpath=', 'fileSrconf=', 'fileFDconf=',
    'stride=', 'scalefactor='])
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
  elif opt in ('-S', '--fileSrconf'):
    fileSrconf = arg
  elif opt in ('-c', '--fileFDconf'):
    fileFDconf = arg
  elif opt in ('-d', '--stride'):
    stride = int(arg)
  elif opt in ('-s', '--scalefactor'):
    scalefactor = float(arg)

print('\nConfigure:')
print('  inputpath:\t', inputpath)
print('  outputpath:\t', outputpath)
print('  fileSrconf:\t', fileSrconf)
print('  fileFDconf:\t', fileFDconf)
print('  stride:\t', stride)
print('  scalefactor:\t', scalefactor, '\n')
#print('\nDemos:\n./figseismo3d.py -l 3 --scalefactor=10\n')
#==============================================================================
confSr = open(fileSrconf, 'r')
linesSr = confSr.read().split('\n')
for i in range(len(linesSr)):
  if linesSr[i].split(' ')[0] == '<anchor_force>':
    Srloc = [float(v) for v in linesSr[i+1].split()][0:2]
    break
confSr.close()
#==============================================================================
confFD = open(fileFDconf, 'r')
linesFD = confFD.read().split('\n')
staloc = []
for j in range(len(linesFD)):
  line = linesFD[j]
  if line.find('#') >= 0: line = line[:line.find('#') - 1]
  if line.split(' ')[0] == 'dims':
    dims = [int(v) for v in line.split() if v.isdigit()][0:3]
  elif line.split(' ')[0] == 'number_of_recv':
    nsta = int(line.split()[2])
  elif line.split('_')[0] == 'recv':
    staloc.append([float(v) for v in line.split()[2:4]])

confFD.close()
Srloc = np.array(Srloc)
staloc = np.array(staloc)
#print("Srloc: ", Srloc.shape)
#print("staloc: ", staloc.shape)
#==============================================================================
r0 = np.zeros(nsta)
for k in range(nsta):
  r0[k] = np.sqrt((staloc[k,0]-Srloc[0])**2+ (staloc[k,1]-Srloc[1])**2)
indx = np.argsort(r0)
r = r0[indx]
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

#print(np.shape(dvar))
#np.set_printoptions(threshold=np.nan)
#print(dvar[1,0,:]) 
#for i in range(1,nt):
#	print('i=',i-1,'\tvalue=',pvar[1,0,i])


#means:
#'1'--->the second point in point-series(start from 0)
#'0'--->the 'Vx'/first item in whichitem
#':'--->total time series

print(np.shape(pvar))
#plt.plot(t,pvar[3,0,:])
#plt.show()

# Integrate to obtain Displacement
Vx = pvar[indx, 0, :]
Vy = pvar[indx, 1, :]
Vz = pvar[indx, 2, :]
Dx = integrate.cumtrapz(Vx, t, axis=1, initial=0)
Dy = integrate.cumtrapz(Vy, t, axis=1, initial=0)
Dz = integrate.cumtrapz(Vz, t, axis=1, initial=0)

# FFT spectrum analysis
#fdx = abs(fft(Dx))
#fdy = abs(fft(Dy))
#fdz = abs(fft(Dz))
#Fs = 1.0/(t[1]-t[0])
#df = Fs / len(t)
#f = df*np.arange(len(t))
#==============================================================================


pname = ['$V_x$', '$V_y$', '$V_z$','$D_x$','$D_y$','$D_z$'] 
xbgl = [t[0], t[- 1]]
ybgl = np.array([- 0.5, - 0.5])
# ==== VELOCITY ====
# -- Vx  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Vx[i]/np.max(abs(Vx[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[0])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')
# -- Vy  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Vy[i]/np.max(abs(Vy[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[1])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')
# -- Vz  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Vz[i]/np.max(abs(Vz[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[2])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')
# -- Dx  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Dx[i]/np.max(abs(Dx[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[3])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')
# -- Dy  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Dy[i]/np.max(abs(Dy[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[4])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')
# -- Dz  --
fig = plt.figure(figsize = (8, 6), dpi = 80)
for i in range(0,nsta,stride):
  plt.plot(t, Dz[i]/np.max(abs(Dz[i]))*100+r[i]/1e3, 'b')
plt.title('Seismogram for ' + pname[5])
plt.ylabel('distance (km)')
plt.xlabel('time (s)')

plt.show()



