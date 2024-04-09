#!/usr/bin/env python3
# -*- conding: utf-8 -*-
"""
Plot figure of 2-D seismogram for FD2Dker.tw Fortran program package.

Created on Thu Aug 17 10:17:25 2017
@author: Tche, USTC
"""

import sys, getopt, os.path
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import math

def print_help_info():
  print('\nUsage: fig4seismo2d.py [-h] | [-l | -p index] [-i inputpath]' +
    ' [-o outputpath] [-c fileconf] [-s scalefactor]')
  print('  [index      ]: seis-line/point number to plot for -l/p.')
  print('  [inputpath  ]: directory of the coordinate data files.')
  print('  [outputpath ]: directory of the seismogram data files.')
  print('  [fileconf   ]: name of the configure file.')
  print('  [scalefactor]: scale factor using to adjust the trace interval.\n')

#==============================================================================

index = 1
inputpath = 'input/'
outputpath = 'output/'
#inputpath = 'ss1_4_input_DP/'
#outputpath = 'ss1_4_output_DP/'
#inputpath = 'ss3_1_input/'
#outputpath = 'ss3_1_output/'
#inputpath = 'ss3_1_input_DP/'
#outputpath = 'ss3_1_output_DP/'
#inputpath = 'ss3_4_input/'
#outputpath = 'ss3_4_output/'
fileconf = 'SeisFD3D.conf'
scalefactor = 1.0
isline = False
LenFD = 3

#==============================================================================

try:
  opts, args = getopt.getopt(sys.argv[1:], 'hl:p:i:o:c:s:', ['help',
    'line=', 'point=', 'inputpath=', 'outputpath=', 'flieconf=',
    'scalefactor='])
except getopt.GetoptError:
  print_help_info()
  sys.exit(2)

for opt, arg in opts:
  if opt in ('-h', '--help'):
    print_help_info()
    sys.exit()
  elif opt in ('-l', '--line'):
    isline = True
    index = int(arg)
  elif opt in ('-p', '--point'):
    isline = False
    index = int(arg)
  elif opt in ('-i', '--inputpath'):
    inputpath = arg
  elif opt in ('-o', '--outputpath'):
    outputpath  = arg
  elif opt in ('-c', '--fileconf'):
    fileconf = arg
  elif opt in ('-s', '--scalefactor'):
    scalefactor = float(arg)

print('\nConfigure:')
if isline:
  print('  line:\t', index)
else:
  print('  point:\t', index)
  index = index - 1
print('  inputpath:\t', inputpath)
print('  outputpath:\t', outputpath)
print('  fileconf:\t', fileconf)
print('  scalefactor:\t', scalefactor, '\n')
print('\nDemos:\n./figseismo3d.py -l 3 --scalefactor=10\n')
#==============================================================================

conf = open(fileconf, 'r')
lines = conf.read().split('\n')

for line in lines:
  if line.find('#') >= 0: line = line[:line.find('#') - 1]
  if line.split(' ')[0] == 'dims':
    dims = [int(v) for v in line.split() if v.isdigit()][0:3]
  elif line.split(' ')[0] == 'ni':
    ni = int(line.split()[2])
  elif line.split(' ')[0] == 'nj':
    nj = int(line.split()[2])
  elif line.split(' ')[0] == 'nk':
    nk = int(line.split()[2])

conf.close()

if isline:#here decide the Line or Point
  iline = index #index is only specified as one number
else:
  iline = 0

#==============================================================================

whitems = ['Vx', 'Vy', 'Vz', 'Txx', 'Tyy', 'Tzz', 'Txy', 'Txz', 'Tyz']

isfirstget = True
isgottime = False
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


#==============================================================================


pname = ['$V_x$', '$V_y$', '$V_z$', r'$\tau_{xx}$', r'$\tau_{yy}$', r'$\tau_{zz}$', 
  r'$\tau_{xy}$', r'$\tau_{xz}$', r'$\tau_{yz}$']
xbgl = [t[0], t[- 1]]
ybgl = np.array([- 0.5, - 0.5])

if isline:
  for iw in range(9):
    traceshift = abs(pvar[:, iw, :]).max()*2/scalefactor #offset when plot
    fig = plt.figure(figsize = (8, 6), dpi = 80)
    for ip in range(npt):
      plt.plot(t, pvar[ip, iw, :] + traceshift*ip, 'b')
#     plt.plot(xbgl, ybgl*traceshift + traceshift*ip, ':k', linewidth = 0.5) #for horinzital line
#   plt.plot(xbgl, ybgl*traceshift + traceshift*npt, ':k', linewidth = 0.5)
    plt.yticks(traceshift*(np.arange(npt + 1.0) - 0.5), np.arange(npt + 1))
    plt.title('Seismogram of No.{} seismo-line for {}'.format(index,
      pname[iw]))
    plt.xlabel('time (s)')
    if iw < 3:
      plt.ylabel('velocity ({:.2f} km/s)'.format(traceshift))
    else:
      plt.ylabel('stress ({:.2f} GPa)'.format(traceshift))
else:
  # ==== VELOCITY ====
  vmin = np.min(pvar[index,0:3,:])
  vmax = np.max(pvar[index,0:3,:])
  #vmax = np.ceil(vmax)
  #vmin = np.floor(vmin)
  fig = plt.figure(figsize = (8, 6), dpi = 80)
  fig.subplots_adjust(hspace=0.5)# adjust distance between small plots
  # -- Vx  --
  plt.subplot(311).plot(t, pvar[index, 0, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[0])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('velocity (km/s)')
  #plt.ylim([vmin,vmax])
  # -- Vy  --
  plt.subplot(312).plot(t, pvar[index, 1, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[1])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('velocity (km/s)')
  #plt.ylim([vmin,vmax])
  # -- Vz  --
  plt.subplot(313).plot(t, pvar[index, 2, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[2])
  plt.xlabel('time (s)')
  plt.ylabel('velocity (km/s)')
  #plt.ylim([vmin,vmax])
  
  # ==== STRESS ====
  vmin = np.min(pvar[index,3:9,:])
  vmax = np.max(pvar[index,3:9,:])
  #vmax = np.ceil(vmax)
  #vmin = np.floor(vmin)
  fig = plt.figure(figsize = (8, 6), dpi = 80)
  fig.subplots_adjust(hspace=0.7) # 10. for with xlabel
  # -- Txx --
  plt.subplot(611).plot(t, pvar[index, 3, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[3])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
  # -- Tyy --
  plt.subplot(612).plot(t, pvar[index, 4, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[4])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
  # -- Tzz --
  plt.subplot(613).plot(t, pvar[index, 5, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[5])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
  # -- Txy --
  plt.subplot(614).plot(t, pvar[index, 6, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[6])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
  # -- Txz --
  plt.subplot(615).plot(t, pvar[index, 7, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[7])
  plt.gca().get_xaxis().set_visible(False)
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
  # -- Tyz --
  plt.subplot(616).plot(t, pvar[index, 8, :], 'b')
  plt.title('Seismogram of No.{} seismo-point for '.format(index+1) + pname[8])
  plt.xlabel('time (s)')
  plt.ylabel('stress (GPa)')
  #plt.ylim([vmin,vmax])
plt.show()
plt.savefig('test.png')


