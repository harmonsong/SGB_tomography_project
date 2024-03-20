#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot figure of 2-D grid for FD2Dker.tw Fortran program package.

Created on Wed Aug 16 19:53:46 2017
@author: Tche, USTC
"""

import sys, getopt, os.path
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

def print_help_info():
  print('\nUsage: figgrid3d.py [-h] | [-i inputpath] [-c fileconf]' +
    ' [-d direction] [-n nprofile] [-x xystep] [-z zstep]')
  print('  [inputpath]: directory of the coordinate data files.')
  print('  [fileconf ]: name of the configure file.')
  print('  [direction]: display the specified direction(x or y).')
  print('  [nprofile ]: number of profiles to be ploted.')
  print('  [xystep    ]: step length to skip on x or y-direction.')
  print('  [zstep    ]: step length to skip on z-direction.\n')
  print('  Demos:\n ./figgrid3d.py --nprofile=100 -x 5 -z 20 -d y \n')

#==============================================================================

inputpath = 'input/'
fileconf = 'SeisFD3D.conf'
direction = 'x'
nprofile = 1
xystep = 1
zstep = 1
LenFD = 3

#==============================================================================

try:
  opts, args = getopt.getopt(sys.argv[1:], 'hi:c:d:n:x:z:', ['help', 'inputpath=',
    'fileconf=', 'direction=', 'nprofile=',  'xystep=', 'zstep='])
except getopt.GetoptError:
  print_help_info()
  sys.exit(2)

for opt, arg in opts:
  if opt in ('-h', '--help'):
    print_help_info()
    sys.exit()
  elif opt in ('-i', '--inputpath'):
    inputpath = arg
  elif opt in ('-c', '--fileconf'):
    fileconf = arg
  elif opt in ('-d', '--direction'):
    direction = arg
  elif opt in ('-n', '--nprofile'):
    nprofile = int(arg)
  elif opt in ('-x', '--xystep'):
    xystep = int(arg)
  elif opt in ('-z', '--zstep'):
    zstep = int(arg)

print('\nConfigure:')
print('  inputpath:\t', inputpath)
print('  fileconf:\t', fileconf)
print('  direction:\t', direction)
print('  nprofile:\t', nprofile)
print('  xystep:\t', xystep)
print('  zstep:\t', zstep, '\n')

#==============================================================================

conf = open(fileconf, 'r')
lines = conf.read().split('\n')

for line in lines:
  if line.find('#') >= 0: 
	  line = line[:line.find('#') - 1]
  if line.split(' ')[0] == 'dims':
    dims = [int(v) for v in line.split() if v.isdigit()][0:3]
  elif line.split(' ')[0] == 'ni':
    ni = int(line.split()[2])
  elif line.split(' ')[0] == 'nj':
    nj = int(line.split()[2])
  elif line.split(' ')[0] == 'nk':
    nk = int(line.split()[2])

conf.close()

#print(dims,ni,nj,nk)

#==============================================================================
if direction == 'x':
  displaydir = 'y'
  xyvar = np.zeros((dims[2]*nk, dims[0]*ni))
  zvar = np.zeros((dims[2]*nk, dims[0]*ni))
  jmod = int((nprofile-1.0)/nj) #profile start from x to nDirection
  j = nprofile - jmod*nj +LenFD
  ijend = ni*dims[0] - 1
  for k in range(dims[2]):
    for i in range(dims[0]):
      filename = '{}coord_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, i, jmod, k)
      if not os.path.isfile(filename):
        raise IOError(filename + ': No such file.')
      crd = nc.Dataset(filename, 'r')
      xyvar[k*nk:(k + 1)*nk, i*ni:(i + 1)*ni] = crd.variables['x'][LenFD:LenFD +
        nk, j, LenFD:LenFD + ni]
      zvar[k*nk:(k + 1)*nk, i*ni:(i + 1)*ni] = crd.variables['z'][LenFD:LenFD +
        nk, j, LenFD:LenFD + ni]
      crd.close()
else:
  displaydir = 'x'
  zvar = np.zeros((dims[2]*nk, dims[1]*nj))
  xyvar = np.zeros((dims[2]*nk, dims[1]*nj))
  imod = int((nprofile-1.0)/ni) #profile start from x to nDirection
  i = nprofile - imod*ni +LenFD #used for specify the profile
  ijend = nj*dims[1] - 1
  for k in range(dims[2]):
    for j in range(dims[1]):
      filename = '{}coord_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, imod, j, k)
      if not os.path.isfile(filename):
        raise IOError(filename + ': No such file.')
      crd = nc.Dataset(filename, 'r')
      #print(crd)
      #print(crd.variables.keys())
      xyvar[k*nk:(k + 1)*nk, j*nj:(j + 1)*nj] = crd.variables['y'][LenFD:LenFD +
        nk, LenFD:LenFD + nj, i]
      zvar[k*nk:(k + 1)*nk, j*nj:(j + 1)*nj] = crd.variables['z'][LenFD:LenFD +
        nk, LenFD:LenFD + nj, i]
      crd.close()
xyvar = xyvar/1.0e3
zvar = zvar/1.0e3

#print(xyvar.shape,zvar.shape)
#==============================================================================

kend = nk*dims[2] - 1
fig = plt.figure(figsize = (7, 7), dpi = 80)

print('Zsize=',np.size(zvar),'Xsize=',np.size(xyvar))
print('Zshape=',np.shape(zvar),'Xshape=',np.shape(xyvar))

for i in range(0, ijend, xystep):
  plt.plot(xyvar[:, i], zvar[:, i], 'k')
plt.plot(xyvar[:, ijend], zvar[:, ijend], 'k')
for k in range(kend, 0, - zstep):
  plt.plot(xyvar[k, :], zvar[k, :], 'k')
plt.plot(xyvar[0, :], zvar[0, :], 'k')
# plt.axes([xyvar.min(), xyvar.max(), zvar.min(), zvar.max()])
plt.axes().set_aspect('equal')
plt.xlabel(displaydir+' distance (km)' )
plt.ylabel('z depth (km)')
plt.title('Physical meshgrid profile No.'+str(nprofile))
plt.show()
# plt.saveas('test.eps')
