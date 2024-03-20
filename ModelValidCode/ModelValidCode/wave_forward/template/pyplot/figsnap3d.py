#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot snapshot of 2-D wavefield for FD2Dker.tw Fortran program package.

Created on Fri Ang 18 20:24:56 2017
@author: Tche, USTC
"""

import sys, time, getopt, os.path
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import matplotlib.animation as anm
import glob as gb
import multiprocessing as mp

def print_help_info():
  print('\nUsage: fig4snap2d.py [-h] | [-i inputpath] [-o outputpath]' +
    ' [-c fileconf] [-s index] [-w whichitem] [-t tsamp] [-n nprofile] [-d direction')
  print('  [inputpath ]: directory of the coordinate files.')
  print('  [outputpath]: directory of the wavefield data files.')
  print('  [fileconf  ]: name of the configure file.')
  print('  [whichitem ]: name of the wavefield component (Vx, Vy, Vz, Txx,' +
    ' Tzz, Txy, Txz or Tyz).')
  print('  [tsamp     ]: sampling paramter (t_start/t_end/t_stride)' + 
    ' of time axis.')
  print('  [index     ]: specify the snap number.')
  print('  [direction ]: specify the showing direction(xz,yz,xy)')
  print('  [nprofile  ]: specify the profile number(absolutely location).\n')
  print('  Demos:\n ./figsnap3d.py -t 1/-1/20 -w Vy --index=2 --nprofile=1 -d xy --Crestrict=600\n')

#==============================================================================

inputpath = 'input/'
outputpath = 'output/'
fileconf = 'SeisFD3D.conf'
index = 1
whichitem = 'Vx'
tsamp = [1, -1, 1]
LenFD = 3
nprofile = 1
direction = 'xz'
Cres = 0.0
outpath = 'picture/'

#==============================================================================

try:
  opts, args = getopt.getopt(sys.argv[1:], 'hi:o:c:s:n:d:w:t:cs:', ['help',
    'inputpath=', 'outputpath=', 'fileconf=', 'index=', 'nprofile=', 'direction=', 'whichitem=',
    'tsamp=','Crestrict='])
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
    outputpath = arg
  elif opt in ('-c', '--fileconf'):
    fileconf = arg
  elif opt in ('-s', '--index'):
    index = int(arg)
  elif opt in ('-n', '--nprofile'):
    nprofile = int(arg)
  elif opt in ('-cs', '--Crestrict'):
    Cres = float(arg)
  elif opt in ('-d', '--direction'):
    direction = arg
  elif opt in ('-w', '--whichitem'):
    whichitem = arg
  elif opt in ('-t', '--tsamp'):
    tsamp = [int(v) for v in arg.split('/')]

print('\nConfigure:')
print('  inputpath:\t', inputpath)
print('  outputpath:\t', outputpath)
print('  fileconf:\t', fileconf)
print('  whichitem:\t', whichitem)
print('  Snap index:\t', index)
print('  Profile index: ',nprofile)
print('  Showing plane: ',direction)
print('  tsamp:\t', tsamp, '\n')
if Cres!=0 :
  print(' Apply colorbar restriction as :', Cres, '\n')

#==============================================================================

sub=np.zeros(12,dtype=int)

conf = open(fileconf, 'r')
lines = conf.read().split('\n')
snapname = 'snap_{:03d}'.format(index)

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
  elif line.split(' ')[0] == snapname:
    sub[0:9] = [int(v) for v in line.split() if v.replace('-', '', 1).isdigit()][0:9]#6 par starts and counts

conf.close()

if 'sub' not in globals():
  raise IOError('{}: No such snapshot number.'.format(index))

#sub['start','count','stride','end']
#start from 1, means index is 0-199 when total = 200
#count = how many point that have already contains start and end point
#stride= step
#ending point location = 200, which index = 199
# start + (count-1)*stride = end-1
# 1+(37-1)*2= 74-1 

# X: 0 3 6 9
# Y: 1 4 7 10
# Z: 2 5 8 11


# -- x --
sub[0] = max(sub[0], 1)
if sub[3] == -1:# '-1' means 'total'
  sub[3] = ni*dims[0] #represents counts
  sub[9] = ni*dims[0] #represents 'ending location (index should -1)
else:
  sub[3] = min(ni*dims[0], sub[3])
  sub[9] = min(ni*dims[0], sub[0]+sub[3])
# -- y --
sub[1] = max(sub[1], 1)
if sub[4] == -1:# '-1' means 'ends'
  sub[4] = nj*dims[1]
  sub[10]= nj*dims[1]
else:
  sub[4] = min(nj*dims[1], sub[4])
  sub[10]= min(nj*dims[1], sub[1]+sub[4])
# -- z --
sub[2] = max(sub[2], 1)
if sub[5] == -1:
  sub[5] = nk*dims[2]
  sub[11]= nk*dims[2]
else:
  sub[5] = min(nk*dims[2], sub[5])
  sub[11]= min(nk*dims[2], sub[2]+sub[5])

#==============================================================================

if whichitem in ['Vx', 'Vy', 'Vz']:
  filepref = 'vel_{:03d}_'.format(index)
  pname = r'$V_' + whichitem[1] + '$'
elif whichitem in ['Txx', 'Tyy', 'Tzz', 'Txy', 'Txz', 'Tyz']:
  filepref = 'sgt_{:03d}_'.format(index)
  pname = r'$\tau_{' + whichitem[1:3] + '}$'
else:
  raise IOError(whichitem + ': Invalid wavefield component.')

#reading time variables
for k in range(dims[2]):
	for j in range(dims[1]):
		for i in range(dims[0]):
			filename = '{}{}mpi{:02d}{:02d}{:02d}_n00001.nc'.format(outputpath, filepref, i, j, k)
			if not os.path.isfile(filename):
				continue
			snp = nc.Dataset(filename, 'r')
			if len(snp.variables['time'])==0 :
				print('Ignore the file ',filename,' since the data is empty')
				continue
			t = snp.variables['time'][:]

if tsamp[1] == -1: tsamp[1] = t.size # t read from snp
tnumber= int( 1 + (tsamp[1]-tsamp[0])/tsamp[2] )
tseries=t[tsamp[0]-1:tsamp[1]:tsamp[2]]

if direction=='xz':
	if nprofile<sub[1]:
		nprofile=sub[1]
	if nprofile> int(sub[1]+(sub[4]-1)*sub[7]):
		nprofile=int(sub[1]+(sub[4]-1)*sub[7])
	nprofile=int((nprofile-sub[1])/sub[7])
	print('the Profile was Modified as',nprofile)
elif direction=='yz':
	if nprofile<sub[0]:
		nprofile=sub[0]
	if nprofile>int(sub[0]+(sub[3]-1)*sub[6]):
		nprofile=int(sub[0]+(sub[3]-1)*sub[6])
	nprofile=int((nprofile-sub[0])/sub[6])
	print('the Profile was Modified as',nprofile)
elif direction=='xy':
	if nprofile<sub[2]:
		nprofile=sub[2]
	if nprofile>int(sub[2]+(sub[5]-1)*sub[8]):
		nprofile=int(sub[2]+(sub[5]-1)*sub[8])
	nprofile=int((nprofile-sub[2])/sub[8])
	print('the Profile was Modified as',nprofile)

if sub[3]==1:
	if direction != 'yz':
		direction='yz'
	print('the ploting plane change to ',direction,' due to data shape')
elif sub[4]==1:
	if direction != 'xz':
		direction='xz'
	print('the ploting plane change to ',direction,' due to data shape')
elif sub[5]==1:
	if direction != 'xy':
		direction='xy'
	print('the ploting plane change to ',direction,' due to data shape')




dvar=np.zeros((tnumber,sub[5],sub[4],sub[3]))
zvar=np.zeros((sub[5],sub[4],sub[3]))
yvar=np.zeros((sub[5],sub[4],sub[3]))
xvar=np.zeros((sub[5],sub[4],sub[3]))
gsubs=np.empty(3,dtype=int)
gsube=np.empty_like(gsubs)
subs=np.empty_like(gsubs)
sube=np.empty_like(gsubs)
subt=np.empty_like(gsubs)
Sin=np.empty_like(gsubs)
Ein=np.empty_like(gsubs)

#reading data and coord
for k in range(dims[2]):
	for j in range(dims[1]):
		for i in range(dims[0]):
			filename = '{}{}mpi{:02d}{:02d}{:02d}_n00001.nc'.format(outputpath, filepref, i, j, k)
			if not os.path.isfile(filename):
				continue
			snp = nc.Dataset(filename, 'r')
			if len(snp.variables['time'])==0 :
				print('Ignore the file ',filename,' since the data is empty')
				continue
			
			gsubs=snp.getncattr('gsubs')
			gsube=snp.getncattr('gsube')
			subs=snp.getncattr('subs')
			sube=snp.getncattr('sube')
			subt=snp.getncattr('subt')

			Sin[0]=int((gsubs[0]-sub[0])/sub[6])
			Sin[1]=int((gsubs[1]-sub[1])/sub[7])
			Sin[2]=int((gsubs[2]-sub[2])/sub[8])
			Ein[0]=int((gsube[0]-sub[0])/sub[6])
			Ein[1]=int((gsube[1]-sub[1])/sub[7])
			Ein[2]=int((gsube[2]-sub[2])/sub[8])
			
			dvar[:,Sin[2]:Ein[2]+1,Sin[1]:Ein[1]+1,Sin[0]:Ein[0]+1]=snp.variables[whichitem][tsamp[0]-1:tsamp[1]:tsamp[2],:,:,:]
			
			snp.close()
	                
			filename = '{}coord_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, i, j, k)
			if not os.path.isfile(filename):
				continue
			crd = nc.Dataset(filename, 'r')
			
			subs=subs+LenFD
			sube=sube+LenFD

			xvar[Sin[2]:Ein[2]+1,Sin[1]:Ein[1]+1,Sin[0]:Ein[0]+1]=\
			crd.variables['x'][subs[2]-1:sube[2]:subt[2],subs[1]-1:sube[1]:subt[1],subs[0]-1:sube[0]:subt[0]]
			yvar[Sin[2]:Ein[2]+1,Sin[1]:Ein[1]+1,Sin[0]:Ein[0]+1]=\
			crd.variables['y'][subs[2]-1:sube[2]:subt[2],subs[1]-1:sube[1]:subt[1],subs[0]-1:sube[0]:subt[0]]
			zvar[Sin[2]:Ein[2]+1,Sin[1]:Ein[1]+1,Sin[0]:Ein[0]+1]=\
			crd.variables['z'][subs[2]-1:sube[2]:subt[2],subs[1]-1:sube[1]:subt[1],subs[0]-1:sube[0]:subt[0]]
			
			crd.close()

xplt=xvar/1.0e3
yplt=yvar/1.0e3
zplt=zvar/1.0e3

#confirm the ploting axes
if direction=='yz':
	print('ploting Y-Z imaging')
	coord=['y','z']
	#FLdim=[dims[1],dims[2],1]
	last=np.squeeze(zplt[:,:,nprofile])
	fast=np.squeeze(yplt[:,:,nprofile])
	dvar=np.squeeze(dvar[:,:,:,nprofile])
elif direction=='xz':
	print('ploting X-Z imaging')
	coord=['x','z']
	#FLdim=[dims[0],dims[2],2]
	last=np.squeeze(zplt[:,nprofile,:])
	fast=np.squeeze(xplt[:,nprofile,:])
	dvar=np.squeeze(dvar[:,:,nprofile,:])
elif direction=='xy':
	print('ploting X-Y imaging')
	coord=['x','y']
	#FLdim=[dims[0],dims[1],3]
	last=np.squeeze(yplt[nprofile,:,:])
	fast=np.squeeze(xplt[nprofile,:,:])
	dvar=np.squeeze(dvar[:,nprofile,:,:])

#==============================================================================




if 1:
  # --- Matlab-like dynamically plot ---
  fig = plt.figure(figsize = (7, 6), dpi = 80)
  for it in range(tnumber):
    plt.clf()
    plt.pcolormesh(fast, last, dvar[it, :, :], cmap = 'seismic', rasterized = True)
    print('time=\n',tseries[it])
    plt.gca().set_aspect('equal')
    plt.xlabel(coord[0]+' distance (km)')
    plt.ylabel(coord[1]+' distance (km)')
    plt.title('{} snapshot of No.{:d} when t = {:.2f} s'.format(pname, index,tseries[it]))
    plt.colorbar()
    if Cres!=0:
      plt.clim([-1.0*Cres,Cres])
    plt.draw()
    plt.pause(0.0001)
  # time.sleep(.2)
else:
  # --- Matplotlib animation plot ---
  imgs = []
  fig = plt.figure(figsize = (7, 6), dpi = 80)
  for it in range(tsamp[0], tsamp[1], tsamp[2]):
    img = plt.pcolormesh(fast, last, dvar[it, :, :], cmap = 'jet', rasterized = True,
      animated = True)
    print('timeIMG=\n',it)
    plt.axes().set_aspect('equal')
    plt.xlabel(coord[0]+' distance (km)')
    plt.ylabel(coord[1]+' distance (km)')
    plt.title('{} snapshot of No.{:d} when t = {:.2f} s'.format(pname, index,
      t[it]))
    #plt.colorbar()
    imgs.append([img])
  # plt.clf()
  ani = anm.ArtistAnimation(fig, imgs, interval = 50, blit = True)#, repeat_delay = 0)

plt.show()



'''



#plot source and receiver
def replot(plt):
	plt.plot(recx,recy,'ko',marksize=3,fillstyle='none')
	plt.plot(srcx,srcy,'w*',marksize=3,fillstyle='none')

#==============================================================================
ts = 0
tt = 1
te = dvar.shape[0]
clims = (-1000,1000)
sc = 1e0

def mpplot(it):
	plt.pcolormesh(fast,last,dvar[it,:,:], cmap='seismic',rasterized=True)
	#rsplot(plt)
	plt.gca().set_aspect('equal')
	plt.xlabel(coord[0]+' distance (km)')
	plt.ylabel(coord[1]+' distance (km)')
	plt.title('{} snapshot of No.{:d} at t = {:.2f} s'.format(pname, index,tseries[it]))
	plt.colorbar()
	if Cres!=0:
		plt.clim([-1.0*Cres,Cres])
	if 0:
		if 0:
			plt.clim(clims)
		else:
			climit=np.max(np.abs(dvar[it,:,:]))/sc
			plt.clim([-climit,climit])
	plt.savefig(outpath+'snap{:03d}.png'.format(it),dpi=300,format='png')
	plt.clf()

#=============================================================================
#parallel plot

if 1:
	for ip in gb.glob(outpath+'*.png'):
		os.remove(ip)
	print('remove all existing pictures and start to save new')

	pool = mp.Pool(20)
	re = pool.map(mpplot, range(ts,te,tt))
	print('picture {}snap{:03d}.png was generated\n'.format(outpath,ts))

'''







