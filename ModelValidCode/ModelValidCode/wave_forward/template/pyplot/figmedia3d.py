#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot figure of 2-D media for FD2Dker.tw Fortran program package.

Created on Wed Aug 16 15:20:23 2017
@author: Tche, USTC
"""

import sys, getopt, os.path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
import numpy as np
import netCDF4 as nc

def print_help_info():
  print('\nUsage: fig4media2d.py [-h] | [-m mapdim] [-d direction] [-i inputpath] [-c fileconf]' +
    ' [-n nprofile] [-w whichitem] [-s scalefactor] [-g grid] [-x xystep] [-z zstep]')
  print('  [mapdim   ]: plot in 2D or 3D.')
  print('  [direction]: extracting from the specified direction(x or y).')
  print('  [inputpath]: directory of the coordinate and media data files.')
  print('  [fileconf ]: name of the configure file.')
  print('  [nprofile ]: number of the specified profile, 2D for single number 3D for mutil number.')
  print('  [whichitem]: name of the elastic parameter (rho, Vp, Vs,lambda or mu) to plot.\n')
  print('  [grid     ]: whether to append the meshgrid.')
  print('  [scalefactor]: coordinate space ratio when plot. (3D only)')
  print('  [xystep   ]: step length to skip on x or y-direction. (2D only)')
  print('  [zstep    ]: step length to skip on z-direction.  (2D only)\n')

#==============================================================================

mapdim = 2
direction = 'x'
inputpath = 'input/'
fileconf = 'SeisFD3D.conf'
nprofile = 1
whichitem = 'rho'
LenFD = 3
grid = 0
scalefactor = 1
xystep = 1
zstep = 1

#==============================================================================

try:
  opts, args = getopt.getopt(sys.argv[1:], 'hm:d:i:c:n:w:s:g:x:z:', ['help', 'mapdim=', 'direction=', 'inputpath=',
    'fileconf=', 'nprofile=', 'whichitem=', 'scalefactor=', 'grid=', 'xystep=', 'zstep='])
except getopt.GetoptError:
  print_help_info()
  sys.exit(2)

for opt, arg in opts:
  if opt in ('-h', '--help'):
    print_help_info()
    sys.exit()
  elif opt in ('-m', '--mapdim'):
    mapdim = int(arg)
  elif opt in ('-d', '--direction'):
    direction = arg
  elif opt in ('-i', '--inputpath'):
    inputpath = arg
  elif opt in ('-c', '--fileconf'):
    fileconf = arg
  elif opt in ('-n', '--nprofile'):
    nprofile = arg
  elif opt in ('-w', '--whichitem'):
    whichitem = arg
  elif opt in ('-s', '--scalefactor'):
    scalefactor = int(arg)
  elif opt in ('-g', '--grid'):
    grid = int(arg)
  elif opt in ('-x', '--xystep'):
    xystep = int(arg)
  else:
    zstep = int(arg)
  #elif opt in ('-z', '--zstep'):
  #  zstep = int(arg)

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

#==============================================================================
np_list=nprofile.split(',',-1)
lengthprofile=len(np_list)
if lengthprofile==1:
	nprof=int(nprofile) #string to int
	if direction=='x':
		if nprof>nj*dims[1] or nprof<1 :
			print('exceeds the Y boundary')
	else:
		if nprof>ni*dims[0] or nprof<1 :
			print('exceeds the X boundary')
else:
	np_list_int=list(map(int,np_list))
	if mapdim==2:
		lengthprofile==1
		nprof=np_list_int[0]
	else:
		nprof=np_list_int
	if direction=='x':
		if nprof>nj*dims[1] or nprof<1 :
			print('exceeds the Y boundary')
	else:
		if nprof>ni*dims[0] or nprof<1 :
			print('exceeds the X boundary')

print('\nConfigure:')
print('  mapdim:\t', mapdim)
print('  direction:\t', direction)
print('  inputpath:\t', inputpath)
print('  fileconf:\t', fileconf)
print('  nprofile:\t', nprof)
print('  whichitem:\t', whichitem)
print('  grid:\t', grid)
print('  dims=',dims,' ni=',ni,' nj=',nj,' nk=',nk)
if mapdim==3:
  	print('  scalefactor:\t', scalefactor,'\n')
else:
	print('  xystep:\t', xystep)
	print('  zstep:\t', zstep, '\n')

print('  Now only valid for 2D case (single slice)')
print('  that means, nprofile should be single and')
print('  scalefactor is also useless!!!')
print('  Demos as: ')
print('  ./figmedia3d.py -m 2 -d x --xystep=5 -w Vs -g 1 -z 20 --nprofile="100"')


if direction=='y':
	displaydir = 'x'
	x = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	y = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	z = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	rho = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	lamda = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	mu = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	pvar = np.zeros((dims[2]*nk, lengthprofile, dims[0]*ni))
	
	
	jidx=int((nprof-1.0)/nj)
	for jth in range(lengthprofile):
		if lengthprofile!=1:
			jmod=jidx[jth]
			js=nprof[jth]-jmod*nj+LenFD
		else:
			jmod=jidx
			js=nprof-jmod*nj+LenFD
		for i in range(dims[0]):
			for k in range(dims[2]):
				filename = '{}media_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, i, jmod, k)
				if not os.path.isfile(filename):
					raise IOError(filename + ': No such file.')
				med = nc.Dataset(filename, 'r')
				rho[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=med.variables['rho'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				lamda[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=med.variables['lambda'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				mu[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=med.variables['mu'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				med.close()
				
				filename = '{}coord_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, i, jmod, k)
				if not os.path.isfile(filename):
					raise IOError(filename + ': No such file.')
				crd = nc.Dataset(filename, 'r')
				x[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=crd.variables['x'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				y[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=crd.variables['y'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				z[k*nk:(k+1)*nk,jth,i*ni:(i+1)*ni]=crd.variables['z'][LenFD:LenFD+nk,js,LenFD:LenFD+ni]
				crd.close()
else:
	displaydir = 'y'
	x = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	y = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	z = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	rho = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	lamda = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	mu = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
	pvar = np.zeros((dims[2]*nk, dims[1]*nj, lengthprofile))
 

	iidx=int((nprof-1.0)/ni)
	for ith in range(lengthprofile):
		if lengthprofile!=1:
			imod=iidx[ith]
			istart=nprof[ith]-imod*ni+LenFD
		else:
			imod=iidx
			istart=nprof-imod*ni+LenFD
		for j in range(dims[1]):
			for k in range(dims[2]):
				filename = '{}media_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, imod, j, k)
				if not os.path.isfile(filename):
					raise IOError(filename + ': No such file.')
				med = nc.Dataset(filename, 'r')
				rho[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=med.variables['rho'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				lamda[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=med.variables['lambda'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				mu[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=med.variables['mu'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				med.close()
				
				filename = '{}coord_mpi{:02d}{:02d}{:02d}.nc'.format(inputpath, imod, j, k)
				if not os.path.isfile(filename):
					raise IOError(filename + ': No such file.')
				crd = nc.Dataset(filename, 'r')
				x[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=crd.variables['x'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				y[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=crd.variables['y'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				z[k*nk:(k+1)*nk,j*nj:(j+1)*nj,ith]=crd.variables['z'][LenFD:LenFD+nk,LenFD:LenFD+nj,istart]
				crd.close()
x=x/1.0e3
y=y/1.0e3
z=z/1.0e3

#==============================================================================
tiny=1e-12
rho=rho+tiny

if whichitem == 'Vp':
  pvar=np.sqrt(np.divide(lamda+2*mu,rho))/1.0e3
  pname = {'title':'Vp', 'unit':'km/s'}
elif whichitem == 'Vs':
  pvar=np.sqrt(np.divide(mu,rho))/1.0e3
  pname = {'title':'Vs', 'unit':'km/s'}
elif whichitem == 'rho':
  pvar=rho/1.0e3
  pname = {'title':r'$\rho$', 'unit':'$kg/m^3$'}
elif whichitem == 'mu':
  pvar=mu/1.0e3
  pname = {'title':r'$\mu$', 'unit':'GPa'}
elif whichitem == 'lamda':
  pvar=lamda/1.0e9
  pname = {'title':r'$\lambda', 'unit':'GPa'}
else:
  raise IOError(whichitem + ': No such elastic parameter.')

fig = plt.figure(figsize = (7, 6), dpi = 80)

#ax = fig.gca(projection='3d')

if mapdim==2:
	if direction=='y':
		xy=np.squeeze(x)
		ijend=ni*dims[0]-1
	else:
		xy=np.squeeze(y)
		ijend=nj*dims[1]-1
	z=np.squeeze(z)
	pvar=np.squeeze(pvar)
	
	plt.pcolormesh(xy,z,pvar,cmap='jet',rasterized=True)
	plt.axes().set_aspect('equal')
	plt.xlabel(displaydir+' distance (km)')
	plt.ylabel('z depth (km)')
	plt.title('Media map for {} ({})'.format(pname['title'], pname['unit']))
	plt.colorbar()
	if grid:
		for i in range(0,ijend,xystep):
			plt.plot(xy[:,i],z[:,i],'k')
		plt.plot(xy[:,ijend],z[:,ijend],'k')
		for k in range(nk*dims[2]-1,0,-zstep):
			plt.plot(xy[k,:],z[k,:],'k')
		plt.plot(xy[0,:],z[0,:],'k')
else:
	if direction=='y':
		for jth in range(lengthprofile):
			x=np.squeeze(x[:,jth,:])
			z=np.squeeze(z[:,jth,:])
			pvar=np.squeeze(pvar[:,jth,:])
			ijend=ni*dims[0]-1
			plt.pcolormesh(x,z,pvar,cmap='jet',rasterized=True)
		if grid:
			for i in range(0,ijend,xystep):
				plt.plot(x[:,i],z[:,i],'k')
			plt.plot(x[:,ijend],z[:,ijend],'k')
			for k in range(nk*dims[2]-1,0,-zstep):
				plt.plot(x[k,:],z[k,:],'k')
			plt.plot(x[0,:],z[0,:],'k')
	else:
		for ith in range(lengthprofile):
			y=np.squeeze(y[:,:,ith])
			z=np.squeeze(z[:,:,ith])
			pvar=np.squeeze(pvar[:,:,ith])
			ijend=nj*dims[1]-1
			plt.pcolormesh(y,z,pvar,cmap='jet',rasterized=True)
		if grid:
			for j in range(0,ijend,xystep):
				plt.plot(y[:,j],z[:,j],'k')
			plt.plot(y[:,ijend],z[:,ijend],'k')
			for k in range(nk*dims[2]-1,0,-zstep):
				plt.plot(y[k,:],z[k,:],'k')
			plt.plot(y[0,:],z[0,:],'k')
			
	plt.axes().set_aspect('equal')
	plt.xlabel(displaydir+' distance (km)')
	plt.ylabel('z depth (km)')
	plt.title('Media map for {} ({})'.format(pname['title'], pname['unit']))
	plt.colorbar()
plt.show()
# plt.savefig('test.eps')
'''	
#else for the 3D surface case, but failed becaused of ax.plot_surface (shape problem) and set_aspect(scalefactor problem)
else:
	if direction=='x':
		print('before:Zshape=',np.shape(z),'Xshape=',np.shape(x),'Yshape=',np.shape(y),'pvarshape=',np.shape(pvar))
		surf=ax.plot_surface(np.squeeze(x[:,0,:]),np.squeeze(y[:,:,0]),np.squeeze(z[0,:,:]),np.squeeze(pvar[0,:,:]),cmap='jet');
		print('Zsize=',np.size(z),'Xsize=',np.size(x),'Ysize=',np.size(y),'pvarsize=',np.size(pvar))
		print('after:Zshape=',np.shape(z),'Xshape=',np.shape(x),'Yshape=',np.shape(y),'pvarshape=',np.shape(pvar))
		for jth in range(1,lengthprofile,1):
			surf=ax.plot_surface(np.squeeze(x[:,jth,:]),np.squeeze(y[:,:,jth,]),np.squeeze(z[jth,:,:]),np.squeeze(pvar[jth,:,:]),cmap='jet');
		if grid:
			for jth in range(lengthprofile):
				surf=ax.plot3D(np.squeeze(x[:,jth,0]),np.squeeze(y[:,jth,0]),\
				np.squeeze(z[:,jth,0]),'k')
				surf=ax.plot3D(np.squeeze(x[:,jth,ni-1]),np.squeeze(y[:,jth,ni-1]),\
				np.squeeze(z[:,jth,ni-1]),'k')
				surf=ax.plot3D(np.squeeze(x[0,jth,:]),np.squeeze(y[0,jth,:]),\
				np.squeeze(z[0,jth,:]),'k')
				surf=ax.plot3D(np.squeeze(x[nk-1,jth,:]),np.squeeze(y[nk-1,jth,:]),\
				np.squeeze(z[nk-1,jth,:]),'k')
		#plt.axes().set_aspect([scalefactor,1,scalefactor])
	else:
		surf=ax.plot_surface(np.squeeze(x[0,:,:]),np.squeeze(y[0,:,:]),\
		np.squeeze(z[0,:,:]),np.squeeze(pvar[0,:,:]),cmap='jet');
		for ith in range(1,lengthprofile,1):
			surf=ax.plot_surface(np.squeeze(x[ith,:,:]),np.squeeze(y[ith,:,:]),\
			np.squeeze(z[ith,:,:]),np.squeeze(pvar[ith,:,:]),cmap='jet');
		if grid:
			for ith in range(lengthprofile):
				surf=ax.plot3D(np.squeeze(x[:,0,ith]),np.squeeze(y[:,0,ith]),\
				np.squeeze(z[:,0,ith]),'k')
				surf=ax.plot3D(np.squeeze(x[:,nj-1,ith]),np.squeeze(y[:,nj-1,ith]),\
				np.squeeze(z[:,nj-1,ith]),'k')
				surf=ax.plot3D(np.squeeze(x[0,:,ith]),np.squeeze(y[0,:,ith]),\
				np.squeeze(z[0,:,ith]),'k')
				surf=ax.plot3D(np.squeeze(x[nk-1,:,ith]),np.squeeze(y[nk-1,:,ith]),\
				np.squeeze(z[nk-1,:,ith]),'k')
		#plt.axes().set_aspect([scalefactor,scalefactor,1])
	ax.set_xlabel(direction+' distance (km)')
	ax.set_zlabel('z depth (km)')
	ax.set_title('Media map for {} ({})'.format(pname['title'], pname['unit']))
	#fig.colorbar(surf)
plt.show()

'''
