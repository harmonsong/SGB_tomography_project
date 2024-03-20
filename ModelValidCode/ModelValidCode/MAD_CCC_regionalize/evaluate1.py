#! /usr/bin/env python3
import numpy as np
import os, getopt, sys,glob
#------------------------------------------------------------
Forwardir = '/home/cjq/ModelValid/Forward/forLove/'
period_str = '3_7s'
dx = 0.5

try:
    opts, args = getopt.getopt(sys.argv[1:], 'p:d:', [
    'period=', 'dx='])
except getopt.GetoptError:
  sys.exit(2)

for opt, arg in opts:
  if opt in ('-p', '--period_str'):
    period_str = arg
  elif opt in ('-d', '--dx'):
    dx = arg

#------------------------------------------------------------
X0 = 108; Y0 = 20
dir_cells = Forwardir + '/dir_' + period_str
cell_names = glob.glob(dir_cells + '/*')
cell_names.sort()
for i in range(len(cell_names)):
#for i in range(1):
    cell_name = cell_names[i].split('/')[-1]
    idx = int(cell_name.split('_')[0])
    idy = int(cell_name.split('_')[-1].split('.')[0])
    lon = idx*dx + X0 - dx/2
    lat = idy*dx + Y0 - dx/2
    print('  cell_name:\t', cell_name)
    print('  idx:\t', idx)
    print('  idy:\t', idy)
    print('  lon:\t', lon)
    print('  lat:\t', lat)
    cell = np.loadtxt(cell_names[i])
    if (cell.shape[0] < 5):
        continue
    ccc_aver = np.average(cell[:,1]) # mean correlation coefficient: CCC

    #calculate mean absolute deviation of traveltime delay: CT_MAD
    lag_aver = np.average(cell[:,0])
    tmp = 0
    for j in range(cell.shape[0]):
        tmp += abs(cell[j,0] - lag_aver)
    ct_mad = tmp / cell.shape[0]
    with open(Forwardir + '/aver_' + period_str + '.txt', 'a+') as f:
        f.write(str(lon) + ' ' + str(lat) + ' ' + str(ct_mad) + ' ' + str(ccc_aver) + '\n')

    print('Done: ' + str(i+1) + '/' + str(len(cell_names)))


#print('  period_str:\t', period_str)
#print('  dir_cells:\t', dir_cells)

