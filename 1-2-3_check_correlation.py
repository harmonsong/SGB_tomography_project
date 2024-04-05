import numpy as np
import os 
import sys, getopt
import yaml
import h5py

with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)
name_CC = 'CC/CC_150_prewhiten/'
dir_CC = os.path.join(dir_CC_workspace,name_CC)
print('dir_CC: ', dir_CC)

d_start = 128
d_len = 31
flag = 0 # 0--ncfs; 1--stack 
try:
    opts, args = getopt.getopt(sys.argv[1:], "s:l:f:", ["start=", "len=", "flag="])
except getopt.GetoptError:
    print('1-2-2_correlation_day.py -s <day start> -l <day len>')
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-s", "--start"):
        d_start = int(arg)
    elif opt in ("-l", "--len"):
        d_len = int(arg)
    elif opt in ("-f", "--flag"):
        flag = int(arg)


if flag == 0:
    print('check correlation: ' + str(d_start) + ' ' + str(d_start + d_len-1))
    d_end = d_start + d_len
    y = 2014
    for d in range(d_start,d_end):
        year = str(y)
        day = "%03d"%d
        outname = os.path.join(dir_CC,year+'-'+day+'.npz')
        try:
            data = np.load(outname)
            ncfs = data['ncfs']
            nan = np.sum(np.isnan(ncfs))
            shape = ncfs.shape
            print('day', d, 'shape: ', shape, ' nan: ', nan)
        except:
            print('no file: ', outname)
            continue
else:
    print('check stack: ')
    try:
        ncffile = h5py.File(dir_CC + 'gather_all.h5','r')
        ncfs = ncffile['ncfs'][:]
        ncffile.close()
        shape = ncfs.shape
        nan = np.sum(np.isnan(ncfs))
        print('shape: ', shape, ' nan: ', nan)
    except:
        print('error in reading: ', dir_CC + 'gather_all.h5')
        pass
