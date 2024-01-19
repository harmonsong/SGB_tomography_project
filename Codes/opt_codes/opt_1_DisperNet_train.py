import numpy as np
import h5py
import os 
import sys
sys.path.append('/home/harmon/data/Codes/FJ-codes/DisperNet_online')
import dispernet


proj_name = 'output_FJSJ_05-01/'
tag = 5

info_basic = np.load(proj_name+'Basic_info.npy',allow_pickle=True).item()
dir_ds = info_basic['dir_ds']
key_ds = info_basic['key_subworks']

# create trainset files
dir_train = 'trainset/ ' + str(tag) + '/'
if not os.path.exists(dir_train):
    os.makedirs(dir_train)
h5file = dir_train + 'h5/'
curvefile = dir_train + 'curve/'
if not os.path.exists(h5file):
    os.makedirs(h5file)
if not os.path.exists(curvefile):
    os.makedirs(curvefile)


# Transform to used h5 files
f_flag = 20
for key in key_ds:    
    data = h5py.File(dir_ds+'/ds_'+key+'.h5', 'r')
    f0 = data['f'][:]
    f = f0[f0<f_flag]
    c = data['c'][:]
    amp = data['ds_linear'][:][0]
    amp = amp[:,f0<f_flag]
    outname = h5file+'ds_'+str(tag)+'_'+key +'.h5'
    data.close()
    dispernet.save2h5(amp, f, c,fileName=outname)

# Run Dispernet
dispernet.App(filePath=h5file, curveFilePath=curvefile,freqSeries=f, trigerMode=False, searchStep=2, cmap='jet', periodCutRate=0.125, semiAutoRange=0.1, autoT=True, url='http://10.20.64.63:8514')

# transfer training
#dispernet.createTrainSet('./trainSetDAS.h5', inputfile, outputfile)