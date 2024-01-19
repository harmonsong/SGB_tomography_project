import numpy as np
import h5py
import os 
import sys
sys.path.append('/home/harmon/data/Codes/FJ-codes/DisperNet_online')
import dispernet

tag = 8

# create trainset files
dir_train = 'trainset/'+str(tag)+'/'
if not os.path.exists(dir_train):
    os.makedirs(dir_train)
h5file = dir_train + 'h5/'
curvefile = dir_train + 'curve/'


# Run Dispernet
#dispernet.App(filePath=h5file, curveFilePath=curvefile,freqSeries=f, trigerMode=False, searchStep=2, cmap='jet', periodCutRate=0.125, semiAutoRange=0.1, autoT=True, url='http://10.20.64.63:8514')

# transfer training
dispernet.createTrainSet(dir_train+'/trainSetDAS_'+str(tag)+'.h5', h5file, curvefile)