import numpy as np
import h5py
import os 
import sys
import shutil
import yaml
sys.path.append('../tools_F-J/toollib_DisperNet_local/')
sys.path.append('../tools_F-J/')
from toollib_standard import plotlib
import dispernet_local_latest as dispernet
import pandas as pd


#%%
flag_project = 1 # 0--regular ; 1--repartition; 2--voronoi
flag_repick = 1 # 0-- no new h5 trans; 1-- yes
flag_forward = 1 # 0-- no forward; 1-- yes
num_refs = 1 # sort with nearst num_refs centroid
flag_partrition = 1 # plot partition
flag_plot_or = 0    # plot non-remove FJ
r_near = 1        # add nearest num_near dispersion picks

#%%
if flag_project == 0:
    file_project = 'a-project.yml'
elif flag_project == 1:
    file_project = 'a-project_repar.yml'
elif flag_project == 2:
    file_project = 'a-project_voro.yml'
    
with open(file_project, 'r', encoding='utf-8') as f:
    proj = yaml.load(f.read(), Loader=yaml.FullLoader)
name_project = proj['name']

#%%
with open('0_config.yml', 'r', encoding='utf-8') as f:
    dir_config = yaml.load(f.read(), Loader=yaml.FullLoader)
dir_project_workspace = dir_config['dir_project_workspace']
dir_CC_workspace = dir_config['dir_CC_workspace']
print('dir_CC_workspace: ', dir_CC_workspace)
print('dir_project_workspace: ', dir_project_workspace)
dir_project = os.path.join(dir_project_workspace, name_project)
print('dir_project: ', dir_project)

#%%
filename = dir_project+'Basic_info.yml'
with open(filename, 'r', encoding='utf-8') as f:
    info_basic = yaml.load(f.read(), Loader=yaml.FullLoader)
filename_bi = dir_project+'Basic_info.npy'
info_basic_bi = np.load(filename_bi, allow_pickle='TRUE').item()      # setting dictionary

#%%
stainfo = pd.read_excel(info_basic['stalistname_all'])
nsta_all = len(stainfo.iloc[:,0])
stalist_all = stainfo['Station'].tolist()
lat_stations_all = stainfo['latitude'].tolist() 
lon_stations_all = stainfo['longitude'].tolist()

# Affine transformation
def compute_affine_transform(original_points, target_points):
    A_matrix = np.array([[original_points[0][0], original_points[0][1], 1, 0, 0, 0],
                         [0, 0, 0, original_points[0][0], original_points[0][1], 1],
                         [original_points[1][0], original_points[1][1], 1, 0, 0, 0],
                         [0, 0, 0, original_points[1][0], original_points[1][1], 1],
                         [original_points[2][0], original_points[2][1], 1, 0, 0, 0],
                         [0, 0, 0, original_points[2][0], original_points[2][1], 1]])

    A1_B1_C1 = np.array([target_points[0][0], target_points[0][1], target_points[1][0], target_points[1][1], target_points[2][0], target_points[2][1]])

    coefficients = np.linalg.solve(A_matrix, A1_B1_C1)

    affine_matrix = np.array([[coefficients[0], coefficients[1], coefficients[2]],
                               [coefficients[3], coefficients[4], coefficients[5]],
                               [0, 0, 1]])

    return affine_matrix
lon_stations_all_new = []
lat_stations_all_new = []
for sta in stalist_all:
    if int(sta[1:3]) <= 60:
        lon_stations_all_new.append(lon_stations_all[stalist_all.index(sta)])
        lat_stations_all_new.append(lat_stations_all[stalist_all.index(sta)])
refs = ['R0101','R6001','R6020']
lon_refs = [lon_stations_all[stalist_all.index(ref)] for ref in refs]
lat_refs = [lat_stations_all[stalist_all.index(ref)] for ref in refs]
loc_refs = np.column_stack([lon_refs,lat_refs])
loc_refs_new = np.array([[0,0],[600,0],[600,600]])
affine_matrix = compute_affine_transform(loc_refs, loc_refs_new)

faults = np.load('clark_faults.npy', allow_pickle='TRUE').item()

#%%
dir_ds = dir_project + info_basic['rdir_ds']
dir_partition = dir_project + info_basic['rdir_partition']
#key_ds  = info_basic['key_subworks_repick']


key_ds = []
"""
#filename = dir_project+info_basic['rdir_inv_BFGS']+'repick_1.txt'
filename = dir_project+'repick_1.txt'
nums = np.loadtxt(filename,dtype='int')  
nums = [277]
#nums = range(0,600)
nums = [str(x) for x in nums]
for key_subwork in info_basic['key_subworks']:
    #print(key_subwork.split('--')[0])
    if key_subwork.split('--')[0] in nums:
        key_ds.append(key_subwork)
"""
key_ds = info_basic['key_subworks']
#print(nums)
#key_ds = info_basic['key_subworks']


"""

key_ds = []
for key in info_basic['key_subworks']:
    num = int(key[3:])
    if num == 14:
        key_ds.append(key)
"""
"""
key_ds = []
for key in info_basic['key_subworks']:
    num = int(key[:2])
    if num >= 32 and num <= 42:
        key_ds.append(key)
"""

#%% 
info_basic['rdir_disp_pick'] = 'disp_pick/'
dir_disp_pick = dir_project + info_basic['rdir_disp_pick']
if not os.path.exists(dir_disp_pick):
    os.makedirs(dir_disp_pick)
with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)

#%%
#r_max = info_basic_bi['r_max']
r_max = {}
for key in key_ds:
    r_max[key] = 10

#%%
fmax = 30
fmin = 1
cmax = 2
cmin = 0.25
inputfile = dir_disp_pick + 'h5/'
outputfile = dir_disp_pick
f = info_basic_bi['f']
c = np.linspace(info_basic['fj_c_min'],info_basic['fj_c_max'],info_basic['fj_c_num'])
if flag_repick == 1:
    if os.path.exists(inputfile):
        os.system('rm -rf '+inputfile)
    os.makedirs(inputfile)
    if not os.path.exists(outputfile):
        os.makedirs(outputfile)
    for key in key_ds:
        data = h5py.File(dir_ds+'/ds_'+str(key)+'.h5', 'r')
        #f = data['f'][:]    
        #c = data['c'][:]
        #amp = data['ds_remove'][:][0][:,f<fmax]
        #amp_or = data['ds_linear'][:][0][:,f<fmax]
        ds_remove = data['ds_remove'][:][0]
        ds_linear = data['ds_linear'][:][0]
        """
        amp = plotlib.smooth_ds(ds_remove)
        amp = plotlib.smooth_ds(ds_remove,level_c  =10)
        amp_or = plotlib.smooth_ds(ds_linear)
        """
        amp = ds_remove
        amp_or = ds_linear
        amp = amp[:,np.logical_and(f>fmin,f<fmax)]
        amp = amp[np.logical_and(c>cmin,c<cmax),:]
        amp_or = amp_or[:,np.logical_and(f>fmin,f<fmax)]
        amp_or = amp_or[np.logical_and(c>cmin,c<cmax),:]
        outname = inputfile+'ds_'+str(key) +'.h5'
        data.close()
        dispernet.save2h5(amp, f[np.logical_and(f>fmin,f<fmax)], c[np.logical_and(c>cmin,c<cmax)],fileName=outname,spectrum_or = amp_or)

#%%
#old_curve_path = 'oldData/disp_data_15-01/'
old_curve_path = 'oldData/disp_data_repar_01-02/'
#old_curve_path = ''
key_olds = []
if os.path.exists(old_curve_path):
    print('exits old curve')
    # read all files in the folder
    files = os.listdir(old_curve_path)
    for file in files:
        if file.find('.xlsx') == -1:
            key_olds.append(file[file.find('_')+1:file.find('curve')])

#%%
key_fund = []
fund_curve_path = dir_project + 'inversion_BGFS/disp_model_fund/nono/'
over_curve_path = dir_project + 'inversion_BGFS/disp_model/nono/'
if flag_forward == 1:
    fund_curve_path = dir_project + info_basic['rdir_inv_BFGS']+'disp_model_fund/'
    if os.path.exists(fund_curve_path):
        print('exits fund curve')
        # read all files in the folder
        files = os.listdir(fund_curve_path)
        for file in files:
            key_fund.append(file[file.find('_')+1:file.find('curve')])
    over_curve_path = dir_project + info_basic['rdir_inv_BFGS']+'disp_model/'

#%%    
"""
fileList = os.listdir(inputfile)
num_fileList = [int(file[file.find('_')+1:file.find('--')]) for file in fileList]
index = np.argsort(num_fileList)
fileList = np.array(fileList)[index].tolist()
"""
# sort fileList based on location
fileList_or = os.listdir(inputfile)
key_ds = [file[file.find('_')+1:file.find('.h5')] for file in fileList_or]
loc_ds = {}
for key in key_ds:
    loc_ds[key] = [0,0]
    filePath = dir_project + info_basic['rdir_partition'] + str(key) +'.txt'
    station,lat,lon = np.loadtxt(filePath,dtype = 'str',unpack=True)
    lat = lat.tolist()
    lon = lon.tolist()
    lat_centroid = np.mean(np.array(lat).astype(float))
    lon_centroid = np.mean(np.array(lon).astype(float))
    loc_ds[key] = [lat_centroid,lon_centroid]

fileList = []
fileList.append(fileList_or[0])
fileList_or.pop(0)
fileList_all = os.listdir(inputfile)

key_subworks_all = info_basic['key_subworks']
lat_centroid_partition_all = []
lon_centroid_partition_all = []
for key in key_subworks_all:
    filepath = dir_partition + str(key) + '.txt'
    stations_this, lat_stations_this, lon_stations_this = np.loadtxt(filepath, dtype='str' , unpack=True)
    lat_centroid_partition_all.append(np.mean(lat_stations_this.astype(float)))
    lon_centroid_partition_all.append(np.mean(lon_stations_this.astype(float)))
lat_centroid_partition_all = np.array(lat_centroid_partition_all)
lon_centroid_partition_all = np.array(lon_centroid_partition_all)

while fileList_or != []:
    key_ref = fileList[-1][fileList[-1].find('_')+1:fileList[-1].find('.h5')]
    loc_ref = loc_ds[key_ref]

    dist = (lat_centroid_partition_all-loc_ref[0])**2+(lon_centroid_partition_all-loc_ref[1])**2
    
    #找到dist中最小的num_refs个值
    key_refs = []
    loc_refs = []
    dist_near = []

    indexs = np.argsort(dist)
    for i in range(num_refs):
        key_refs.append(key_subworks_all[indexs[i]])
        loc_refs.append([lat_centroid_partition_all[indexs[i]],lon_centroid_partition_all[indexs[i]]])
    #print(key_refs)
    
    for file in fileList_or:
        key = file[file.find('_')+1:file.find('.h5')]
        loc = loc_ds[key]
        dist_this = 0
        for i in range(num_refs):
            dist_this += np.sqrt((loc[0]-loc_refs[i][0])**2+(loc[1]-loc_refs[i][1])**2)
        dist_near.append(dist_this)


    index = np.argmin(dist_near)
    fileList.append(fileList_or[index])
    fileList_or.pop(index)

#%%
# Run Dispernet
#dispernet.App(filePath=inputfile, curveFilePath=outputfile,freqSeries=f, trigerMode=False, searchStep=2, cmap='jet', periodCutRate=0.12, semiAutoRange=0.1, autoT=True, url='http://10.20.64.63:8514')
print('fmax = '+str(fmax))
v_min = 0.01
#dispernet.App(r_flag = r_max,vmin = v_min,oldfile=old_curve_path,oldkeys= key_olds,fundfile = fund_curve_path,overfile = over_curve_path,fundkeys = key_fund,filePath=inputfile, curveFilePath=outputfile,freqSeries=f[f<fmax], trigerMode=False, searchStep=2, cmap='jet', periodCutRate=0.8, semiAutoRange=0.1, autoT=True, url='http://10.20.64.63:8514')
dispernet.App(info_basic,lon_stations_all,lat_stations_all,fileList,r_near = r_near,matrix = affine_matrix,faults = faults,file_project = file_project,flag_plot_or=flag_plot_or,flag_plot_partrition=flag_partrition,vmin = v_min,oldfile=old_curve_path,oldkeys= key_olds,fundfile = fund_curve_path,overfile = over_curve_path,fundkeys = key_fund,filePath=inputfile, curveFilePath=outputfile,freqSeries=f[f<fmax], trigerMode=False, searchStep=2, cmap='jet', periodCutRate=0.2, semiAutoRange=0.03, autoT=True, url='http://10.20.64.63:8514')

# transfer training
#dispernet.createTrainSet('./trainSetDAS.h5', inputfile, outputfile)