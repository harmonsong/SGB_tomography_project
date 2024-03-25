# %%
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
import shutil

# cd /shdisk/rem2/Harmon/F-J/San/project_repartition_v4.0/output_repar_v9.2_01-01/inv_dispernet/
# mpirun -np 60 /home/songsh/git_repo/F-J/tools_F-J/toollib_DisbaCode/inversion.py --num_init 10  --config config_inv.yml --key select_50_60

# %% [markdown]
# ### Parameters

# %% [markdown]
# Vs:  Mordret 2019
# Vp/Vs: Zigone 2019 (3-4)     (0-100m)
# Density: Zigone 2019 (2-2.5) (0-100m)
# 
# 
# Vp: Share 2020 (0.6-2.6) (0-70m)

# %% [markdown]
flag_project = 1 # 0--regular; 1--repartition
# %%
if flag_project == 0:
    file_project = 'a-project.yml'
elif flag_project == 1:
    file_project = 'a-project_repar.yml'
elif flag_project == 2:
    file_project = 'a-project_voro.yml'
    
with open(file_project, 'r', encoding='utf-8') as f:
    proj = yaml.load(f.read(), Loader=yaml.FullLoader)
name_project = proj['name']

#name_project = 'project/output_FJSJ_16-01/'               
#name_project = 'project_repartrition/repartrition_01-03/'               
#name_project = 'project_voronoi/voronoi_01-03/'         

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

# %%
fmin = 1
fmax = 30
nf = 291
info_basic['BFGS_for_fmin'] = fmin
info_basic['BFGS_for_fmax'] = fmax
info_basic['BFGS_for_nf'] = nf


# %%
old_curve_path = 'Structure/disp_data_5-15/'
key_olds = []
if os.path.exists(old_curve_path):
    # read all files in the folder
    files = os.listdir(old_curve_path)
    for file in files:
        key_olds.append(file[3:8])

# %%
rdir_inv = 'inversion_BFGS/'
dir_inv = dir_project + rdir_inv
info_basic['rdir_inv_BFGS'] = rdir_inv
dir_file = dir_inv + 'initial/'
dir_image = dir_project + info_basic['dir_image'] + 'initial_model/'
if os.path.exists(dir_image) == False:
    os.makedirs(dir_image)
if os.path.exists(dir_inv) == False:
    os.makedirs(dir_inv)
if os.path.exists(dir_file) == False:
    os.makedirs(dir_file)

#%%
#将config_inv.yml从San_Jansinto复制到
outname_config = dir_inv+'config_inv.yml'
outname_config_fund = dir_inv+'config_inv_fund.yml'
if os.path.exists(outname_config)==False:
    shutil.copyfile('config_inv.yml', dir_inv+'config_inv.yml')
if os.path.exists(outname_config_fund)==False:
    shutil.copyfile('config_inv_fund.yml', dir_inv+'config_inv_fund.yml')

# %%
V0 = 0.3
alpha = 0.25
#beta = 3
#A1 = 0.4
#A2 = 0.1
d0 = 0.0

# %%
def Vs_model(V0,alpha,d0,d):
    Vs = V0 * ((d*1e3+1)**alpha - (d0*1e3+1)**alpha + 1)
    return Vs
        

# %%
def Vs_model_linear(V0,k,d):
    return V0 + k*d
        

# %%
def brocher(z, vs):
    vp = np.zeros(len(z))
    rho = np.zeros(len(z))
    for i in range(len(vs)):
        vp[i] = (0.9409 + 2.0947 * vs[i] - 0.8206 * vs[i]**2 +
              0.2683 * vs[i]**3 - 0.0251 * vs[i]**4)
        rho[i] = (1.6612 * vp[i] - 0.4721 * vp[i]**2 + 0.0671 * vp[i]**3
               - 0.0043 * vp[i]**4 + 0.000106 * vp[i]**5)
    return vp,rho


def gardner(z, vs):
    vp = np.zeros(len(z))
    rho = np.zeros(len(z))

    vp_vs_ratio = 1.7321
    vp = vp_vs_ratio * vs
    rho = 1.741 * vp ** 0.25

    return vp,rho


def user_defined(z, vs):
    vp = np.zeros(len(z))
    rho = np.zeros(len(z))
    vp_vs_ratio = 1.5
    vp = vp_vs_ratio * vs
    for i in range(len(vs)):
        vp[i] = (0.9409 + 2.0947 * vs[i] - 0.8206 * vs[i]**2 +
              0.2683 * vs[i]**3 - 0.0251 * vs[i]**4)
        rho[i] = (1.6612 * vp[i] - 0.4721 * vp[i]**2 + 0.0671 * vp[i]**3
               - 0.0043 * vp[i]**4 + 0.000106 * vp[i]**5)
        
    return vp,rho

# %%
#N = 50
#dz = 0.006
#N = 200
#dz = 0.0015  
N = 100
dz = 0.003
flag_relation = 1
info_basic['inv_BFGS_layers'] = N
info_basic['inv_BFGS_dz'] = dz
info_basic['inv_BFGS_empirical'] = flag_relation

# %%
layers = np.linspace(1,N,N)
depths = np.zeros(N)
for i in range(N):
    depths[i] = d0 + i*dz
Vp = np.zeros(N)
Vs = np.zeros(N)
rho = np.zeros(N)
for i in range(N):
    d = depths[i]
    Vs[i] = Vs_model(V0,alpha,d0,d)
if flag_relation == 1:
    Vp,rho = brocher(depths,Vs)
elif flag_relation == 2:
    Vp,rho = gardner(depths,Vs)
elif flag_relation == 3:
    Vp,rho = user_defined(depths,Vs)

# %%
def plot_initial(Vp,Vs,rho,depths,tag):
    plt.figure(figsize=(8,6))
    plt.plot(Vp,depths,'r',label='Vp')
    plt.plot(Vs,depths,'b',label='Vs')
    plt.plot(rho,depths,'g',label='rho')
    plt.legend()
    # 翻转y轴
    plt.gca().invert_yaxis()
    plt.xlim(0,4)
    plt.savefig(dir_image+'initial_model_'+str(tag)+'.png')
    #plt.close()

# %%
# write a function to write the initial model as txt file, which seperate by space，with four significant digits
def write_initial_model(dir_file,layers,depths,Vp,Vs,rho,tag):
    with open(dir_file+'/initial_model_'+str(tag)+'.txt','w') as f:
        for i in range(len(layers)):
            f.write('{} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(int(layers[i]),depths[i],Vp[i],Vs[i],rho[i]))

# %%
tag = 5
flag_plot = 0
write_initial_model(dir_file,layers,depths,Vp,Vs,rho,tag)
if flag_plot == 1:
    plot_initial(Vp,Vs,rho,depths,tag)

with open(dir_project+'Basic_info.yml', 'w', encoding='utf-8') as f:
   yaml.dump(data=info_basic, stream=f, allow_unicode=True)