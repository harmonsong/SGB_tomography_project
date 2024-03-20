#!/usr/bin/env python
from __future__ import division
import numpy as np
from shapely.geometry.polygon import Polygon
from shapely.geometry import MultiPoint, Point
from matplotlib.colors import ListedColormap
import pygmt, concurrent.futures
from multiprocessing import Manager
import os,sys
from scipy.io import loadmat
#=======================================================================================================
allsta = np.loadtxt('/home/cjq/ModelValid/Forward/SCB_locat_v2.dat').tolist()
P = MultiPoint(allsta).convex_hull
periods = ['3_7s', '5_10s', '8_15s', '12_20s']
labels_mad = ['(a1)', '(b1)', '(c1)','(d1)']
labels_ccc = ['(a2)', '(b2)', '(c2)','(d2)']
#periods = ['17_25s', '23_30s', '28_36s']
#labels_mad = ['(a1)', '(b1)', '(c1)']
#labels_ccc = ['(a2)', '(b2)', '(c2)']
#--------------------mesh grid--------------------------------------------------------------
lon_lat_size = 0.05
grid_lon = np.linspace(108,123,num=int((123-108)/lon_lat_size+1))
grid_lat = np.linspace(20,36,num=int((36-20)/lon_lat_size+1))
loncun, latcun = np.meshgrid(grid_lon, grid_lat) #lat*lon
lonflatten = loncun.flatten().reshape(loncun.flatten().shape[0],1)
latflatten = latcun.flatten().reshape(latcun.flatten().shape[0],1)

#lonlatMadCcc finally obtain eight conclumns, that is, lon, lat, lag_Shen, ccc_Shen, lag_Chen, ccc_Chen, lag_Gao, ccc_Gao
lonlatMadCcc = np.concatenate((lonflatten,latflatten),axis=1)
#--------------------prepare CT_MAD and CCC for specific period arange----------------------------------
def select_data(lonlatvel):
    if P.contains(Point(lonlatvel[0],lonlatvel[1])):
        plotlist.append(np.array(lonlatvel))

for period in periods:
    MADdict = loadmat('/home/cjq/ModelValid/Forward/forLove/' + 'lag_' + period + '.mat')
    cccdict = loadmat('/home/cjq/ModelValid/Forward/forLove/' + 'ccc_' + period + '.mat')
    MAD = MADdict['lag_' + str(periods.index(period))] #size: 321*301(lat*lon)
    ccc = cccdict['ccc_' + str(periods.index(period))] #size: 321*301(lat*lon)
    MADflatten = MAD.flatten().reshape(MAD.flatten().shape[0],1)
    cccflatten = ccc.flatten().reshape(ccc.flatten().shape[0],1)
    lonlatMadCcc = np.concatenate((lonlatMadCcc, MADflatten), axis=1)
    lonlatMadCcc = np.concatenate((lonlatMadCcc, cccflatten), axis=1) 
#print("lonlatMadCcc: ", lonlatMadCcc.shape)
#print("ccc: ", ccc['ccc_Gao'].shape)

#--------------select plotlist based on envelope of stations
m = Manager()
plotlist = m.list()

with concurrent.futures.ProcessPoolExecutor(2) as executor:
    executor.map(select_data, lonlatMadCcc)
plotlist = np.array(plotlist)
#print("plotlist: ", plotlist.shape)

#---------------plot MAD and ccc distribution -------------------------------------------
region = [108,123,20,36]
projection = "M?"
fig1 = pygmt.Figure() # plot MAD distribution
with fig1.subplot(nrows=1,ncols=4,figsize=("16c","5c"),margins=['0.25c']):
    for col in range(len(periods)):
        with fig1.set_panel(panel=[0,col]):
            #pygmt.config(FONT_LABEL="5p,Helvetica,red")
            pygmt.config(FONT_ANNOT_PRIMARY="5p,Helvetica")
            fig1.basemap(region=region,projection=projection,frame=["WSen","afg"])
            fig1.grdimage(grid="topo30.grd",cmap="topo.cpt",shading="+d")
            fig1.coast(region=region,projection=projection,area_thresh=[0/5000],resolution='f',water="lightblue")
            #fig1.plot(data="CN-border-La.dat",region=region,projection=projection,pen="0.2p")
            #fig1.plot(data="boundary.dat",region=region,projection=projection,pen="1p")
            #fig1.plot(data="L_NSGL.txt",region=region,projection=projection,pen="1.5p,-")
            pygmt.makecpt(cmap="jet",series=[0.5, 3.0],background=True, output="test.cpt")
            pygmt.xyz2grd(x=plotlist[:,0],y=plotlist[:,1],z=plotlist[:,2+col*2], outgrid='test.nc',region=region, spacing='0.05/0.05')
            fig1.grdimage(grid="test.nc",cmap="test.cpt",nan_transparent=True, transparency=20)
            #fig1.text(x=108,y=36.5,text=labels_mad[col],region=region, projection=projection, angle=0, justify="BR",font="6P,Helvetica-Narrow-Bold,black",fill='white',pen="black",clearance="0.05c/0.05c")

        os.remove("test.nc")
        print(periods[col], " MAD Done!")

fig1.colorbar(cmap="test.cpt", position="JCB+o0.3c/0.3c+w5c/0.2c", frame=["a0.4f0.2",'x+l"@%25%@:7p:@;black;MAD(dt/dist) (s/100km)@%%@::@;;"'])
#fig1.show()
fig1.savefig("MAD.png", dpi=400)
os.remove("test.cpt")
#------------------------------------------------------------------------------------------------------------------------------------
fig2 = pygmt.Figure() # plot ccc distribution
with fig2.subplot(nrows=1,ncols=4,figsize=("16c","5c"),margins=['0.25c']):
    for col in range(len(periods)):
        with fig2.set_panel(panel=[0,col]):
            #pygmt.config(FONT_LABEL="5p,Helvetica,red")
            pygmt.config(FONT_ANNOT_PRIMARY="5p,Helvetica")
            fig2.basemap(region=region,projection=projection,frame=["WSen","afg"])
            fig2.grdimage(grid="topo30.grd",cmap="topo.cpt",shading="+d")
            fig2.coast(region=region,projection=projection,area_thresh=[0/5000],resolution='f',water="lightblue")
            #fig2.plot(data="CN-border-La.dat",region=region,projection=projection,pen="0.2p")
            #fig2.plot(data="boundary.dat",region=region,projection=projection,pen="1p")
            #fig2.plot(data="L_NSGL.txt",region=region,projection=projection,pen="1.5p,-")
            pygmt.makecpt(cmap="jet",series=[0.6, 0.8],background=True, output="test.cpt")
            pygmt.xyz2grd(x=plotlist[:,0],y=plotlist[:,1],z=plotlist[:,2+col*2+1], outgrid='test.nc',region=region, spacing='0.05/0.05')
            fig2.grdimage(grid="test.nc",cmap="test.cpt",nan_transparent=True, transparency=20)
            #fig2.text(x=108,y=36.5,text=labels_mad[col],region=region, projection=projection, angle=0, justify="BR",font="6P,Helvetica-Narrow-Bold,black",fill='white',pen="black",clearance="0.05c/0.05c")

        os.remove("test.nc")
        print(periods[col], " ccc Done!")

fig2.colorbar(cmap="test.cpt", position="JCB+o0.3c/0.3c+w5c/0.2c", frame=["a0.05f0.01",'x+l"@%25%@:7p:@;black;Mean(cc coef) @%%@::@;;"'])
#fig2.show()
fig2.savefig("ccc.png", dpi=400)
os.remove("test.cpt")
