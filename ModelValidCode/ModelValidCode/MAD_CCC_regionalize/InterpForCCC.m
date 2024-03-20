clc
clear
%% import original Vs model
aver0 = load('/home/cjq/ModelValid/Forward/forLove/aver_3_7s.txt');
aver1 = load('/home/cjq/ModelValid/Forward/forLove/aver_5_10s.txt');
aver2 = load('/home/cjq/ModelValid/Forward/forLove/aver_8_15s.txt');
aver3 = load('/home/cjq/ModelValid/Forward/forLove/aver_12_20s.txt');
%% set meshgrid for Vs model
lon0 = aver0(:,1);
lat0 = aver0(:,2);
ccc0 = aver0(:,4);

lon1 = aver1(:,1);
lat1 = aver1(:,2);
ccc1 = aver1(:,4);

lon2 = aver2(:,1);
lat2 = aver2(:,2);
ccc2 = aver2(:,4);

lon3 = aver3(:,1);
lat3 = aver3(:,2);
ccc3 = aver3(:,4);

%% interp with kriging interpolation
latmin = 20; latmax =36;
lonmin = 108; lonmax = 123;
new_lon_lat_size = 0.05;
xlon = lonmin:new_lon_lat_size:lonmax;
ylat = latmin:new_lon_lat_size:latmax;
[newlon,newlat] = meshgrid(xlon,ylat); %lat*lon
ccc_0 = zeros(length(ylat),length(xlon));
ccc_1 = zeros(length(ylat),length(xlon));
ccc_2 = zeros(length(ylat),length(xlon));
ccc_3 = zeros(length(ylat),length(xlon));

v_0 = variogram([lon0 lat0],ccc0,'plot',false,'maxdist',21);
v_1 = variogram([lon1 lat1],ccc1,'plot',false,'maxdist',21);
v_2 = variogram([lon2 lat2],ccc2,'plot',false,'maxdist',21);
v_3 = variogram([lon3 lat3],ccc3,'plot',false,'maxdist',21);

[dum,dum,dum,vstruct_0] = variogramfit(v_0.distance,v_0.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_1] = variogramfit(v_1.distance,v_1.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_2] = variogramfit(v_2.distance,v_2.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_3] = variogramfit(v_3.distance,v_3.val,[],[],[],'model','stable');

[ccc_0,Zvar] = kriging(vstruct_0,lon0,lat0,ccc0,newlon,newlat);
[ccc_1,Zvar] = kriging(vstruct_1,lon1,lat1,ccc1,newlon,newlat);
[ccc_2,Zvar] = kriging(vstruct_2,lon2,lat2,ccc2,newlon,newlat);
[ccc_3,Zvar] = kriging(vstruct_3,lon3,lat3,ccc3,newlon,newlat);
%% save
save('ccc_3_7s.mat', 'ccc_0')
save('ccc_5_10s.mat', 'ccc_1')
save('ccc_8_15s.mat', 'ccc_2')
save('ccc_12_20s.mat', 'ccc_3')
%% plot
figure()
subplot(141)
imagesc(newlon(1,:),newlat(:,1),ccc_0); axis image; axis xy
caxis([0.6,0.8])
subplot(142)
imagesc(newlon(1,:),newlat(:,1),ccc_1); axis image; axis xy
caxis([0.6,0.8])
subplot(143)
imagesc(newlon(1,:),newlat(:,1),ccc_2); axis image; axis xy
caxis([0.6,0.8])
subplot(144)
imagesc(newlon(1,:),newlat(:,1),ccc_3); axis image; axis xy
caxis([0.6,0.8])
