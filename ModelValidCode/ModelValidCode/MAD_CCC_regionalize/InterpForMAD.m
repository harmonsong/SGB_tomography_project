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
lag0 = aver0(:,3);

lon1 = aver1(:,1);
lat1 = aver1(:,2);
lag1 = aver1(:,3);

lon2 = aver2(:,1);
lat2 = aver2(:,2);
lag2 = aver2(:,3);

lon3 = aver3(:,1);
lat3 = aver3(:,2);
lag3 = aver3(:,3);

%% interp with kriging interpolation
latmin = 20; latmax =36;
lonmin = 108; lonmax = 123;
new_lon_lat_size = 0.05;
xlon = lonmin:new_lon_lat_size:lonmax;
ylat = latmin:new_lon_lat_size:latmax;
[newlon,newlat] = meshgrid(xlon,ylat); %lat*lon
lag_0 = zeros(length(ylat),length(xlon));
lag_1 = zeros(length(ylat),length(xlon));
lag_2 = zeros(length(ylat),length(xlon));
lag_3 = zeros(length(ylat),length(xlon));

v_0 = variogram([lon0 lat0],lag0,'plot',false,'maxdist',21);
v_1 = variogram([lon1 lat1],lag1,'plot',false,'maxdist',21);
v_2 = variogram([lon2 lat2],lag2,'plot',false,'maxdist',21);
v_3 = variogram([lon3 lat3],lag3,'plot',false,'maxdist',21);

[dum,dum,dum,vstruct_0] = variogramfit(v_0.distance,v_0.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_1] = variogramfit(v_1.distance,v_1.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_2] = variogramfit(v_2.distance,v_2.val,[],[],[],'model','stable');
[dum,dum,dum,vstruct_3] = variogramfit(v_3.distance,v_3.val,[],[],[],'model','stable');

[lag_0,Zvar] = kriging(vstruct_0,lon0,lat0,lag0,newlon,newlat);
[lag_1,Zvar] = kriging(vstruct_1,lon1,lat1,lag1,newlon,newlat);
[lag_2,Zvar] = kriging(vstruct_2,lon2,lat2,lag2,newlon,newlat);
[lag_3,Zvar] = kriging(vstruct_3,lon3,lat3,lag3,newlon,newlat);
%% save
save('lag_3_7s.mat', 'lag_0')
save('lag_5_10s.mat', 'lag_1')
save('lag_8_15s.mat', 'lag_2')
save('lag_12_20s.mat', 'lag_3')
%% plot
figure()
subplot(141)
imagesc(newlon(1,:),newlat(:,1),lag_0); axis image; axis xy
caxis([0.5,3.0])
subplot(142)
imagesc(newlon(1,:),newlat(:,1),lag_1); axis image; axis xy
caxis([0.5,3.0])
subplot(143)
imagesc(newlon(1,:),newlat(:,1),lag_2); axis image; axis xy
caxis([0.5,3.0])
subplot(144)
imagesc(newlon(1,:),newlat(:,1),lag_3); axis image; axis xy
caxis([0.5,3.0])
