clc
clear
%% import original Vs model
VsModStruct = load('/home/cjq/ModelValid/Velmod/Chen_Vs_3D_smooth.mat');
VsMod = getfield(VsModStruct,'Vs'); %size: lat*lon*dep
%% set meshgrid for Vs model
latmin = 20; latmax = 36;
lonmin = 108; lonmax = 123;
lon_lat_size = 0.5;
dep_size = [3,5,10];

grid_lon = lonmin:lon_lat_size:lonmax;
grid_lat = latmin:lon_lat_size:latmax;
grid_dep0 = 0:dep_size(1):75;
grid_dep1 = 75:dep_size(2):100;
grid_dep2 = 100:dep_size(3):150;
grid_dep = [grid_dep0, grid_dep1(2:end), grid_dep2(2:end)];
fprintf('grid_lat = %2d\n',length(grid_lat))
fprintf('grid_lon = %2d\n',length(grid_lon))
fprintf('grid_dep = %2d\n',length(grid_dep))
[lon_cun,lat_cun] = meshgrid(grid_lon,grid_lat); %lat*lon*dep
%% interp first along horizontal direction with kriging interpolation
new_lon_lat_size = lon_lat_size/20;
xlon = lonmin:new_lon_lat_size:lonmax;
ylat = latmin:new_lon_lat_size:latmax;
[newlon,newlat] = meshgrid(xlon,ylat); %lat*lon
Vs_h = zeros(length(ylat),length(xlon),length(grid_dep));

for i=1:length(grid_dep)
    vs=VsMod(:,:,i);
%     subplot(2,2,1)
%     imagesc(newlon(1,:),newlat(:,1),VsMod(:,:,1)); axis image; axis xy
%     title('Original Vs model')
%     subplot(2,2,2)
    v = variogram([lon_cun(:) lat_cun(:)],vs(:),'plot',false,'maxdist',21);
%     title('variogram')
    [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
    [Zhat,Zvar] = kriging(vstruct,lon_cun(:),lat_cun(:),vs(:),newlon,newlat);
    Vs_h(:,:,i) = Zhat;
    fprintf('%2d Done!\n',i)
%     subplot(2,2,3)
%     imagesc(newlon(1,:),newlat(:,1),Zhat); axis image; axis xy
%     title('Vs model after the kriging interpolation')
%     subplot(2,2,4)
%     contour(newlon,newlat,Zvar); axis image
%     title('kriging variance')
end
%% interp then along depth direction 
newdep_size = dep_size/5;
newdep0 = 0:newdep_size(1):depBound(1);
newdep1 = depBound(1):newdep_size(2):depBound(2);
newdep2 = depBound(2):newdep_size(3):depBound(3);
newdep = [newdep0, newdep1(2:end),newdep2(2:end)];
Vs_v = zeros(length(ylat),length(xlon),length(newdep));
for ix = 1: length(xlon)
    for iy = 1:length(ylat)
        vs0 = Vs_h(iy,ix,:);
        Vs_v(iy,ix,:) = interp1(grid_dep,vs0(:),newdep,'spline');
    end
end
%% plot Vs model after interpolation
Vs_full = Vs_v*1000;
[tmplon,tmpdep] = meshgrid(xlon,newdep); %lon*dep
tmpVs = squeeze(Vs_full(:,100,:));
imagesc([108 123],[0 150],tmpVs');
colormap(flip(jet));
xlabel('Longitude')
ylabel('Depth (km)')
title('Vs model after interpolation')
