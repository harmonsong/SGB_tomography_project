clc
clear
%% import original Vs model
VsModStruct = load('/home/cjq/ModelValid/Velmod/Chen_Vsh_3D_smooth.mat');
VsMod = getfield(VsModStruct,'Vs'); %size: lat*lon*dep
%% set meshgrid for Vs model
latmin = 20; latmax = 36;
lonmin = 108; lonmax = 123;
lon_lat_size = 0.05;
dep_size = 3;

grid_lon = lonmin:lon_lat_size:lonmax;
grid_lat = latmin:lon_lat_size:latmax;
grid_dep = 0:dep_size:75;
fprintf('grid_lat = %2d\n',length(grid_lat))
fprintf('grid_lon = %2d\n',length(grid_lon))
fprintf('grid_dep = %2d\n',length(grid_dep))
[lon_cun,lat_cun] = meshgrid(grid_lon,grid_lat); %lat*lon*dep
%% interp first along horizontal direction
new_lon_lat_size = lon_lat_size/2;
xlon = lonmin:new_lon_lat_size:lonmax;
ylat = latmin:new_lon_lat_size:latmax;
[newlon,newlat] = meshgrid(xlon,ylat); %lat*lon
Vs_h = zeros(length(ylat),length(xlon),length(grid_dep));
for j = 1:length(grid_dep)
    Vs_h(:,:,j) = interp2(lon_cun, lat_cun, VsMod(:,:,j),newlon, newlat);
    fprintf('%2d Done!\n',j)
end

%% interp then along depth direction
newdep_size = dep_size/5;
newdep = 0:newdep_size:75;
Vs_v = zeros(length(ylat),length(xlon),length(newdep));
for ix = 1: length(xlon)
    for iy = 1:length(ylat)
        vs0 = Vs_h(iy,ix,:);
        Vs_v(iy,ix,:) = interp1(grid_dep,vs0(:),newdep,'spline');
    end
end
%% plot Vs model after interpolation
Vs_full = Vs_v*1000;
save("Vsh_full_75km.mat", 'Vs_full')
%%
[tmplon,tmpdep] = meshgrid(xlon,newdep); %lon*dep
tmpVs = squeeze(Vs_full(600,:,:));
imagesc([108 123],[0 75],tmpVs');
colormap(flip(jet));
xlabel('Longitude')
ylabel('Depth (km)')
title('Vs model after interpolation')
