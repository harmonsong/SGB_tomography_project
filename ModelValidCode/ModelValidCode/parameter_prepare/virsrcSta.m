clear; clc;
% load('virsrcpos.dat')
% Stalon = virsrcpos(:,1);
% Stalat = virsrcpos(:,2);
Staloc = load("sta_virsrc.mat");
Stalon = Staloc.lon_sta';
Stalat = Staloc.lat_sta';
%%
file = load("Vs_inter_smooth.mat");
Vs_full = file.vs_inter;
lon = file.lon_sta;
lat = file.lat_sta;
mz = file.z;

Clat1 = min(lat); Clat2 = max(lat);
Xlon = min(lon); Ylat = min(lat);
lonunit = 2*pi*6359.752/360*cos((Clat1+Clat2)/2/180*pi)*1e3; % unit length of longitude                                    
latunit = 2*pi*6378.137/360*1e3; % unit length of latitude                        

nsta = length(Stalon);

staposx = zeros(nsta,1);
staposy = zeros(nsta,1);

staposx = (Stalon - Xlon)*lonunit;
staposy = (Stalat - Ylat)*latunit;
%% plot
scatter(staposx/1000,staposy/1000)

%% write
filename = 'virsrc.dat';
fileid = fopen(filename,'w');
for i = 1:nsta
	sta = [staposx(i); staposy(i)];
    fprintf(fileid,'%-f	%-f 4.5e3 0.0   1.0e+22   0.0 0.0 1.0\n',sta(1), sta(2));
end
fclose(fileid);

disp('haha');
