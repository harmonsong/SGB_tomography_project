clear; clc;
% load('Staloc.dat')
% Stalon = Staloc(:,1);
% Stalat = Staloc(:,2);
Staloc = load("sta_in.mat");
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
hold on
plot([staposx(400)/1000, staposx(100)/1000], [staposy(400)/1000, staposy(100)/1000])
grid minor
%% write
filename = 'stapos.dat';
fileid = fopen(filename,'w');
for i = 1:nsta
	sta = [staposx(i); staposy(i)];
% 	fprintf(fileid,'recv_%3.3d = %-f	%-f 	9000.0E3	%s\n',i,aft(1), aft(2), Staname{i});
    fprintf(fileid,'recv_%3.3d = %-f	%-f 	9000.0E3 \n',i,sta(1), sta(2));
end
fclose(fileid);

disp('haha');
