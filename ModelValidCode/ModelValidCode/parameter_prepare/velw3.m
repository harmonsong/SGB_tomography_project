clear; clc;

storepath='./pictures/';

% %vel: size(lat*lon*dep)
% load('newdep.mat');
% load('Vs_full.mat');
% mx = linspace(0,1500e3,601);%lon
% my = linspace(0,1800e3,641);%lat
% mz = -1*newdep*1e3; %dep
file = load("Vs_inter.mat");
Vs_full = file.vs_inter;
mx = file.x;
my = file.y;
mz = -file.z;

dir_file = '../wave_forward/template/parfile/';

%% get Vall(vs), unit km/s
    %vall nx ny nz,  lon,lat,dep
    %Vs_full lat*lon*dep
Vall = permute(Vs_full, [1,2,3]);
nz=length(mz);ny=length(my);nx=length(mx);
%% construct velocity array

if Vall(end,1,1) > 100
	Vall = Vall/1e3; %transfer to km/s
end

Vs = Vall;
clear Vall;


Vp = 0.9409+2.0947.*Vs-0.8206.*Vs.^2+0.2683.*Vs.^3-0.0251.*Vs.^4;
rho  = (1.6612.*Vp-0.4721.*Vp.^2+0.0671.*Vp.^3-0.0043.*Vp.^4+0.000106.*Vp.^5); % in grams/cubic     centimeter

Vs = Vs*1e3;
Vp = Vp*1e3;
rho = rho*1e3;

velmax = max(max(max(Vp)));
velmin = min(min(min(Vs)));



%% write

disp('start to write~');


%write fortran
ncid = netcdf.create([dir_file,'Vsim_fortran_new.nc'],'netcdf4');%rewrite

dimid = ones(1,3);%fast mid low
dimid(1) = netcdf.defDim(ncid,'x',nx);
dimid(2) = netcdf.defDim(ncid,'y',ny);
dimid(3) = netcdf.defDim(ncid,'depth',nz);

varid = ones(1,7);
varid(1) = netcdf.defVar(ncid, 'x', 'NC_FLOAT', dimid(1));
varid(2) = netcdf.defVar(ncid, 'y', 'NC_FLOAT', dimid(2));
varid(3) = netcdf.defVar(ncid, 'depth', 'NC_FLOAT', dimid(3));

varid(4) = netcdf.defVar(ncid, 'depth2sealevel', 'NC_FLOAT', dimid(3));
varid(5) = netcdf.defVar(ncid, 'Vp', 'NC_FLOAT', dimid);
varid(6) = netcdf.defVar(ncid, 'Vs', 'NC_FLOAT', dimid);
varid(7) = netcdf.defVar(ncid, 'rho', 'NC_FLOAT', dimid);

varid5 = netcdf.getConstant('Global');
netcdf.putAtt(ncid, varid5, 'sealevel', 0.0);
varid6 = netcdf.getConstant('Global');
netcdf.putAtt(ncid, varid6, 'velmax', velmax);
varid7 = netcdf.getConstant('Global');
netcdf.putAtt(ncid, varid7, 'velmin', velmin);

netcdf.endDef(ncid);

netcdf.putVar(ncid, varid(1), [0],[nx],[1], mx);
netcdf.putVar(ncid, varid(2), [0],[ny],[1], my);
netcdf.putVar(ncid, varid(3), [0],[nz],[1], mz);
netcdf.putVar(ncid, varid(4), [0],[nz],[1], -mz);

%(nz,nx,ny)
% VAR = permute(Vp, [2,3,1]);
% netcdf.putVar(ncid, varid(5), VAR);
% VAR = permute(Vs, [2,3,1]);
% netcdf.putVar(ncid, varid(6), VAR);
% VAR = permute(rho, [2,3,1]);
% netcdf.putVar(ncid, varid(7), VAR);

netcdf.putVar(ncid, varid(5), Vp);
netcdf.putVar(ncid, varid(6), Vs);
netcdf.putVar(ncid, varid(7), rho);

netcdf.close(ncid)

fprintf('pass for_output\n');






