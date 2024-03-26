clear; clc;

%%
% nx=1500; %EW
% ny=1800; %NS       
% dx=1e3; %
% dy=1e3;

dir_file = '../wave_forward/template/parfile/';

file = load("Vs_inter.mat");
mz = file.z';

nx = 1000;
ny = 1000;
dx = 0.6;
dy = 0.6;
Xlen = nx*dx;%1500KM
Ylen = ny*dy;%1600KM
Zdep = 200;%here doesnot need to exceed bounds
%%
gridX = zeros(nx,ny);
gridY = zeros(nx,ny);
mx = linspace(0,Xlen,nx); 
my = linspace(0,Ylen,ny);
% mx = [0:nx-1]*dx;
% my = [0:ny-1]*dy;
[gridX,gridY]=meshgrid(mx,my);
gridX = gridX';
gridY = gridY';

ninterface = 2; %nlayer = ninterface-1
zlayer= zeros(ninterface);
zlayer(1)=0; 
%zlayer(2)=-15e3; 
%zlayer(3)=-25e3;
zlayer(2)=-Zdep;

disp('start write');


%write for fortran
filename = [dir_file,'/gridvmap_F.dat'];
fp1 = fopen(filename,'w');

fprintf(fp1,'%d %d %d\n',nx,ny,ninterface);

for j=1:ny
	for i=1:nx
		fprintf(fp1,'%f %f',gridX(i,j),gridY(i,j));
		for k=1:ninterface
			fprintf(fp1,' %f',zlayer(k));
		end
		fprintf(fp1,'\n');
	end
end

fclose(fp1);


%write for fortran X&Y
filename = [dir_file,'/grid_FX.dat'];
fp2 = fopen(filename,'w');
fprintf(fp2,'%d\n',nx);
for i=1:nx
	fprintf(fp2,'%f\n',mx(i));
end
fclose(fp2);

filename = [dir_file,'/grid_FY.dat'];
fp3 = fopen(filename,'w');
fprintf(fp3,'%d\n',ny);
for i=1:ny
	fprintf(fp3,'%f\n',my(i));
end
fclose(fp3);


disp('hahah');

