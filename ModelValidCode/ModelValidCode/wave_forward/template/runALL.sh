#!/bin/bash
#BSUB -L /bin/bash
#BSUB -J RAwave
#BSUB -q medium
#BSUB -n 240
#BSUB -R "span[ptile=40]"
#BSUB -oo sc.out
#BSUB -eo sc.err




# load nessary modules

module load ~/rupture
#module load intel/2018.4
#module load mpi/mpich/3.3-intel-18.4 
#module load netcdf-fortran/4.4.4-intel-18.4

# run your program
ulimit -s unlimited
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpirun -np $LSB_DJOB_NUMPROC /work/ess-wangp/FD3DtopoEw/src/bin/seis3d_grid_mpi SeisFD3D.conf
mpirun -np $LSB_DJOB_NUMPROC /work/ess-wangp/FD3DtopoEw/src/bin/seis3d_metric_mpi SeisFD3D.conf
mpirun -np $LSB_DJOB_NUMPROC /work/ess-wangp/FD3DtopoEw/src/bin/seis3d_media_mpi SeisFD3D.conf
#/work/ess-wangp/FD3DtopoEw/src/bin/seis3d_source SeisFD3D.conf
#/work/ess-wangp/FD3DtopoEw/src/bin/seis3d_station SeisFD3D.conf
#mpirun -np $LSB_DJOB_NUMPROC /work/ess-wangp/FD3DtopoEw/src/bin/seis3d_wave_mpi SeisFD3D.conf
