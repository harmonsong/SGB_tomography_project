#!/bin/bash
for((i=1;i<=60;i++))
do
    {
        cd src$i/
        /work/ess-wangp/FD3DtopoEw/src/bin/seis3d_source SeisFD3D.conf
        pwd
        cd ../
    }
done
