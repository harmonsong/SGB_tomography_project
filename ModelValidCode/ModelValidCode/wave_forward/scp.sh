#!/bin/bash
for ((i=1;i<=60;i++))
#for ((i=$1;i<=$1;i++))
do
{
    cd src$i/
    scp -P 77 -r output SeisFD3D.conf harmon@10.16.21.134:~/data/git_repo/SGB_tomography_project/project_repartition_v4.0/output_repar_v9.5_01--10-16Hz/ModelValidate/src$i/
    cd ../
    pwd
}
done
