#!/bin/bash
for((i=1;i<=60;i++))
do
    {
        #cp -r src1/pyplot src$i/pyplot
        #cp src1/sc.p4.mplstyle src$i/
        #cp  src1/pyplot/misfit_cal2.py src$i/pyplot/
        #cp  src1/pyplot/figsgf7.py src$i/pyplot/
        #cp  src1/runALL.sh src$i/
        cp -rv template src$i
    }
done
