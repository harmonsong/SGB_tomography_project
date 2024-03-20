#!/bin/bash
for((i=1;i<=60;i++))
do
    {
        cd src$i/
        pwd
        ./pyplot/figsgf6.py -s $i -d 5
        ./pyplot/misfit_cal.py -s $i -d 1
        cd ../
    }
done
