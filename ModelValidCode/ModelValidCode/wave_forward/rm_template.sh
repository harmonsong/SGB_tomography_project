#!/bin/bash
for ((i=1;i<=63;i++))
do
{   
    cd src$i/
    pwd
    rm -rf input/
    rm -rf bin/
    rm -rf checkpoint/
    rm -rf parfile/
    rm -rf pyplot/
    rm *.dat
    rm *.log
    rm *.mplstyle
    rm *.sh
    rm *.py
    cd ../
}
done