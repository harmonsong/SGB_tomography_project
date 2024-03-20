#!/bin/bash
for((i=2;i<=$1;i++))
do
    {
        xx=`cat virsrc.dat|head -n $i|tail -n+$i`
        nline=`grep -n '<anchor_force>' ./src$i/SeisSource.conf|awk -F: '{print $1}'`
        #echo $xx
        #echo $((nline+1))
        sed -i ''"$((nline+1))"'c '"$xx"'' ./src$i/SeisSource.conf
    }
done
