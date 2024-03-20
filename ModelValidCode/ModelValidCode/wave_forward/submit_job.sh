#!/bin/bash
for((SEQ=1;SEQ<=10;SEQ++))
do

        cd src$SEQ
	pwd
	bsub<runALL.sh

	cd ../

	
done
