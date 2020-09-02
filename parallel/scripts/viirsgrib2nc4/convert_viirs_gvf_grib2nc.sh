#!/bin/bash

execdir=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/viirs_gvf_test

base=${1}
for dir in ${base}/*; do
    echo "*******************************************************"
    echo ${dir}
    echo "*******************************************************"
    if [ -d ${dir} ]; then 
	cd ${dir}
	ls -t GVF-WKL-GLB_v2r3_* | xargs -I {} --max-procs 10 $execdir/viirsgrib2nc4.py -vf {}
	cd ${base}
    fi
done
