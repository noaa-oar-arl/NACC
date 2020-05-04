#!/bin/ksh
#BSUB -J jnacc_test_v1
#BSUB -o /gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/NACC/serial/scripts/jnacc_test.out
#BSUB -e /gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/NACC/serial/scripts/jnacc_test.err
##BSUB -q debug
#BSUB -q dev
##BSUB -extsched "CRAYLINUX[]" -R "1*{select[craylinux && !vnode]} + 576*{select[craylinux && vnode] span [ptile=24]}"
#BSUB -M 3000
##BSUB -W 01:00
#BSUB -W 00:10
#BSUB -P CMAQ-T2O
#BSUB -extsched 'CRAYLINUX[]'

#%include <head.h>
#%include <envir-xc40.h>
#
module load PrgEnv-intel/5.2.82
module load cray-netcdf/4.3.2

export PROCS=1
export NODES=1
export IOBUF_PARAMS="128M:verbose"

./run-nacc-fv3.ksh
