#!/bin/bash
#SBATCH --partition=contrib            # submit   to the debug partition
#SBATCH --qos=qtong            # submit   to the debug partition
#SBATCH --job-name=nacc-test          # name the job
#SBATCH --output=nacc-test-%j.out     # write stdout/stderr   to named file
#SBATCH --error=nacc-test-%j.err
#SBATCH --time=0-00:10:00             # Run for max of 00 hrs, 10 mins, 00 secs
#SBATCH --nodes=2                    # Request N nodes
#SBATCH --ntasks=2                    # Request n tasks
##SBATCH --cpus-per-task=12            # Request n cores per node
#SBATCH --mem-per-cpu=8GB             # Request nGB RAM per core

#Load more than one library at a time.
#export LMOD_EXPERT=1

#Module settings
module load gnu9 openmpi4
module load netcdf-c
module load ioapi/3.2-spack 

#Set number of nacc times  = processors, and # of nodes
NTIMES=2
export NODES=2

APPL=aqm.srw-lam
InMetDir=/groups/ESS/pcampbe8/srw-lam-test
InGeoDir=/groups/ESS/pcampbe8/nacc_geofiles
InVIIRSDir_GVF=./
InVIIRSDir_LAI=./
OutDir=/groups/ESS/pcampbe8/srw-lam-test/nacc_output_parallel
ProgDir=/groups/ESS/pcampbe8/NACC/parallel/src


if [ ! -s $InMetDir ]; then
  echo "No such input directory $InMetDir"
  exit 1
fi

if [ ! -s $InGeoDir ]; then
  echo "No such input directory $InGeoDir"
  exit 1
fi

if [ ! -d $OutDir ]; then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if [ $? != 0 ]; then
    echo "Failed to make output directory, $OutDir"
    exit 1
  fi
fi

if [ ! -d $ProgDir ]; then
  echo "No such program directory $ProgDir"
  exit 1
fi

cd $OutDir
cat>namelist.mcip<<!
&FILENAMES
  file_gd    = 'GRIDDESC'
  file_mm    = '$InMetDir/dynf','.nc'
  file_sfc   = '$InMetDir/phyf','.nc'
  file_geo   = '$InGeoDir/gfs.t12z.geo.06.nc' !needed for valid LAI/LANDUSEF for CMAQ
  file_viirs_gvf = '$InVIIRSDir_GVF/'
  file_viirs_lai = '$InVIIRSDir_LAI/'
  ioform     =  1
 &END

 &USERDEFS
  inmetmodel =  4
  dx_out      =  12000
  dy_out      =  12000
!  met_cen_lat_in =  0.0
!  met_cen_lon_in =  0.0
  lpv        =  0
  lwout      =  1
  luvbout    =  1
  ifdiag_pbl = .FALSE.
  ifviirs_gvf = .FALSE.
  ifviirs_lai = .FALSE.
  iffengsha_dust = .FALSE.
  ifbioseason = .FALSE.
  ifcanopy    = .FALSE.
  eradm       = 6378388.0 ! maybe needed for slightly different SRW-LAM rearth
  mcip_start = "2019-06-15-00:00:00.0000"
  mcip_end   = "2019-06-15-01:00:00.0000"
  intvl      =  60
  coordnam   = "FV3_RPO"
  grdnam     = "FV3_CONUS"
  ctmlays    =  1.000000, 0.994670, 0.990479, 0.985679, 0.980781,
              0.975782, 0.970684, 0.960187, 0.954689, 0.936895,
              0.930397, 0.908404, 0.888811, 0.862914, 0.829314,
              0.786714, 0.735314, 0.645814, 0.614214, 0.582114,
              0.549714, 0.511711, 0.484394, 0.451894, 0.419694,
              0.388094, 0.356994, 0.326694, 0.297694, 0.270694,
              0.245894, 0.223694, 0.203594, 0.154394, 0.127094, 0.000000
  cutlay_collapx = 0
  btrim      =  -1
  lprt_col   =  0
  lprt_row   =  0
  ntimes     = 2
  projparm = 2., 33.,45., -97., -97., 40.
  domains = -2508000., -1716000., 12000., 12000., 442, 265
 &END

 &WINDOWDEFS
  x0         =  1
  y0         =  1
  ncolsin    =  442
  nrowsin    =  265
 &END
!

export IOAPI_CHECK_HEADERS=T
export EXECUTION_ID=$PROG

export GRID_BDY_2D=${APPL}.grdbdy2d.ncf
export GRID_CRO_2D=${APPL}.grdcro2d.ncf
export GRID_DOT_2D=${APPL}.grddot2d.ncf
export MET_BDY_3D=${APPL}.metbdy3d.ncf
export MET_CRO_2D=${APPL}.metcro2d.ncf
export MET_CRO_3D=${APPL}.metcro3d.ncf
export MET_DOT_3D=${APPL}.metdot3d.ncf
export LUFRAC_CRO=${APPL}.lufraccro.ncf
export SOI_CRO=${APPL}.soicro.ncf
export MOSAIC_CRO=${APPL}.mosaiccro.ncf

rm -f *.ncf 

#Slurm Parallel
#srun -n$NTIMES $ProgDir/mcip.exe
module load ucx libfabric
#export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
prun $ProgDir/mcip.exe
