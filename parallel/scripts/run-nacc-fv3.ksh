#!/bin/ksh -l
#BSUB -J nacc-test
#BSUB -o jnacc_par.out1
#BSUB -e jnacc_par.err1
#BSUB -q debug
##BSUB -q dev
##BSUB -extsched "CRAYLINUX[]" -R "1*{select[craylinux && !vnode]} + 576*{select[craylinux && vnode] span [ptile=24]}"
#BSUB -M 3000
##BSUB -W 01:00
#BSUB -W 00:10
#BSUB -P CMAQ-T2O
#BSUB -extsched 'CRAYLINUX[]'
#BSUB -cwd .

#Module settings
module load PrgEnv-intel/5.2.82
module load cray-netcdf/4.3.2
module load cray-mpich/7.2.0

#Set number of nacc times  = processors, and # of nodes
NTIMES=73
export NODES=12

APPL=aqm.t12z
InMetDir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fv3gfs16_testdata
InGeoDir=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/nacc_geofiles
InVIIRSDir_GVF=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/viirs_gvf_test/grib2
InVIIRSDir_LAI=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/viirs_lai_test/
OutDir=/gpfs/hps2/ptmp/Patrick.C.Campbell/fv3gfs16_testdata/nacc_output_parallel_canopy_shade
ProgDir=/gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/NACC_canopy_shade/parallel/src

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
  file_mm    = '$InMetDir/gfs.t12z.atmf','.nc'
  file_sfc   = '$InMetDir/gfs.t12z.sfcf','.nc'
  file_geo   = '$InGeoDir/gfs.t12z.geo.01.canopy.nc'
  file_viirs_gvf = '$InVIIRSDir_GVF/GVF-WKL-GLB_v2r3_j01_s20200824_e20200830_c202008311235100.grib2.nc'
  file_viirs_lai = '$InVIIRSDir_LAI/VIIRS_VNP15A2H.001_20190829.nc'
  ioform     =  1
 &END

 &USERDEFS
  inmetmodel =  3
  dx_out      =  12000
  dy_out      =  12000
  met_cen_lat_in =  0.0
  met_cen_lon_in =  0.0
  lpv        =  0
  lwout      =  1
  luvbout    =  1
  ifdiag_pbl = .FALSE.
  ifviirs_gvf = .FALSE. 
  ifviirs_lai = .FALSE.
  iffengsha_dust = .FALSE.
  ifbioseason = .FALSE.
  ifcanopy    = .FALSE.
  mcip_start = "2020-01-12-12:00:00.0000"
  mcip_end   = "2020-01-15-13:00:00.0000"
  intvl      =  60
  coordnam   = "FV3_RPO"
  grdnam     = "FV3_CONUS"
  ctmlays    =  1.000000, 0.995253, 0.990479, 0.985679, 0.980781,
              0.975782, 0.970684, 0.960187, 0.954689, 0.936895, 
              0.930397, 0.908404, 0.888811, 0.862914, 0.829314, 
              0.786714, 0.735314, 0.645814, 0.614214, 0.582114, 
              0.549714, 0.511711, 0.484394, 0.451894, 0.419694, 
              0.388094, 0.356994, 0.326694, 0.297694, 0.270694, 
              0.245894, 0.223694, 0.203594, 0.154394, 0.127094, 0.000000
  cutlay_collapx = 22
  btrim      =  -1
  lprt_col   =  0
  lprt_row   =  0
  ntimes     = $NTIMES
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

# LSF Parallel
aprun -n$NTIMES -m8000 $ProgDir/mcip.exe
