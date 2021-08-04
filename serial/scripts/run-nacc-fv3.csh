#!/bin/csh -f
set APPL=aqm.t12z
set InMetDir=/groups/ESS/pcampbe8/fv3gfs16_testdata
set InGeoDir=/groups/ESS/pcampbe8/nacc_geofiles
set InVIIRSDir=$InMetDir
set OutDir=/groups/ESS/pcampbe8/fv3gfs16_testdata/nacc_output
set ProgDir=/groups/ESS/pcampbe8/NACC/serial/src

if ( ! -s $InMetDir ) then
  echo "No such input directory $InMetDir"
  exit 1
endif

if ( ! -s $InGeoDir ) then
  echo "No such input directory $InGeoDir"
  exit 1
endif

if ( ! -d $OutDir ) then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if ( $? != 0 ) then
    echo "Failed to make output directory, $OutDir"
    exit 1
  endif
endif

if ( ! -d $ProgDir ) then
  echo "No such program directory $ProgDir"
  exit 1
endif

cd $OutDir

cat>namelist.mcip<<!
&FILENAMES
  file_gd    = 'GRIDDESC'
  file_mm    = '$InMetDir/gfs.t12z.atmf','.nc'
  file_sfc   = '$InMetDir/gfs.t12z.sfcf','.nc'
  file_geo   = '$InGeoDir/gfs.t12z.geo.01.nc'
  file_viirs_gvf = '$InVIIRSDir/GVF-WKL-GLB_v2r3_j01_s20200824_e20200830_c202008311235100.grib2.nc'
  file_viirs_lai = '$InVIIRSDir/VIIRS_VNP15A2H.001_20190829.nc'
  ioform     =  1
 &END

 &USERDEFS
  inmetmodel =  3
  dx_in      =  12000
  dy_in      =  12000
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
  ntimes     =  73
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

setenv IOAPI_CHECK_HEADERS T

setenv GRID_BDY_2D ${APPL}.grdbdy2d.ncf
setenv GRID_CRO_2D ${APPL}.grdcro2d.ncf
setenv GRID_DOT_2D ${APPL}.grddot2d.ncf
setenv MET_BDY_3D  ${APPL}.metbdy3d.ncf
setenv MET_CRO_2D  ${APPL}.metcro2d.ncf
setenv MET_CRO_3D  ${APPL}.metcro3d.ncf
setenv MET_DOT_3D  ${APPL}.metdot3d.ncf
setenv LUFRAC_CRO  ${APPL}.lufraccro.ncf
setenv SOI_CRO     ${APPL}.soicro.ncf
setenv MOSAIC_CRO  ${APPL}.mosaiccro.ncf

rm -f *.ncf

# Serial
$ProgDir/mcip.exe

# LSF Serial
#aprun -n${PROCS} -N${NODES} $ProgDir/mcip.exe
