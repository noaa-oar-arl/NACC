#!/bin/csh -f 

#------------------------------------------------------------------------------#
#  The Community Multiscale Air Quality (CMAQ) system software is in           #
#  continuous development by various groups and is based on information        #
#  from these groups: Federal Government employees, contractors working        #
#  within a United States Government contract, and non-Federal sources         #
#  including research institutions.  These groups give the Government          #
#  permission to use, prepare derivative works of, and distribute copies       #
#  of their work in the CMAQ system to the public and to permit others         #
#  to do so.  The United States Environmental Protection Agency                #
#  therefore grants similar permission to use the CMAQ system software,        #
#  but users are requested to provide copies of derivative works or            #
#  products designed to operate in the CMAQ system to the United States        #
#  Government without restrictions as to use by others.  Software              #
#  that is used with the CMAQ system but distributed under the GNU             #
#  General Public License or the GNU Lesser General Public License is          #
#  subject to their copyright restrictions.                                    #
#------------------------------------------------------------------------------#

#=======================================================================
#
#  Script:  run.mcip
#  Purpose: Runs Models-3/CMAQ Meteorology-Chemistry Interface
#           Processor.  Part of the US EPA's Models-3/CMAQ system.
#  Method:  In UNIX/Linux:  run.mcip >&! mcip.log
#  Revised: 20 Sep 2001  Original version.  (T. Otte)
#           18 Oct 2001  Added CoordName to user definitions.  Deleted
#                        script variable DomIdMM5.  Added Fortran link
#                        for GRIDDESC file.  Moved namelist output to
#                        WorkDir, and mmheader output to OutDir.  Added
#                        user variables I0, J0, NCOLS, and NROWS for
#                        MCIP windowing.  (T. Otte)
#           29 Jan 2002  Added new namelist for file names.  Generalized
#                        the end-of-namelist delimiter.  (T. Otte)
#           27 Feb 2002  Removed minimum size for windows.  (T. Otte)
#           19 Mar 2002  Changed default grid cell for printing.
#                        (T. Otte)
#           11 Jun 2003  Clarified instructions on use of BTRIM and
#                        setting I0 and J0 for windowing option.
#                        Removed GRIDBDY2D, GRIDBDY3D, and METBDY2D
#                        from output.  (T. Otte)
#           01 Jul 2004  Restored GRIDBDY2D to output.  (T. Otte)
#           29 Nov 2004  Added TERRAIN option for input to get
#                        fractional land use from MM5 preprocessor.
#                        (T. Otte)
#           26 May 2005  Changed I0 and J0 to Y0 and X0 to make code
#                        more general.  Removed "_G1" from environment
#                        variables for output files.  Created two new
#                        user options for calculating dry deposition
#                        velocities.  Added capability to process more
#                        than five input meteorology files in a single
#                        MCIP run.  (T. Otte)
#           27 Feb 2006  Updated automated namelist generator for
#                        Linux on Mac (assumed to be) using the XLF
#                        compiler.  (T. Otte)
#           24 Jul 2007  Added option to bypass dry deposition velocity
#                        calculations in MCIP so that they can be done
#                        inline in the CCTM.  Eliminated options to
#                        use RADM (Wesely) dry deposition, eliminated
#                        multiple versions of M3Dry (Pleim) dry
#                        deposition, and eliminated options and to
#                        recalculate PBL and radiation fields in MCIP.
#                        (T. Otte)
#           27 May 2008  Added optional namelist variable to override
#                        earth radius default from MM5 and WRF.  
#                        (T. Otte)
#                        Added variables to support GOES satellite
#                        cloud processing (InSatDir, InSatFile, LSAT).
#                        Requires additional data and preprocessing
#                        package available from University of Alabama
#                        at Huntsville.  Contributed by University of
#                        Alabama at Huntsville.  (A. Biazar and T. Otte)
#           23 Dec 2008  Added optional namelist variable to override
#                        default setting for reference latitude for
#                        WRF Lambert conformal projection.  (T. Otte)
#           19 Mar 2010  Added namelist variable option to compute
#                        and output potential vorticity.  Added namelist
#                        variable option to output vertical velocity
#                        predicted by meteorological model.  Allow
#                        output from WRF Preprocessing System (WPS)
#                        routine, GEOGRID, to provide fractional land
#                        use output if it is unavailable in WRF output.
#                        Add user option to output u- and v-component
#                        winds on C-staggered grid.  (T. Otte)
#           09 Sep 2010  Removed option to generate dry deposition
#                        velocities in MCIP.  (T. Otte)
#           07 Sep 2011  Corrected minor typos in error-checking (as
#                        identified by Debra Baker, Univ. of Maryland).
#                        Updated disclaimer.  (T. Otte)
#           31 May 2012  Changed comment about MAX_MM to be consistent
#                        with the change to the code.  (T. Otte)
#           16 Mar 2018  Added new optional output files for land use,
#                        soil, and mosaic output.  Now delete rather
#                        than overwrite existing MCIP output files.
#                        (T. Spero)
#           18 Dec 2018  Removed support for MM5v3 input.  Added runtime
#                        option to choose output format.  Removed option
#                        to turn off static output.  (T. Spero)
#           20 Jun 2019  Removed layer collapsing.  Changed LUVCOUT to
#                        to LUVBOUT to make the default output for
#                        u- and v-component winds on the Arakawa-C
#                        staggering.  The Arakawa-B staggering is now
#                        optional (additional fields), and the Arakawa-C
#                        staggering is the default.  (T. Spero)
#=======================================================================

#-----------------------------------------------------------------------
# Set identification for input and output files.
#
#   APPL       = Application Name (tag for MCIP output file names)
#   CoordName  = Coordinate system name for GRIDDESC
#   GridName   = Grid Name descriptor for GRIDDESC
#   InMetDir   = Directory that contains input meteorology files
#   InGeoDir   = Directory that contains input WRF "GEOGRID" file to
#                provide fractional land-use categories if "LANDUSEF"
#                was not included in the WRFOUT files.
#   OutDir     = Directory to write MCIP output files
#   ProgDir    = Directory that contains the MCIP executable
#   WorkDir    = Working Directory for Fortran links and namelist
#-----------------------------------------------------------------------

#source $CMAQ_HOME/config_cmaq.csh

set APPL       = aqm.t12z
set CoordName  = WRF_RPO    # 16-character maximum
set GridName   = WRF_CONUS        # 16-character maximum

#set DataPath   = $CMAQ_DATA
set InMetDir   = /gpfs/hps2/ptmp/Patrick.C.Campbell/wrf_nacc_test
set InGeoDir   = $InMetDir
set OutDir     = /gpfs/hps2/ptmp/Patrick.C.Campbell/wrf_nacc_test/output
set ProgDir    = /gpfs/hps3/emc/naqfc/noscrub/Patrick.C.Campbell/CMAQ_REPO/PREP/mcip/src
set WorkDir    = $OutDir

#-----------------------------------------------------------------------
# Set name(s) of input meteorology file(s)
#
#   File name(s) must be set inside parentheses since "InMetFiles" is
#   a C-shell script array.  Multiple file names should be space-
#   delimited.  Additional lines can be used when separated by a
#   back-slash (\) continuation marker.  The file names can be as
#   they appear on your system; MCIP will link the files in by a
#   Fortran unit number and the explicit name via a namelist.  The
#   files must be listed in chronological order.  The maximum number
#   of input meteorology files must be less than or equal to the number
#   in MAX_MM in file_mod.F (default is 367).
#
#   Example:
#     set InMetFiles = ( $InMetDir/wrfout_d01_date1 \
#                        $InMetDir/wrfout_d01_date2 )
#
#-----------------------------------------------------------------------

set InMetFiles = ( $InMetDir/wrfout_d01_2017-12-01_00:00:00 \
                   $InMetDir/wrfout_d01_2017-12-01_01:00:00 \
                   $InMetDir/wrfout_d01_2017-12-01_02:00:00 )

set IfGeo      = "F"
set InGeoFile  = $InGeoDir/geo_em_d01.nc

#-----------------------------------------------------------------------
# Additional NACC Settings that are required to specify model
# Currently Supported models:
# WRF     = 2
# FV3-GFS = 3
set IN_MET_MODEL = 2

# Set exact number of times processed in NACC equal to MCIP 
# time period below
set N_TIMES = 2
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Set user control options.
#
#   LPV:     0 = Do not compute and output potential vorticity
#            1 = Compute and output potential vorticity
#
#   LWOUT:   0 = Do not output vertical velocity
#            1 = Output vertical velocity
#
#   LUVBOUT: 0 = Do not output u- and v-component winds on B-grid
#            1 = Output u- and v-component winds on B-grid (cell corner)
#                in addition to the C-grid (cell face) output
#-----------------------------------------------------------------------

set LPV     = 0
set LWOUT   = 0
set LUVBOUT = 1

#-----------------------------------------------------------------------
# Set run start and end date.  (YYYY-MO-DD-HH:MI:SS.SSSS)
#   MCIP_START:  First date and time to be output [UTC]
#   MCIP_END:    Last date and time to be output  [UTC]
#   INTVL:       Frequency of output [minutes]
#-----------------------------------------------------------------------

set MCIP_START = 2017-12-01-01:00:00.0000  # [UTC]
set MCIP_END   = 2017-12-01-02:00:00.0000  # [UTC]

set INTVL      = 60 # [min]


#-----------------------------------------------------------------------
# Set CTM layers.  Should be in descending order starting at 1 and
# # ending with 0.  There is currently a maximum of 100 layers allowed.
# # To use all of the layers from the input meteorology without
# # collapsing (or explicitly specifying), set CTMLAYS = -1.0.
# #-----------------------------------------------------------------------
#
set CTMLAYS = "-1.0"
#
# #set CTMLAYS = "1.000, 0.995, 0.990, 0.980, 0.960, 0.940, 0.910, 0.860, \
# #               0.800, 0.740, 0.650, 0.550, 0.400, 0.200, 0.000"
#
#set CTMLAYS =  "1.000000, 0.995253, 0.990479, 0.985679, 0.980781, 0.975782, 0.970684, 0.960187, 0.954689, 0.936895, 0.930397, 0.908404, 0.888811, 0.862914, 0.829314, 0.786714, 0.735314, 0.645814, 0.614214, 0.582114,0.549714, 0.511711, 0.484394, 0.451894, 0.419694, 0.388094, 0.356994, 0.326694, 0.297694, 0.270694, 0.245894, 0.223694, 0.203594, 0.154394, 0.127094, 0.000000"

#-----------------------------------------------------------------------
# Choose output format.
#   1 = Models-3 I/O API
#   2 = netCDF
#-----------------------------------------------------------------------

set IOFORM = 1

#-----------------------------------------------------------------------
# Set number of meteorology "boundary" points to remove on each of four
# horizontal sides of MCIP domain.  This affects the output MCIP domain
# dimensions by reducing meteorology domain by 2*BTRIM + 2*NTHIK + 1,
# where NTHIK is the lateral boundary thickness (in BDY files), and the
# extra point reflects conversion from grid points (dot points) to grid
# cells (cross points).  Setting BTRIM = 0 will use maximum of input
# meteorology.  To remove MM5 lateral boundaries, set BTRIM = 5.
#
# *** If windowing a specific subset domain of input meteorology, set
#     BTRIM = -1, and BTRIM will be ignored in favor of specific window
#     information in X0, Y0, NCOLS, and NROWS.
#-----------------------------------------------------------------------

set BTRIM = 0

#-----------------------------------------------------------------------
# Define MCIP subset domain.  (Only used if BTRIM = -1.  Otherwise,
# the following variables will be set automatically from BTRIM and
# size of input meteorology fields.)
#   X0:     X-coordinate of lower-left corner of full MCIP "X" domain
#           (including MCIP lateral boundary) based on input MM5 domain.
#           X0 refers to the east-west dimension.  Minimum value is 1.
#   Y0:     Y-coordinate of lower-left corner of full MCIP "X" domain
#           (including MCIP lateral boundary) based on input MM5 domain.
#           Y0 refers to the north-south dimension.  Minimum value is 1.
#   NCOLS:  Number of columns in output MCIP domain (excluding MCIP
#           lateral boundaries).
#   NROWS:  Number of rows in output MCIP domain (excluding MCIP
#           lateral boundaries).
#-----------------------------------------------------------------------

set X0    =  13
set Y0    =  94
set NCOLS =  89
set NROWS = 104

#-----------------------------------------------------------------------
# Set coordinates for cell for diagnostic prints on output domain.
# If coordinate is set to 0, domain center cell will be used.
#-----------------------------------------------------------------------

set LPRT_COL = 0
set LPRT_ROW = 0

#-----------------------------------------------------------------------
# Optional:  Set WRF Lambert conformal reference latitude.
#            (Handy for matching WRF grids to existing MM5 grids.)
#            If not set, MCIP will use average of two true latitudes.
# To "unset" this variable, set the script variable to "-999.0".
# Alternatively, if the script variable is removed here, remove it
# from the setting of the namelist (toward the end of the script).
#-----------------------------------------------------------------------

set WRF_LC_REF_LAT = 36.316

#=======================================================================
#=======================================================================
# Set up and run MCIP.
#   Should not need to change anything below here.
#=======================================================================
#=======================================================================

set PROG = mcip

date

#-----------------------------------------------------------------------
# Make sure directories exist.
#-----------------------------------------------------------------------

if ( ! -d $InMetDir ) then
  echo "No such input directory $InMetDir"
  exit 1
endif

if ( ! -d $OutDir ) then
  echo "No such output directory...will try to create one"
  mkdir -p $OutDir
  if ( $status != 0 ) then
    echo "Failed to make output directory, $OutDir"
    exit 1
  endif
endif

if ( ! -d $ProgDir ) then
  echo "No such program directory $ProgDir"
  exit 1
endif

#-----------------------------------------------------------------------
# Make sure the input files exist.
#-----------------------------------------------------------------------

if ( $IfGeo == "T" ) then
  if ( ! -f $InGeoFile ) then
    echo "No such input file $InGeoFile"
    exit 1
  endif
endif

foreach fil ( $InMetFiles )
  if ( ! -f $fil ) then
    echo "No such input file $fil"
    exit 1
  endif
end

#-----------------------------------------------------------------------
# Make sure the executable exists.
#-----------------------------------------------------------------------

if ( ! -f $ProgDir/${PROG}.exe ) then
  echo "Could not find ${PROG}.exe"
  exit 1
endif

#-----------------------------------------------------------------------
# Create a work directory for this job.
#-----------------------------------------------------------------------

if ( ! -d $WorkDir ) then
  mkdir -p $WorkDir
  if ( $status != 0 ) then
    echo "Failed to make work directory, $WorkDir"
    exit 1
  endif
endif

cd $WorkDir

#-----------------------------------------------------------------------
# Set up script variables for input files.
#-----------------------------------------------------------------------

if ( $IfGeo == "T" ) then
  if ( -f $InGeoFile ) then
    set InGeo = $InGeoFile
  else
    set InGeo = "no_file"
  endif
else
  set InGeo = "no_file"
endif

set FILE_GD  = $OutDir/GRIDDESC

#-----------------------------------------------------------------------
# Create namelist with user definitions.
#-----------------------------------------------------------------------

set MACHTYPE = `uname`
if ( ( $MACHTYPE == "AIX" ) || ( $MACHTYPE == "Darwin" ) ) then
  set Marker = "/"
else
  set Marker = "&END"
endif

cat > $WorkDir/namelist.${PROG} << !
 &FILENAMES
  file_gd    = "$FILE_GD"
  file_mm    = "$InMetFiles[1]",
!

if ( $#InMetFiles > 1 ) then
  @ nn = 2
  while ( $nn <= $#InMetFiles )
    cat >> $WorkDir/namelist.${PROG} << !
               "$InMetFiles[$nn]",
!
    @ nn ++
  end
endif

if ( $IfGeo == "T" ) then
cat >> $WorkDir/namelist.${PROG} << !
  file_geo   = "$InGeo"
!
endif

cat >> $WorkDir/namelist.${PROG} << !
  ioform     =  $IOFORM
 $Marker
 &USERDEFS
  inmetmodel =  $IN_MET_MODEL
  lpv        =  $LPV
  lwout      =  $LWOUT
  luvbout    =  $LUVBOUT
  mcip_start = "$MCIP_START"
  mcip_end   = "$MCIP_END"
  intvl      =  $INTVL
  coordnam   = "$CoordName"
  ctmlays    =  $CTMLAYS
  grdnam     = "$GridName"
  btrim      =  $BTRIM
  lprt_col   =  $LPRT_COL
  lprt_row   =  $LPRT_ROW
  ntimes     =  $N_TIMES
  wrf_lc_ref_lat = $WRF_LC_REF_LAT
 $Marker
 &WINDOWDEFS
  x0         =  $X0
  y0         =  $Y0
  ncolsin    =  $NCOLS
  nrowsin    =  $NROWS
 $Marker
!

#-----------------------------------------------------------------------
# Set links to FORTRAN units.
#-----------------------------------------------------------------------

rm fort.*
if ( -f $FILE_GD ) rm -f $FILE_GD

ln -s $FILE_GD                   fort.4
ln -s $WorkDir/namelist.${PROG}  fort.8

set NUMFIL = 0
foreach fil ( $InMetFiles )
  @ NN = $NUMFIL + 10
  ln -s $fil fort.$NN
  @ NUMFIL ++
end

#-----------------------------------------------------------------------
# Set output file names and other miscellaneous environment variables.
#-----------------------------------------------------------------------

setenv IOAPI_CHECK_HEADERS  T
setenv EXECUTION_ID         $PROG

setenv GRID_BDY_2D          $OutDir/GRIDBDY2D_${APPL}.nc
setenv GRID_CRO_2D          $OutDir/GRIDCRO2D_${APPL}.nc
setenv GRID_DOT_2D          $OutDir/GRIDDOT2D_${APPL}.nc
setenv MET_BDY_3D           $OutDir/METBDY3D_${APPL}.nc
setenv MET_CRO_2D           $OutDir/METCRO2D_${APPL}.nc
setenv MET_CRO_3D           $OutDir/METCRO3D_${APPL}.nc
setenv MET_DOT_3D           $OutDir/METDOT3D_${APPL}.nc
setenv LUFRAC_CRO           $OutDir/LUFRAC_CRO_${APPL}.nc
setenv SOI_CRO              $OutDir/SOI_CRO_${APPL}.nc
setenv MOSAIC_CRO           $OutDir/MOSAIC_CRO_${APPL}.nc

if ( -f $GRID_BDY_2D ) rm -f $GRID_BDY_2D
if ( -f $GRID_CRO_2D ) rm -f $GRID_CRO_2D
if ( -f $GRID_DOT_2D ) rm -f $GRID_DOT_2D
if ( -f $MET_BDY_3D  ) rm -f $MET_BDY_3D
if ( -f $MET_CRO_2D  ) rm -f $MET_CRO_2D
if ( -f $MET_CRO_3D  ) rm -f $MET_CRO_3D
if ( -f $MET_DOT_3D  ) rm -f $MET_DOT_3D
if ( -f $LUFRAC_CRO  ) rm -f $LUFRAC_CRO
if ( -f $SOI_CRO     ) rm -f $SOI_CRO
if ( -f $MOSAIC_CRO  ) rm -f $MOSAIC_CRO

if ( -f $OutDir/mcip.nc      ) rm -f $OutDir/mcip.nc
if ( -f $OutDir/mcip_bdy.nc  ) rm -f $OutDir/mcip_bdy.nc

#-----------------------------------------------------------------------
# Execute MCIP.
#-----------------------------------------------------------------------

$ProgDir/${PROG}.exe

if ( $status == 0 ) then
  rm fort.*
  exit 0
else
  echo "Error running $PROG"
  exit 1
endif
