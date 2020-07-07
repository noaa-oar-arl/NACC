!------------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in           !
!  continuous development by various groups and is based on information        !
!  from these groups: Federal Government employees, contractors working        !
!  within a United States Government contract, and non-Federal sources         !
!  including research institutions.  These groups give the Government          !
!  permission to use, prepare derivative works of, and distribute copies       !
!  of their work in the CMAQ system to the public and to permit others         !
!  to do so.  The United States Environmental Protection Agency                !
!  therefore grants similar permission to use the CMAQ system software,        !
!  but users are requested to provide copies of derivative works or            !
!  products designed to operate in the CMAQ system to the United States        !
!  Government without restrictions as to use by others.  Software              !
!  that is used with the CMAQ system but distributed under the GNU             !
!  General Public License or the GNU Lesser General Public License is          !
!  subject to their copyright restrictions.                                    !
!------------------------------------------------------------------------------!

SUBROUTINE rdfv3 (mcip_now,nn)

!-------------------------------------------------------------------------------
! Name:     Read WRFv2 and WRFv3 (Eulerian Mass Core) Output
! Purpose:  Reads incoming WRFv2 and WRFv3 output files for use in MCIP.
! Notes:    Adapted from S.-B. Kim's get_wrf.F in WCIP.
! Revised:  31 Mar 2005  Original version.  (T. Otte)
!           15 Jul 2005  Modified variable retrievals so that the code will
!                        stop if a variable is not found.  Corrected print
!                        statement for sample output.  (T. Otte)
!           30 Jul 2007  Corrected error in processing incremental precipitation
!                        for first MCIP output time period when the first WRF
!                        WRF output time is not used by MCIP.  Added reads for
!                        fractional land use, leaf area index, aerodynamic and
!                        stomatal resistances, inverse Monin-Obukhov length, and
!                        soil moisture, temperature, and type, if those fields
!                        are available.  Removed read for emissivity.  Allowed
!                        for roughness length to be filled from a lookup table
!                        if it is not available in WRF output.  (T. Otte)
!           14 May 2008  Read TSLB (layer 1) if TSK is unavailable.  Change
!                        static data tables for roughness length to allow for
!                        up to 33 categories of USGS, and added error-checking
!                        when ZNT is set from the lookup table.  Corrected
!                        season assignment for lookup table for the Southern
!                        Hemisphere.  Check LAI, RA, and RSTOM to ensure that
!                        there are non-zero values in the fields, if they
!                        exist.  If the values of RA and/or RSTOM are all 0.0,
!                        reset IFRESIST flag so that they will be calculated
!                        later.  If LAI is in output but is 0.0, set LAI to
!                        realistic values for NOAH LSM.  Added 2-m mixing ratio
!                        (Q2) and turbulent kinetic energy (TKE), if available.
!                        Changed read on vegetation fraction to preferentially
!                        use VEGF_PX rather than VEGFRA for Pleim-Xiu land-
!                        surface model.  Changed algorithm to find "valid"
!                        data to require time difference to be < TTOL rather
!                        than <= TTOL.  Added urban fraction (FRC_URB),
!                        urban roughness length (Z0C_URB2D), and urban Monin-
!                        Obukhov length (XXXC_URB) for MET_UCMCALL=1.  Added
!                        error checking to ensure that WRF files used in this
!                        MCIP run are from the same simulation so that
!                        incremental precipitation totals in RN and RC are
!                        processed correctly.  (T. Otte)
!           29 Oct 2009  Cleaned up file opening and logging for WRF I/O API,
!                        particularly when the WRF headers of new files are
!                        checked, to prevent condition with too many files open
!                        for long simulations.  Changed MET_UCMCALL to
!                        MET_URBAN_PHYS, and allowed for variable to be set to
!                        be greater than 1.  Capture potential temperature
!                        (THETA) and Coriolis (CORIOLIS) when potential
!                        vorticity is needed.  Changed method of computing
!                        latitude, longitude, and map-scale factor arrays to
!                        be more general; removed subroutine GRIDGEOMETRY.
!                        Added default roughness length values for NCLD-MODIS,
!                        SiB, and MODIS-NOAH.  Increased MAX_TIMES to 1000 to
!                        enable processing of longer data sets.  Removed
!                        DUM2D_D.  Added latitude, longitude, and map-scale
!                        factors on U and V faces.  Allow output from WRF
!                        Preprocessing System (WPS) routine, GEOGRID, to
!                        provide fractional land use output if it is unavailable
!                        in WRF output.  Removed Z0C_URB2D.  Corrected units
!                        for U10 and V10 in log file.  Changed error condition
!                        to warning condition if LAI is set to zero on input
!                        and LSM other than NOAH was used.  Changed reads of
!                        fractional land use and roughness length so that they
!                        are only performed if those fields are known to exist.
!                        Changed real-number comparisons of maximum values from
!                        "equivalences" to "less than tolerances".  (T. Otte)
!           12 Feb 2010  Removed unused variables COMM and SYSDEP_INFO, and
!                        removed unused format 9600.  (T. Otte)
!           18 Mar 2010  Added CDFID as an input argument for subroutine
!                        CHKWRFHDR.  Changed all calls to netCDF routines to use
!                        the Fortran interface rather than the C interface.
!                        Changed input arguments for routines in WRF_NETCDF
!                        from FILENAME to CDFID to minimize I/O.  (T. Otte)
!           23 Dec 2010  Improved support for long MCIP runs from long WRF
!                        runs by increasing MAX_TIMES to 9999.  Also added
!                        missing "close" command for incoming WRF files.
!                        Added sea ice.  Added support for precipitation
!                        tipping bucket option in WRF.  Changed latitude and
!                        longitude calculations for polar stereograhic
!                        projection to interpolations.  (T. Otte)
!           31 Aug 2011  Changed name of module FILE to FILES to avoid conflict
!                        with F90 protected intrinsic.  Updated netCDF commands
!                        to F90, and improved error handling.  Replaced calls
!                        to GET_TIMES_CDF with explicit netCDF functions.
!                        Changed F77 character declarations to F90 standard.
!                        Changed DATA statements to parameters.  Changed
!                        arguments to 19-character elements for GETH_IDTS.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           21 Nov 2011  Corrected error in tipping bucket precipitation
!                        calculation.  (T. Otte)
!           21 Aug 2012  Added MET_PCP_INCR to accommodate WRFv3.3 option to
!                        output incremental precipitation.  (T. Otte)
!           11 Sep 2012  Added handling for 40-category 2006 NLCD-MODIS land
!                        use classification as "NLCD40".  Added alternate name
!                        for 50-category 2001 NLCD-MODIS land use classification
!                        as "NLCD50".  Changed SFZ0NLCSUM and SFZ0NLCWIN to
!                        SFZ0NLCD50SUM and SFZ0NLCD50WIN.  Added analogous
!                        arrays for NLCD40.  Updated values for roughness length
!                        for both NLCD-MODIS classifications using tables from
!                        WRFv3.4 module_sf_pxlsm_data.F.  Added read of array
!                        LANDMASK to be used with runs that used Pleim-Xiu LSM.
!                        (T. Otte)
!           26 Nov 2014  Corrected formatting error in error-handling F9400,
!                        and corrected four locations of this routine that now
!                        reference F9400 incorrectly.  Removed requirement to
!                        have FRC_URB available when urban canopy model is used
!                        in WRF.  (T. Spero)
!           10 Apr 2015  If 3D resolved cloud fraction is in the WRF output,
!                        collect that field to pass through to output.
!                        (T. Spero)
!           25 Aug 2015  Changed latent heat flux from QFX to LH.  Fill THETA
!                        and moisture flux (QFX) for IFMOLACM.  If Pleim-Xiu
!                        land-surface model is used, realign soil categories
!                        to be consistent with WRF documentation.  (T. Spero)
!           08 Sep 2015  Commented out realignment of soil categories for
!                        Pleim-Xiu land-surface model because CMAQ cannot
!                        handle this yet.  (T. Spero)
!           17 Sep 2015  Changed IFMOLACM to IFMOLPX.  (T. Spero)
!           30 Oct 2015  Changed WRITE statements for printing sampled data to
!                        log file to eliminate warning messages.  (T. Spero)
!           22 Nov 2016  Changed urban model variable FRC_URB to FRC_URB2D to
!                        be consistent with its use in WRF.  (T. Spero)
!           21 Apr 2017  Updated SFZ0 for MODIS so that category 21 is "Lake".
!                        (T. Spero)
!           16 Mar 2018  Corrected the settings for II and JJ in the loop for
!                        calculating dot-point latitude and longitude for polar
!                        stereographic WRF projections.  Moved TTOL to
!                        MCIPPARM_MOD, and changed its local name to TTOL_SEC.
!                        Corrected error in print statement for WRF variable
!                        CLDFRA.  Added SNOWH to output.  Created a minimum
!                        value for rainfall in order to avoid underflow
!                        condition.  Corrected minor error in array mapping in
!                        rain buckets in unused column and row.  Added
!                        LUFRAC2, MOSCATIDX, LAI_MOS, RA_MOS, RS_MOS, TSK_MOS,
!                        ZNT_MOS, and DUM3D_M to support NOAH Mosaic land-
!                        surface model.  Added DZS to capture soil layers, and
!                        added 3D soil arrays, SOIT3D and SOIM3D.  Added
!                        WSPDSFC and XLAIDYN for Noah.  (T. Spero)
!           27 Jun 2018  Changed name of module with netCDF IO to broaden its
!                        usage.  Removed local aliases for dimensions of input
!                        meteorological fields.  (T. Spero)
!           14 Sep 2018  Changed condition to enable hybrid vertical coordinate
!                        from WRF.  Removed support for MM5v3 input.  (T. Spero)
!           16 Oct 2018  Corrected error in computing precipitation amounts when
!                        the tipping bucket is used and less than 0.5 mm of
!                        precipitation accumulated during the same hour that the
!                        bucket tips; corrects erroneous precipitation spikes.
!                        Corrected error in array mapping for precipitation on
!                        initial time step in the outermost row and column
!                        (dummy cells not used by CMAQ).
!                        (C. Nolte and T. Spero)
!           23 Nov 2018  Changed local usages of NX, NY, and NZ to MET_NX,
!                        MET_NY, and MET_NZ to avoid confusion with generic
!                        usages of those variables for global dimensions in
!                        netCDF output.  (T. Spero)
!           18 Jun 2019  Added new surface variables with PX LSM that can
!                        improve dust simulation in CCTM.  Added optional
!                        variables from KF convective scheme with radiative
!                        feedbacks.  (T. Spero)
!           24 Feb 2020  Adapted for FV3GFSv16 at NOAA-ARL (P. C. Campbell)
!           24 Feb 2020  Added horiz LCC interpolation and wind rotation 
!                        Y. Tang and P. C. Campbell)
!-------------------------------------------------------------------------------

  USE date_pack
  USE files
  USE metinfo
  USE metvars
  USE coord
  USE mcipparm
  USE netcdf_io
  USE netcdf
  USE m3utilio

  IMPLICIT NONE

  INTEGER, SAVE                     :: cdfid, cdfid2
  INTEGER                           :: cdfidg
  INTEGER                           :: dimid
  INTEGER                           :: dimids     ( nf90_max_var_dims )
  REAL,    SAVE,      ALLOCATABLE   :: dum1d      ( : )
  REAL,    SAVE,      ALLOCATABLE   :: dum2d      ( : , : )
  INTEGER, SAVE,      ALLOCATABLE   :: dum2d_i    ( : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum2d_u    ( : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum2d_v    ( : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_l    ( : , : , : )
  INTEGER, SAVE,      ALLOCATABLE   :: dum3d_li   ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_m    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_p    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_s    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_t    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_v    ( : , : , : )
  REAL,    SAVE,      ALLOCATABLE   :: dum3d_w    ( : , : , : )
  CHARACTER(LEN=19)                 :: endseas
  LOGICAL, SAVE                     :: first      = .TRUE.
  CHARACTER(LEN=256)                :: fl
  CHARACTER(LEN=256)                :: flg
  CHARACTER(LEN=32)                 :: date_init
  LOGICAL                           :: gotfaces   = .TRUE.
  LOGICAL                           :: gotseaice
  LOGICAL                           :: gotznt
  INTEGER                           :: i,varid
  INTEGER                           :: id_data
  INTEGER                           :: idts_end
  INTEGER                           :: idts_start
  INTEGER                           :: idtsec
  LOGICAL                           :: iffl
  CHARACTER(LEN=64)                 :: ifmt1
  CHARACTER(LEN=64)                 :: ifmt1a
  CHARACTER(LEN=64)                 :: ifmt2
  CHARACTER(LEN=64)                 :: ifmt3
  CHARACTER(LEN=64)                 :: ifmt4
  CHARACTER(LEN=64)                 :: ifmt5
  INTEGER                           :: ii
  INTEGER                           :: im1
  INTEGER                           :: it
  INTEGER, SAVE                     :: it_start
  INTEGER                           :: itm1
  INTEGER                           :: j
  INTEGER                           :: jj
  INTEGER                           :: jm1
  INTEGER                           :: k,kk
  INTEGER                           :: k1
  INTEGER                           :: k2
  INTEGER                           :: lent
  INTEGER, intent(in)               :: nn
  REAL,               EXTERNAL      :: mapfac_lam
  REAL,               EXTERNAL      :: mapfac_merc
  REAL,               EXTERNAL      :: mapfac_ps
  CHARACTER(LEN=24),  INTENT(INOUT) :: mcip_now
  CHARACTER(LEN=24)                 :: mcip_previous,mcip_rd,mcip_next
  INTEGER                           :: m1count    = 1
  INTEGER, SAVE                     :: mmcount    = 1
  INTEGER, SAVE                     :: n_times
  LOGICAL, SAVE                     :: newfile    = .TRUE.
  LOGICAL                           :: newfilem1  = .TRUE.
  INTEGER                           :: nxm
  INTEGER                           :: nym
  INTEGER                           :: nzp
  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'RDFV3'
  INTEGER                           :: rcode, rcode2
  REAL,               PARAMETER     :: rdovcp     = 2.0 / 7.0
  REAL,               PARAMETER     :: smallnum   = 1.0e-7
  CHARACTER(LEN=19)                 :: startseas
  CHARACTER(LEN=2)                  :: str1
  CHARACTER(LEN=2)                  :: str2
  CHARACTER(LEN=3)                  :: str3
!  CHARACTER(LEN=19),SAVE,ALLOCATABLE:: times      ( : )
  REAL,SAVE,ALLOCATABLE             :: times      ( : ), atmp(:,:), utmp(:,:)
  REAL                              :: xoff
  REAL                              :: xxin
  REAL                              :: yoff
  REAL                              :: yyin
  double precision  :: rdtime
  ! Define roughness length as functions of land use and season in case
  ! it is not available in WRF output.

  REAL, PARAMETER :: sfz0oldsum ( 13 ) = &  ! summer [cm]
    (/ 50.0,  15.0,  12.0,  50.0,  50.0,  40.0,  0.01, 20.0,   &
       10.0,  10.0,   5.0,  50.0,  15.0 /)

  REAL, PARAMETER :: sfz0oldwin ( 13 ) = &  ! winter [cm]
    (/ 50.0,   5.0,  10.0,  50.0,  50.0,  40.0,  0.01, 20.0,   &
       10.0,  10.0,   5.0,  50.0,  15.0 /)

  REAL, PARAMETER :: sfz0modsum ( 33 ) = &  ! summer [cm]
    (/ 50.0,  50.0,  50.0,  50.0,  50.0,   5.0,  6.0,   5.0,   &
       15.0,  12.0,  30.0,  15.0,  80.0,  14.0,  0.1,   1.0,   &
        0.01, 30.0,  15.0,  10.0,   0.01, 80.0,  80.0,  80.0,  &
       80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0 /)

  REAL, PARAMETER :: sfz0modwin ( 33 ) = &  ! winter [cm]
    (/ 50.0,  50.0,  50.0,  50.0,  20.0,   1.0,   1.0,   1.0,  &
       15.0,  50.0,  30.0,   5.0,  80.0,   5.0,   0.1,   1.0,  &
        0.01, 10.0,  30.0,  15.0,   0.01, 80.0,  80.0,  80.0,  &
       80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0 /)

  REAL, PARAMETER :: sfz0nlcd50sum ( 50 ) = &  ! summer [cm]
    (/  0.1,   1.2,  30.0,  40.0,  60.0, 100.0,   5.0,   5.0,  &
      100.0, 100.0, 100.0,  10.0,  15.0,   7.0,   7.0,   5.0,  &
        5.0,   5.0,   7.0,  10.0,  55.0,  80.0,  30.0,  60.0,  &
       30.0,  11.0,  11.0,  11.0,   5.0,   5.0,   0.1, 100.0,  &
       90.0, 100.0, 100.0, 100.0,  15.0,  15.0,  25.0,  15.0,  &
        7.0,  20.0,  10.0,  80.0,  30.0,   1.2,   5.0,   0.1,  &
        0.1,   0.1 /)

  REAL, PARAMETER :: sfz0nlcd50win ( 50 ) = &  ! winter [cm]
    (/  0.1,   1.2,  30.0,  40.0,  60.0, 100.0,   5.0,   5.0,  &
      100.0, 100.0, 100.0,  10.0,  15.0,   7.0,   7.0,   5.0,  &
        5.0,   5.0,   7.0,  10.0,  55.0,  80.0,  30.0,  60.0,  &
       30.0,  11.0,  11.0,  11.0,   5.0,   5.0,   0.1, 100.0,  &
       90.0, 100.0, 100.0, 100.0,  15.0,  15.0,  25.0,  15.0,  &
        7.0,  20.0,  10.0,  80.0,  30.0,   1.2,   5.0,   0.1,  &
        0.1,   0.1 /)

  REAL, PARAMETER :: sfz0nlcd40sum ( 40 ) = &  ! summer [cm]
    (/100.0,  90.0, 100.0, 100.0, 100.0,  30.0,  15.0,  25.0,  &
       15.0,   7.0,  20.0,  10.0,  80.0,  30.0,   1.2,   5.0,  &
        0.1,   0.1,   0.1,   0.1,   0.1,   1.2,  30.0,  40.0,  &
       60.0, 100.0,   5.0, 100.0, 100.0, 100.0,  10.0,  15.0,  &
        7.0,   7.0,   5.0,   5.0,   7.0,  10.0,  55.0,  11.0 /)

  REAL, PARAMETER :: sfz0nlcd40win ( 40 ) = &  ! winter [cm]
    (/100.0,  90.0, 100.0, 100.0, 100.0,  30.0,  15.0,  25.0,  &
       15.0,   7.0,  20.0,  10.0,  80.0,  30.0,   1.2,   5.0,  &
        0.1,   0.1,   0.1,   0.1,   0.1,   1.2,  30.0,  40.0,  &
       60.0, 100.0,   5.0, 100.0, 100.0, 100.0,  10.0,  15.0,  &
        7.0,   7.0,   5.0,   5.0,   7.0,  10.0,  55.0,  11.0 /)

  REAL, PARAMETER :: sfz0sibsum ( 16 ) = &  ! summer [cm]
    (/ 50.0,  50.0,  40.0,  50.0,  50.0,  15.0,  12.0,  12.0,  &
       12.0,  10.0,  10.0,  15.0,  20.0,  12.0,   0.01,  5.0 /)

  REAL, PARAMETER :: sfz0sibwin ( 16 ) = &  ! winter [cm]
    (/ 50.0,  50.0,  40.0,  50.0,  50.0,  15.0,  10.0,  10.0,  &
       10.0,  10.0,  10.0,   5.0,  20.0,  10.0,   0.01,  5.0 /)

  REAL, PARAMETER :: sfz0usgssum ( 33 ) = &  ! summer [cm]
    (/ 80.0,  15.0,  10.0,  15.0,  14.0,  20.0,  12.0,   5.0,  &
        6.0,  15.0,  50.0,  50.0,  50.0,  50.0,  50.0,   0.01, &
       20.0,  40.0,   1.0,  10.0,  30.0,  15.0,  10.0,   5.0,  &
        1.0,  15.0,   1.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0 /)

  REAL, PARAMETER :: sfz0usgswin ( 33 ) = &  ! winter [cm]
    (/ 80.0,   5.0,   2.0,   5.0,   5.0,  20.0,  10.0,   1.0,  &
        1.0,  15.0,  50.0,  50.0,  50.0,  50.0,  20.0,   0.01, &
       20.0,  40.0,   1.0,  10.0,  30.0,  15.0,   5.0,   5.0,  &
        1.0,  15.0,   1.0,  80.0,  80.0,  80.0,  80.0,  80.0,  80.0 /)

!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f6000 = "(1x, a, 1x, f12.4, 2x, a)"
  CHARACTER(LEN=256), PARAMETER :: f6100 = "(1x, a, 1x, i12,   2x, a)"

  CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   LOOKING FOR INPUT MET AT TIME ', a, &
    & /, 1x, '***   NO MORE INPUT FV3 FILES', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9200 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   LOOKING FOR INPUT MET AT TIME ', a, &
    & /, 1x, '***   INPUT FILE NUMBER ', i3, ' IS BLANK', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   LOOKING FOR INPUT MET AT TIME ', a, &
    & /, 1x, '***   COULD NOT FIND FILE ', a, &
    & /, 1x, '***   FILE MAY NOT EXIST', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   RCODE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9410 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING NCF ID FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9420 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR INQUIRING ABOUT VAR IN FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9430 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING DIMS FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9500 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNKNOWN LAND USE CLASSIFICATION SYSTEM', &
    & /, 1x, '***   LAND USE SOURCE = ', a, &
    & /, 1x, '***   HIGHEST INDEX FOUND = ', i4, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9700 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNABLE TO SET ZNT FROM LOOKUP TABLE', &
    & /, 1x, '***   LAND USE SOURCE = ', a, &
    & /, 1x, '***   NUMBER OF LAND USE CATEGORIES = ', i3, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9800 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNABLE TO BLEND ', a, ' FOR UCM', &
    & /, 1x, '***   UNKNOWN LAND USE SOURCE', &
    & /, 1x, '***   LAND USE SOURCE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9900 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING FV3 NETCDF FILE', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9950 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING FV3 NETCDF FILE', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9975 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   DID NOT FIND ARRAY ', a, &
    & /, 1x, '***   WILL DEFINE FROM OTHER FIELDS LATER', &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Interfaces for FV3GFS getxyindex, horizontal interpolation, and wind rotation
! 
!-------------------------------------------------------------------------------

  INTERFACE

    SUBROUTINE getxyindex (xlat,xlon,xi,yj,tlat,tlon,ix,jy)
    IMPLICIT NONE
    REAL,               INTENT(IN)    :: xlat
    REAL,               INTENT(IN)    :: xlon
    REAL,               INTENT(IN)    :: tlat ( : )
    REAL,               INTENT(IN)    :: tlon ( : )
    INTEGER,            INTENT(IN)    :: ix, jy
    REAL,               INTENT(OUT)   :: xi
    REAL,               INTENT(OUT)   :: yj

    END SUBROUTINE getxyindex

  END INTERFACE

  INTERFACE

    SUBROUTINE myinterp (ain,met_nx,met_ny,aout,xindex,yindex,iout,jout,iflag)
    IMPLICIT NONE
    REAL,               INTENT(IN)    :: ain        ( : , : )
    REAL,               INTENT(IN)    :: xindex     ( : , : )
    REAL,               INTENT(IN)    :: yindex     ( : , : )
    INTEGER,            INTENT(IN)    :: met_nx, met_ny
    INTEGER,            INTENT(IN)    :: iout, jout
    INTEGER,            INTENT(IN)    :: iflag
    REAL,               INTENT(OUT)   :: aout       ( : , : )

    END SUBROUTINE myinterp

  END INTERFACE

  INTERFACE

    SUBROUTINE windrotation (ua,va,londot,ix,jy,reflon,reflat)
    IMPLICIT NONE
    REAL,               INTENT(INOUT) :: ua         ( : , : )
    REAL,               INTENT(INOUT) :: va         ( : , : )
    REAL,               INTENT(IN)    :: londot     ( : , : )
    INTEGER,            INTENT(IN)    :: ix, jy
    REAL,               INTENT(IN)    :: reflon, reflat

    END SUBROUTINE windrotation

  END INTERFACE

!-------------------------------------------------------------------------------
! Define additional staggered grid dimensions. (***No staggered FV3 dimensions,e.g., nxm=met_nx***)
!-------------------------------------------------------------------------------

  nxm = met_nx 
  nym = met_ny 
  nzp = met_nz 

!-------------------------------------------------------------------------------
! Set up print statements.
!-------------------------------------------------------------------------------

  k1 = met_nz / 5
  k2 = MOD(met_nz, 5)

  WRITE ( str1, '(i2)' ) k1 - 1
  WRITE ( str2, '(i2)' ) k2

  IF ( (k1 - 1) > 0 ) THEN
    IF ( k2 > 0 ) THEN
      ifmt1 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)),/,10x," &
        & // str2 // "(1x,f12.4))"
    ELSE
      ifmt1 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)))"
    ENDIF
  ELSE
    IF ( k2 > 0 ) THEN
      ifmt1 = "(/,1x,a,5(1x,f12.4),/,11x," // str2 // "(1x,f12.4))"
    ELSE
      ifmt1 = "(/,1x,a,5(1x,f12.4))"
    ENDIF
  ENDIF

  k1 = (nzp) / 5
  k2 = MOD(nzp, 5)

  WRITE ( str1, '(i2)' ) k1 - 1
  WRITE ( str2, '(i2)' ) k2

  IF ( (k1 - 1) > 0 ) THEN
    IF ( k2 > 0 ) THEN
      ifmt1a = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)),/,10x," &
        & // str2 // "(1x,f12.4))"
    ELSE
      ifmt1a = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)))"
    ENDIF
  ELSE
    IF ( k2 > 0 ) THEN
      ifmt1a = "(/,1x,a,5(1x,f12.4),/,10x," // str2 // "(1x,f12.4))"
    ELSE
      ifmt1a = "(/,1x,a,5(1x,f12.4))"
    ENDIF
  ENDIF

  k1 = nummetlu / 5
  k2 = MOD(nummetlu, 5)

  WRITE ( str1, '(i2)' ) k1 - 1
  WRITE ( str2, '(i2)' ) k2

  IF ( (k1 - 1) > 0 ) THEN
    IF ( k2 > 0 ) THEN
      ifmt2 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)),/,10x," &
        & // str2 // "(1x,f12.4))"
      ifmt3 = "(/,1x,a,5(i12,1x)," // str1 // "(/,10x,5(1x,i12)),/,10x," &
        & // str2 // "(1x,i12))"
    ELSE
      ifmt2 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)))"
      ifmt3 = "(/,1x,a,5(i12,1x)," // str1 // "(/,10x,5(1x,i12)))"
    ENDIF
  ELSE
    IF ( k2 > 0 ) THEN
      ifmt2 = "(/,1x,a,5(1x,f12.4),/,10x," // str2 // "(1x,f12.4))"
      ifmt3 = "(/,1x,a,5(i12,1x),/,10x," // str2 // "(1x,i12))"
    ELSE
      ifmt2 = "(/,1x,a,5(1x,f12.4))"
      ifmt3 = "(/,1x,a,5(i12,1x))"
    ENDIF
  ENDIF

  k1 = nummosaic / 5
  k2 = MOD(nummosaic, 5)

  WRITE ( str1, '(i2)' ) k1 - 1
  WRITE ( str2, '(i2)' ) k2

  IF ( (k1 - 1) > 0 ) THEN
    IF ( k2 > 0 ) THEN
      ifmt4 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)),/,10x," &
        & // str2 // "(1x,f12.4))"
    ELSE
      ifmt4 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)))"
    ENDIF
  ELSE
    IF ( k2 > 0 ) THEN
      ifmt4 = "(/,1x,a,5(1x,f12.4),/,10x," // str2 // "(1x,f12.4))"
    ELSE
      ifmt4 = "(/,1x,a,5(1x,f12.4))"
    ENDIF
  ENDIF

  k1 = met_ns / 5
  k2 = MOD(met_ns, 5)

  WRITE ( str1, '(i2)' ) k1 - 1
  WRITE ( str2, '(i2)' ) k2

  IF ( (k1 - 1) > 0 ) THEN
    IF ( k2 > 0 ) THEN
      ifmt5 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)),/,10x," &
        & // str2 // "(1x,f12.4))"
    ELSE
      ifmt5 = "(/,1x,a,5(1x,f12.4)," // str1 // "(/,10x,5(1x,f12.4)))"
    ENDIF
  ELSE
    IF ( k2 > 0 ) THEN
      ifmt5 = "(/,1x,a,5(1x,f12.4),/,11x," // str2 // "(1x,f12.4))"
    ELSE
      ifmt5 = "(/,1x,a,5(1x,f12.4))"
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Allocate necessary variables.
!-------------------------------------------------------------------------------
  IF ( .NOT. ALLOCATED ( dum3d ) )  &
    ALLOCATE ( dum3d (met_nx, met_ny, met_nz ) ) ! 3D, N-S flux pts, full lvls


  IF ( .NOT. ALLOCATED ( dum2d   ) )  &
    ALLOCATE ( dum2d   (met_nx, met_ny)      )        ! 2D, cross points
  IF ( .NOT. ALLOCATED ( dum2d_i ) )  &
    ALLOCATE ( dum2d_i (met_nx, met_ny)      )        ! 2D integer, cross points
  IF ( .NOT. ALLOCATED ( dum2d_u ) )  &
    ALLOCATE ( dum2d_u (met_nx, met_ny ) )         ! 2D, E-W flux pts
  IF ( .NOT. ALLOCATED ( dum2d_v ) )  &
    ALLOCATE ( dum2d_v (met_nx, met_ny  ) )        ! 2D, N-S flux pts
  IF ( .NOT. ALLOCATED ( dum3d_l ) )  &
    ALLOCATE ( dum3d_l (met_nx, met_ny, nummetlu ) )  ! 3D, cross points, lu
  IF ( .NOT. ALLOCATED ( dum3d_li ) ) &
    ALLOCATE ( dum3d_li (met_nx, met_ny, nummetlu ) ) ! 3D, cross points, lu int
  IF ( .NOT. ALLOCATED ( dum3d_m ) )  &
    ALLOCATE ( dum3d_m (met_nx, met_ny, nummosaic) )  ! 3D, cross pts in mosaic cat
  IF ( .NOT. ALLOCATED ( dum3d_s ) )  &
    ALLOCATE ( dum3d_s (met_nx, met_ny, met_ns ) )    ! 3D, cross points, soil lvls
!  IF ( .NOT. ALLOCATED ( dum3d_w ) )  &
!    ALLOCATE ( dum3d_w (met_nx, met_ny, met_nz ) )    ! 3D, cross points, full lvls
!  IF ( .NOT. ALLOCATED ( dum3d_t ) )  &
!    ALLOCATE ( dum3d_t (met_nx, met_ny, met_nz ) )    ! 3D, cross points, full lvls
  if(.not.allocated(atmp)) allocate(atmp(ncols_x,nrows_x))
  if(.not.allocated(utmp)) allocate(utmp(ncols_x+1,nrows_x+1))


!-------------------------------------------------------------------------------
! If not processing the first output time of the WRF run (and if not using the
! incremental precipitation option available in WRFv3.2+), retrieve accumulated
! precipitation totals from time increment before first MCIP step so that
! first incremental precipitation "rates" can be computed.  This step ensures
! that the "hold" values for convective and non-convective precipitation are
! correctly set with last accumulated total.
!-------------------------------------------------------------------------------

  gotseaice = .FALSE.
!FV3 does not contain the lat/lon of cell faces
  gotfaces = .FALSE.

  if (first) then
   allocate(xindex(ncols_x,nrows_x), yindex(ncols_x,nrows_x), xuindex(ncols_x+1,nrows_x+1), &
    yuindex(ncols_x+1,nrows_x+1),xvindex(ncols_x+1,nrows_x+1),yvindex(ncols_x+1,nrows_x+1), &
    xdindex(ncols_x+1,nrows_x+1),ydindex(ncols_x+1,nrows_x+1))


    ! Compute distance from origin (at reflat, standlon) to domain center, and
    ! store in MET_XXCTR and MET_YYCTR.  Then calculate latitude, longitude,
    ! and map-scale factors using offset distance of given grid point from
    ! center of domain.

      if (gdtyp_gd.eq.lamgrd3) then 

        xoff = -2.0  ! extend one more point from dot-point center value to cover boundary
        yoff = -2.0  

        met_tru1=sngl(p_alp_gd)
	met_tru2=sngl(p_bet_gd)
	met_proj_clon=sngl(p_gam_gd)
        met_ref_lat=sngl(ycent_gd)
	
	met_xxctr=sngl(xorig_gd)  ! in meter
	met_yyctr=sngl(yorig_gd)

	
        DO j = 1, nrows_x+1
          DO i = 1, ncols_x+1

            xxin = met_xxctr + (FLOAT(i) + xoff) * met_resoln
            yyin = met_yyctr + (FLOAT(j) + yoff) * met_resoln

            CALL xy2ll_lam (xxin, yyin, met_tru1, met_tru2, met_proj_clon,  &
                            met_ref_lat, latdot(i,j), londot(i,j))

            mapdot(i,j) = mapfac_lam (latdot(i,j), met_tru1, met_tru2)
	    call getxyindex(latdot(i,j),londot(i,j),xdindex(i,j),ydindex(i,j),fv3lat,fv3lon,met_nx,met_ny)
	    
          ENDDO
        ENDDO

        xoff=-1.5
	yoff=-1.5
        DO j = 1, nrows_x
          DO i = 1, ncols_x

            xxin = met_xxctr + (FLOAT(i) + xoff) * met_resoln
            yyin = met_yyctr + (FLOAT(j) + yoff) * met_resoln

            CALL xy2ll_lam (xxin, yyin, met_tru1, met_tru2, met_proj_clon,  &
                            met_ref_lat, latcrs(i,j), loncrs(i,j))

            mapcrs(i,j) = mapfac_lam (latcrs(i,j), met_tru1, met_tru2)

            call getxyindex(latcrs(i,j),loncrs(i,j),xindex(i,j),yindex(i,j),fv3lat,fv3lon,met_nx,met_ny)
	    
          ENDDO
        ENDDO


        IF ( .NOT. gotfaces ) THEN  ! get lat, lon, map-scale factor on faces

          xoff = -2.0  ! U-face: no offset in X from dot-point center value
          yoff = -1.5  ! U-face: 0.5-cell offset in Y from dot-point center value

          DO j = 1, nrows_x+1  ! use all Y to fill array; last row outside domain
            DO i = 1, ncols_x+1

              xxin = met_xxctr + (FLOAT(i) + xoff) * met_resoln
              yyin = met_yyctr + (FLOAT(j) + yoff) * met_resoln

              CALL xy2ll_lam (xxin, yyin, met_tru1, met_tru2, met_proj_clon,  &
                              met_ref_lat, latu(i,j), lonu(i,j))

              mapu(i,j) = mapfac_lam (latu(i,j), met_tru1, met_tru2)
 	     
	      call getxyindex(latu(i,j),lonu(i,j),xuindex(i,j),yuindex(i,j),fv3lat,fv3lon,met_nx,met_ny)

            ENDDO
          ENDDO

          xoff = -1.5  ! V-face: 0.5-cell offset in X from dot-point center value
          yoff = -2.0  ! V-face: no offset in Y from dot-point center value

          DO j = 1, nrows_x+1
            DO i = 1, nrows_x+1  ! use all X to fill array; last col outside domain

              xxin = met_xxctr + (FLOAT(i) + xoff) * met_resoln
              yyin = met_yyctr + (FLOAT(j) + yoff) * met_resoln

              CALL xy2ll_lam (xxin, yyin, met_tru1, met_tru2, met_proj_clon,  &
                              met_ref_lat, latv(i,j), lonv(i,j))

              mapv(i,j) = mapfac_lam (latv(i,j), met_tru1, met_tru2)

	      call getxyindex(latv(i,j),lonv(i,j),xvindex(i,j),yvindex(i,j),fv3lat,fv3lon,met_nx,met_ny)


            ENDDO
          ENDDO

        ENDIF

    endif

  endif

!-------------------------------------------------------------------------------------
! Open FV3GFS files and check headers
!-------------------------------------------------------------------------------------


!open 3d atm file
  write(str3,'(i3.3)')nn-1
  rcode = nf90_open (trim(file_mm(1))//str3//trim(file_mm(2)), nf90_nowrite,cdfid)

  IF ( rcode /= nf90_noerr ) THEN
   print*,'error open ATM file', nn,str3,trim(file_mm(1))//str3//trim(file_mm(2))
   call graceful_stop (pname)
  endif 

!open 2d sfc file
  rcode2 = nf90_open (trim(file_sfc(1))//str3//trim(file_sfc(2)), nf90_nowrite,cdfid2)

  IF ( rcode2 /= nf90_noerr ) THEN
   print*,'error open SFC file', nn,str3,trim(file_sfc(1))//str3//trim(file_sfc(2))
   call graceful_stop (pname)
  endif 

!dimension check
  rcode = nf90_get_att (cdfid, nf90_global, 'im', ii)
  IF ( rcode /= nf90_noerr ) THEN
   write(*,*)'error get im ATM file ',str3
   call graceful_stop (pname)
  else 
   if(ii.ne.met_nx) then
    print*,'inconsistent ATM x-dimension ',ii,met_nx,str3
    CALL graceful_stop (pname)
   endif
  endif
  
  rcode2 = nf90_get_att (cdfid2, nf90_global, 'im', ii)
  IF ( rcode2 /= nf90_noerr ) THEN
   write(*,*)'error SFC file ',str3
   call graceful_stop (pname)
  else 
   if(ii.ne.met_nx) then
    print*,'inconsistent SFC x-dimension ',ii,met_nx,str3
    CALL graceful_stop (pname)
   endif
  endif
  rcode = nf90_get_att(cdfid, nf90_global, 'jm', jj)
  IF ( rcode /= nf90_noerr ) THEN
   write(*,*)'error get jm in ATM file ',str3
   call graceful_stop (pname)
  else 
   if(jj.ne.met_ny) then
    print*,'inconsistent ATM y-dimension ',jj,met_ny,str3
    CALL graceful_stop (pname)
   endif
  endif
  
  rcode2 = nf90_get_att(cdfid2, nf90_global, 'jm', jj)
  IF ( rcode2 /= nf90_noerr ) THEN
   write(*,*)'error get jm in SFC file ',str3
   call graceful_stop (pname)
  else 
   if(jj.ne.met_ny) then
    print*,'inconsistent SFC y-dimension ',jj,met_ny,str3
    CALL graceful_stop (pname)
   endif
  endif

! Extract FV3 start date and time information.
!-------------------------------------------------------------------------------
  rcode = nf90_inq_varid (cdfid, 'time', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'time',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF 

  rcode = nf90_get_att (cdfid, varid, 'units', date_init)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SIMULATION_START_DATE',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF  
  mcip_rd = date_init(13:31) // '.0000'
  mcip_rd(11:11)='-'
  rcode=nf90_get_var(cdfid,varid,rdtime)
  IF ( rcode /= nf90_noerr ) THEN
   write(*,*)'error getting time in ATM file',str3
   CALL graceful_stop (pname)
  ENDIF
  if(nn-1.ne.int(rdtime)) then
    print*,'time inconsistent ',nn,rdtime
    stop
  endif
      

  CALL geth_newdate (mcip_next, mcip_rd, int(rdtime)*intvl*60)
  if(nn.eq.1.and.mcip_next.ne.mcip_now) then
   write(*,*)'time mismatch in ATM file ',mcip_now,mcip_next,mcip_rd,date_init,rdtime
   CALL graceful_stop (pname)
  ENDIF
  mcip_now=mcip_next  ! update mcip_now for each task   
  rcode2 = nf90_inq_varid (cdfid2, 'time', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'time',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  rcode2 = nf90_get_att (cdfid2, varid, 'units', date_init)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SIMULATION_START_DATE',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF  
  mcip_rd = date_init(13:31) // '.0000'
  mcip_rd(11:11)='-'  
  rcode2=nf90_get_var(cdfid2,varid,rdtime)
  IF ( rcode /= nf90_noerr ) THEN
   write(*,*)'error getting time in SFC file',str3
   CALL graceful_stop (pname)
  ENDIF

  CALL geth_newdate (mcip_next, mcip_rd, int(rdtime)*intvl*60)
  if(mcip_next.ne.mcip_now) then
   write(*,*)'time mismatch in SFC file ',mcip_now,mcip_next,mcip_rd,date_init,rdtime
   CALL graceful_stop (pname)
  ENDIF   

!-------------------------------------------------------------------------------
! Read FV3 data for this domain.
!-------------------------------------------------------------------------------
  it=1
  CALL get_var_3d_real_cdf (cdfid, 'ugrd', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,utmp,xdindex,ydindex,ncols_x+1,nrows_x+1,2)  ! put it into Dot point for later rotation
    kk=met_nz-k+1                                            ! flip to bottom up
    ua(1:ncols_x+1,1:nrows_x+1,kk) = utmp(1:ncols_x+1,1:nrows_x+1)
   enddo 
   WRITE (*,ifmt1) 'ugrd     ', (ua(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'ugrd', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'vgrd', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
  do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,utmp,xdindex,ydindex,ncols_x+1,nrows_x+1,2)
    
    kk=met_nz-k+1
    va(1:ncols_x+1,1:nrows_x+1,kk) = utmp(1:ncols_x+1,1:nrows_x+1)
     call windrotation(ua(1:ncols_x+1,1:nrows_x+1,kk),va(1:ncols_x+1,1:nrows_x+1,kk),londot(1:ncols_x+1,1:nrows_x+1), &
                            ncols_x+1,nrows_x+1,sngl(p_gam_gd),0.5*sngl(p_alp_gd+p_bet_gd))
   enddo
   WRITE (*,ifmt1) 'vgrd     ', (va(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'vgrd', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'dzdt', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    wa(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1a) 'dzdt     ', (wa(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'dzdt', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  IF ( .NOT. ALLOCATED ( dum1d ) )  &
    ALLOCATE ( dum1d (met_nz ) )    ! 3D, cross points, half lvls
 
  CALL get_var_1d_real_cdf (cdfid, 'pfull', dum1d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
  phalf = dum1d(met_nz:1:-1) *100.0  !FV3 mb (hPa) to Pa
  ELSE
    WRITE (*,f9400) TRIM(pname), 'phalf', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
   DEALLOCATE (dum1d)

  IF ( .NOT. ALLOCATED ( dum1d ) )  &
    ALLOCATE ( dum1d (met_nz + 1 ) )    ! 3D, cross points, full lvls

  CALL get_var_1d_real_cdf (cdfid, 'phalf', dum1d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
  pfull = dum1d(met_nz+1:1:-1) *100.0  !FV3 mb (hPa) to Pa
  ELSE
    WRITE (*,f9400) TRIM(pname), 'pfull', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
   DEALLOCATE (dum1d)
   
!Calculate sigma from pressure layers
  
  sigmaf = (pfull - pfull(met_nz+1)) / (pfull(1) - pfull(met_nz+1))
  sigmah = (phalf - phalf(met_nz)) / (phalf(1) - phalf(met_nz))
  
  CALL get_var_3d_real_cdf (cdfid, 'dpres', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    dpres(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1a) 'dpres      ', (dpres(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'dpres', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'delz', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    delz(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1a) 'delz      ', (delz(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'delz', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'tmp', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    ta(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1a) 'tmp      ', (ta(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'tmp', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'spfh', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    qva(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1) 'spfh     ', (qva(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'spfh', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'clwmr', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    qca(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
   WRITE (*,ifmt1) 'clwmr    ', (qca(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'clwmr', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_3d_real_cdf (cdfid, 'rwmr', dum3d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
   do k=1,met_nz
    call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    kk=met_nz-k+1
    qra(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
   enddo
    WRITE (*,ifmt1) 'rwmr     ', (qra(lprt_metx,lprt_mety,k),k=1,met_nz)
  ELSE
    WRITE (*,f9400) TRIM(pname), 'rwmr', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'icmr', rcode)
  IF ( rcode == nf90_noerr ) THEN
    CALL get_var_3d_real_cdf (cdfid, 'icmr', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     do k=1,met_nz
      call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      kk=met_nz-k+1
      qia(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
     enddo
      WRITE (*,ifmt1) 'icmr     ', (qia(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'icmr', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    qia(:,:,:) = 0.0
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'snmr', rcode)
  IF ( rcode == nf90_noerr ) THEN
    CALL get_var_3d_real_cdf (cdfid, 'snmr', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     do k=1,met_nz
      call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      kk=met_nz-k+1
      qsa(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
     enddo
     WRITE (*,ifmt1) 'snmr     ', (qsa(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'snmr', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    qsa(:,:,:) = 0.0
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'grle', rcode)
  IF ( rcode == nf90_noerr ) THEN
    CALL get_var_3d_real_cdf (cdfid, 'grle', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     do k=1,met_nz
      call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      kk=met_nz-k+1
      qga(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
     enddo
      WRITE (*,ifmt1) 'grle     ', (qga(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'grle', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    qga(:,:,:) = 0.0
  ENDIF

!FV3 does not have TKE
  IF ( ( iftke ) .AND. ( iftkef ) ) THEN  ! TKE on full-levels
    CALL get_var_3d_real_cdf (cdfid, 'TKE', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     do k=1,met_nz
      call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      kk=met_nz-k+1
      tke(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
     enddo
      WRITE (*,ifmt1a) 'TKE      ', (tke(lprt_metx,lprt_mety,k),k=1,nzp)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'TKE', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ELSE IF ( ( iftke ) .AND. ( .NOT. iftkef ) ) THEN  ! TKE on half-layers
    CALL get_var_3d_real_cdf (cdfid, 'TKE_MYJ', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      do k=1,met_nz
       call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
       kk=met_nz-k+1
       tke(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
      enddo
      WRITE (*,ifmt1) 'TKE_MYJ  ', (tke(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'TKE_MYJ', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  IF ( ifcld3d ) THEN  ! 3D resolved cloud fraction
    CALL get_var_3d_real_cdf (cdfid, 'cld_amt', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     do k=1,met_nz
      call myinterp(dum3d(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      kk=met_nz-k+1
      cldfra(1:ncols_x,1:nrows_x,kk) = atmp(1:ncols_x,1:nrows_x)
     enddo
     WRITE (*,ifmt1a) 'cld_amt  ', (cldfra(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'cld_amt', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  IF ( ift2m ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'tmp2m', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
     t2(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
     WRITE (*,f6000) 'tmp2m      ', t2(lprt_metx, lprt_mety), 'K'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'tmp2m', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  IF ( ifq2m ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'spfh2m', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      q2(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'spfh2m   ', q2(lprt_metx, lprt_mety), 'kg kg-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'spfh2m', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  IF ( ifw10m ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'ugrd10m', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      u10(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'ugrd10m  ', u10(lprt_metx, lprt_mety), 'm s-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'ugrd10m', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
    CALL get_var_2d_real_cdf (cdfid2, 'vgrd10m', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      v10(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      call windrotation(u10(1:ncols_x,1:nrows_x),v10(1:ncols_x,1:nrows_x),loncrs(1:ncols_x,1:nrows_x), & 
                        ncols_x,nrows_x,sngl(p_gam_gd),0.5*sngl(p_alp_gd+p_bet_gd))
      WRITE (*,f6000) 'vgrd10m  ', v10(lprt_metx, lprt_mety), 'm s-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'vgrd10m', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'pressfc', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    psa(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'pressfc       ', psa(lprt_metx, lprt_mety), 'Pa'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'pressfc', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

!Assume at surface the FV3 geopotential height (gpm) = geometric height (m)
  CALL get_var_2d_real_cdf (cdfid2, 'orog', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    terrain(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'orog      ', terrain(lprt_metx, lprt_mety), 'm'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'orog', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

    ! incremental ave precip taken directly from FV3 (avoid IF block)

    CALL get_var_2d_real_cdf (cdfid2, 'cprat_ave', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      WHERE ( dum2d < smallnum )
        dum2d = 0.0
      ENDWHERE
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
!      if(nn.eq.1) then
!       raincon(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!      else
!       raincon(1:ncols_x,1:nrows_x) = amax1(0., atmp(1:ncols_x,1:nrows_x)*nn-  &
!          rcold(1:ncols_x,1:nrows_x)*(nn-1))
!      endif
!      rcold(1:ncols_x,1:nrows_x)=atmp(1:ncols_x,1:nrows_x)
      	  
      !Convert mass precip rate in FV3 (kg/m2/s) to column amount (cm/hour)
!      raincon(1:ncols_x,1:nrows_x)=raincon(1:ncols_x,1:nrows_x)/997.0 * 100.0*3600.
      
      raincon(1:ncols_x,1:nrows_x)= atmp(1:ncols_x,1:nrows_x)/997.0 * 100.0*3600.
      
      WRITE (*,f6000) 'cprat_ave ', raincon(lprt_metx, lprt_mety), 'cm'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'cprat_ave', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'prate_ave', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      WHERE ( dum2d < smallnum )
        dum2d = 0.0
      ENDWHERE
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
!      if(nn.eq.1) then
!       rainnon(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!      else
!       rainnon(1:ncols_x,1:nrows_x) = amax1(0., atmp(1:ncols_x,1:nrows_x)*nn-  &
!          rnold(1:ncols_x,1:nrows_x)*(nn-1))
!      endif
!      rnold(1:ncols_x,1:nrows_x)=atmp(1:ncols_x,1:nrows_x)
      	  
      !Convert mass precip rate in FV3 (kg/m2/s) to column amount (cm/hour)
!      rainnon(1:ncols_x,1:nrows_x)=rainnon(1:ncols_x,1:nrows_x)/997.0 * 100.0*3600.
       rainnon(1:ncols_x,1:nrows_x)=atmp(1:ncols_x,1:nrows_x)/997.0 * 100.0*3600.
      WRITE (*,f6000) 'prate_ave  ', rainnon(lprt_metx, lprt_mety), 'cm'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'prate_ave', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'dswrf', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    rgrnd(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'dswrf   ', rgrnd(lprt_metx, lprt_mety), 'W m-2'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'dswrf', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'dlwrf', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    glw(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'dlwrf      ', glw(lprt_metx, lprt_mety), 'W m-2'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'dlwrf', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF


  CALL get_var_2d_real_cdf (cdfid2, 'vtype', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    IF ( MAXVAL(dum2d) > nummetlu ) THEN
      WRITE (*,f9500) TRIM(pname), met_lu_src, MAXVAL(dum2d)
      CALL graceful_stop (pname)
    ENDIF
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    landuse(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6100) 'vtype ', landuse(lprt_metx, lprt_mety), 'category'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'vtype', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

    WHERE ( INT(landuse) == met_lu_water_fv3 )  ! FV3 water = 0
      landuse = met_lu_water ! MODIS IGBP water = 17
    ELSEWHERE  ! land
      landuse = landuse
    END WHERE


  CALL get_var_2d_real_cdf (cdfid2, 'land', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    landmask(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'land ', landmask(lprt_metx, lprt_mety), 'category'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'land', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'shtfl', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    hfx(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'shtfl      ', hfx(lprt_metx, lprt_mety), 'W m-2'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'shtfl', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'lhtfl', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    lh(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'lhtfl       ', lh(lprt_metx, lprt_mety), 'W m-2'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'lhtfl', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'fricv', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    ust(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'fricv      ', ust(lprt_metx, lprt_mety), 'm s-1'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'fricv', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  !M-O length not in FV3GFSv16
  IF ( ifmol ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'RMOL', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      mol(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'RMOL', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
    IF ( met_urban_phys >= 1 ) THEN  ! UCM used; get MOL above urban canopy
      CALL get_var_2d_real_cdf (cdfid2, 'XXXC_URB', dum2d, it, rcode)
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      IF ( rcode == nf90_noerr ) THEN  ! blend urban M-O length with RMOL
        IF ( ( met_lu_src(1:4) == 'USGS' ) .AND.  &
             ( MAXVAL(landuse)  >  24    ) ) THEN  ! 33-category USGS/NLCD
          DO j = 1, nym
            DO i = 1, nxm
              IF ( ( landuse(i,j) ==  1 ) .OR. ( landuse(i,j) == 31 ) .OR.  &
                   ( landuse(i,j) == 32 ) .OR. ( landuse(i,j) == 33 ) ) THEN
                mol(i,j) = atmp(i,j)  ! XXXC_URB is not inverted
              ENDIF
            ENDDO
          ENDDO
        ELSE IF ( met_lu_src(1:4) == 'USGS' ) THEN  ! 24-category USGS
          DO j = 1, nym
            DO i = 1, nxm
              IF ( landuse(i,j) == 1 ) THEN  ! urban
                mol(i,j) = atmp(i,j)  ! XXXC_URB is not inverted
              ENDIF
            ENDDO
          ENDDO
        ELSE
          WRITE (*,f9800) TRIM(pname), 'XXXC_URB (URBAN MOL)', met_lu_src(1:4)
          CALL graceful_stop (pname)
        ENDIF
      ELSE
!~~~    Just use RMOL to fill Monin-Obukhov length without extra urban field
      ENDIF
    ENDIF
    WRITE (*,f6000) 'MOL      ', mol(lprt_metx, lprt_mety), 'm'
  ENDIF
! No P-X physics in FV3. :)
  IF ( ifmolpx ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'QFX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      qfx(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'QFX      ', qfx(lprt_metx, lprt_mety), 'kg m-2 s-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'QFX ', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'hpbl', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    zpbl(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'hpbl     ', zpbl(lprt_metx, lprt_mety), 'm'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'hpbl', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
! No stomatal resistance in FV3 yet.  :(
  IF ( ifresist ) THEN

    CALL get_var_2d_real_cdf (cdfid2, 'RA', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      ra(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      IF ( ABS(MAXVAL(ra)) < smallnum ) THEN
        ifresist = .FALSE.
      ENDIF
      WRITE (*,f6000) 'RA       ', ra(lprt_metx, lprt_mety), 's m-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'RA', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'RS', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      rstom(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      IF ( ABS(MAXVAL(rstom)) < smallnum ) THEN
        ifresist = .FALSE.
      ENDIF
      WRITE (*,f6000) 'RS       ', rstom(lprt_metx, lprt_mety), 's m-1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'RS', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  ENDIF

    IF ( ifclayf ) THEN
      IF ( ifclayfwrfout ) THEN  ! clayf in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'CLAY_FRAC', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           clayf(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!        IF ( ABS(MAXVAL(clayf)) < smallnum ) THEN
!          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
!            clayf(:,:) = 0.1
!          ENDIF
!        ENDIF
          WRITE (*,ifmt2) 'CLAYF      ',(clayf(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'CLAYF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! clay fraction in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'CLAY_FRAC', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          clayf(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! CLAYF check over  water, set as negative numbers for improved error checking
!          WHERE ( (INT(landmask) == 0) .OR. (clayf > 1.0) ) ! FV3 water = 0 or frac > 1, set CLAYF < 0.0
!           clayf = -1.0
!          END WHERE
          WRITE (*,ifmt2) 'CLAYF ', clayf(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'CLAYF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF

    IF ( ifsandf ) THEN
      IF ( ifsandfwrfout ) THEN  ! sandf in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'SAND_FRAC', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           sandf(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!        IF ( ABS(MAXVAL(sandf)) < smallnum ) THEN
!          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
!            sandf(:,:) = 0.1
!          ENDIF
!        ENDIF
          WRITE (*,ifmt2) 'SANDF      ',(sandf(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'SANDF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! sand fraction in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'SAND_FRAC', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          sandf(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! SANDF check over  water, set as negative numbers for improved error checking
!          WHERE ( (INT(landmask) == 0) .OR. (sandf > 1.0) ) ! FV3 water = 0 or frac > 1, set SANDF < 0.0
!           sandf = -1.0
!          END WHERE
          WRITE (*,ifmt2) 'SANDF ', sandf(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'SANDF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF
 
    IF ( ifdrag ) THEN
      IF ( ifdragwrfout ) THEN  ! drag in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'DRAG_PART', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           drag(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!        IF ( ABS(MAXVAL(drag)) < smallnum ) THEN
!          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
!            drag(:,:) = 1.0e-6
!          ENDIF
!        ENDIF
          WRITE (*,ifmt2) 'DRAG      ',(drag(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'DRAG', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! sand fraction in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'DRAG_PART', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          drag(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! DRAG check over  water, set as negative numbers for improved error checking
!          WHERE ( (INT(landmask) == 0) ) ! FV3 water = 0 and DRAG < 0.0
!           drag = -1.0
!          END WHERE
          WRITE (*,ifmt2) 'DRAG ', drag(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'DRAG', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF

    IF ( ifssm ) THEN
      IF ( ifssmwrfout ) THEN  ! ssm in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'SSM', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           ssm(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!        IF ( ABS(MAXVAL(ssm)) < smallnum ) THEN
!          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
!            ssm(:,:) = 1.0e-6
!          ENDIF
!        ENDIF
          WRITE (*,ifmt2) 'SSM      ',(ssm(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'SSM', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! ssm in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'SSM', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          ssm(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! SSM check over  water, set as negative numbers for improved error checking
!          WHERE ( (INT(landmask) == 0) ) ! FV3 water = 0 and CLAYF < 0.0
!           ssm = -1.0
!          END WHERE
          WRITE (*,ifmt2) 'SSM ', ssm(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'SSM', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF

    IF ( ifuthr ) THEN
      IF ( ifuthrwrfout ) THEN  ! uthr in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'UTHRES', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           uthr(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
!        IF ( ABS(MAXVAL(ssm)) < smallnum ) THEN
!          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
!            ssm(:,:) = 1.0e-6
!          ENDIF
!        ENDIF
          WRITE (*,ifmt2) 'UTHR      ',(uthr(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'UTHR', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! uthr in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'UTHRES', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          uthr(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! SSM check over  water, set as negative numbers for improved error checking
!          WHERE ( (INT(landmask) == 0) ) ! FV3 water = 0 and CLAYF < 0.0
!           ssm = -1.0
!          END WHERE
          WRITE (*,ifmt2) 'UTHR ', uthr(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'UTHR', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF


    IF ( iflai ) THEN
      IF ( iflaiwrfout ) THEN  ! leaf area index in FV3 history file
        CALL get_var_2d_real_cdf (cdfid2, 'LAI', dum2d, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
           call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           lai(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
        IF ( ABS(MAXVAL(lai)) < smallnum ) THEN
          IF ( met_soil_lsm == 2 ) THEN  ! NOAH LSM
            lai(:,:) = 4.0
          ENDIF
        ENDIF
          WRITE (*,ifmt2) 'LAI      ',(lai(lprt_metx,lprt_mety))
        ELSE
          WRITE (*,f9400) TRIM(pname), 'LAI', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
      ELSE  ! leaf area index in GEOGRID file from WPS
       flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_2d_real_cdf (cdfidg, 'LAI', dum2d, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
          call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
          lai(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
          ! Another LAI Check in case LAI=0 over land for processed satellite inputs
          ! Set to average LAI=4 for the representative pixels.
          WHERE ( (INT(landmask) == 1) .AND. (lai <= 0.0))  ! FV3 land = 1 and LAI = 0.0
           lai = 4.0
          END WHERE
          WRITE (*,ifmt2) 'LAI ', lai(lprt_metx,lprt_mety)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'LAI', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF

  IF ( ifwr ) THEN
    CALL get_var_2d_real_cdf (cdfid2, 'cnwat', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
      wr(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'cnwat   ', wr(lprt_metx, lprt_mety), 'kg m-2'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'cnwat', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF
  IF ( ifveg ) THEN
      CALL get_var_2d_real_cdf (cdfid2, 'veg', dum2d, it, rcode)
      IF ( rcode == nf90_noerr ) THEN
        call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
        veg(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)*0.01
        WRITE (*,f6000) 'veg   ', veg(lprt_metx, lprt_mety), 'fraction'
      ELSE
        WRITE (*,f9400) TRIM(pname), 'veg', TRIM(nf90_strerror(rcode))
        CALL graceful_stop (pname)
      ENDIF
  ENDIF

  IF ( ifsoil ) THEN

    CALL get_var_2d_real_cdf (cdfid2, 'sotyp', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      isltyp(1:ncols_x,1:nrows_x) = int(atmp(1:ncols_x,1:nrows_x))
      !Fix for isltyp = 0 for water in FV3GFS16 (SLTYP = 0 not allowed in CMAQ)
      !Set to isltyp = 14 as in for the WRFv4 16-category soil types.
      WHERE ( isltyp == 0 )
        isltyp = 14
      ENDWHERE
      WRITE (*,f6100) 'sotyp   ', isltyp(lprt_metx, lprt_mety), 'category'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'sotyp', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

!    Note the top two soil layers in FV3GFSv16 are at 0-10 cm and 10-40 cm
!    Noah LSM
!    Will need CMAQ adjustment to 0-1 cm and 1-10 cm.
    CALL get_var_2d_real_cdf (cdfid2, 'soilw1', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      wg(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      soim3d(1:ncols_x,1:nrows_x,1) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'soilw1  ', wg(lprt_metx, lprt_mety), 'm3 m-3'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilw1', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'soilw2', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      w2(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'soilw2  ', w2(lprt_metx, lprt_mety), 'm3 m-3'
      
      soim3d(1:ncols_x,1:nrows_x,2) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilw2', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'soilw3', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      soim3d(1:ncols_x,1:nrows_x,3) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilw3', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
 
    CALL get_var_2d_real_cdf (cdfid2, 'soilw4', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
     call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
     soim3d(1:ncols_x,1:nrows_x,4) = atmp(1:ncols_x,1:nrows_x)
     WRITE (*,ifmt5) 'soim3d    ', (soim3d(lprt_metx,lprt_mety,k),k=1,met_ns)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilw4', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    !    Note the top two soil layers in FV3GFSv16 are at 0-10 cm and 10-40 cm
    !    Noah LSM
    CALL get_var_2d_real_cdf (cdfid2, 'soilt1', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      soilt1(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'soilt1  ', soilt1(lprt_metx, lprt_mety), 'K'
      
      soit3d(1:ncols_x,1:nrows_x,1) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilt1', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'soilt2', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      soilt2(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'soilt2  ', soilt2(lprt_metx, lprt_mety), 'K'
      
      soit3d(1:ncols_x,1:nrows_x,2) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilt2', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'soilt3', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      soit3d(1:ncols_x,1:nrows_x,3) = atmp(1:ncols_x,1:nrows_x)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilt3', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'soilt4', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      soit3d(1:ncols_x,1:nrows_x,4) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,ifmt5) 'soit3d    ', (soit3d(lprt_metx,lprt_mety,k),k=1,met_ns)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'soilt4', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'tmpsfc', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    groundt(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'tmpsfc      ', groundt(lprt_metx, lprt_mety), 'K'
  ELSE
    IF ( ifsoil ) THEN
      groundt(:,:) = soilt1(:,:)
      WRITE (*,f6000) 'tmpsfc      ', groundt(lprt_metx, lprt_mety), 'K'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'tmpsfc', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  CALL get_var_2d_real_cdf (cdfid2, 'albdo_ave', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    albedo(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)*0.01
    WRITE (*,f6000) 'albdo_ave   ', albedo(lprt_metx, lprt_mety), 'fraction'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'albdo_ave', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  IF ( first ) THEN
    IF ( iflufrc ) THEN
      IF ( ifluwrfout ) THEN  ! land use fractions in WRF history file
        CALL get_var_3d_real_cdf (cdfid2, 'LANDUSEF', dum3d_l, it, rcode)
        IF ( rcode == nf90_noerr ) THEN
          do k=1,nummetlu
           call myinterp(dum3d_l(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           lufrac(1:ncols_x,1:nrows_x,k) = atmp(1:ncols_x,1:nrows_x)
          enddo
          WRITE (*,ifmt2) 'LANDUSEF ', (lufrac(lprt_metx,lprt_mety,k),k=1,nummetlu)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'LANDUSEF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        IF ( iflu2wrfout ) THEN  ! land use fractions (ranked) in WRF file
          CALL get_var_3d_real_cdf (cdfid2, 'LANDUSEF2', dum3d_l, it, rcode)
          IF ( rcode == nf90_noerr ) THEN
           do k=1,nummetlu
            call myinterp(dum3d_l(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
            lufrac2(1:ncols_x,1:nrows_x,k) = atmp(1:ncols_x,1:nrows_x)
           enddo
            WRITE (*,ifmt2) 'LANDUSEF2', (lufrac2(lprt_metx,lprt_mety,k),k=1,nummetlu)
          ELSE
            WRITE (*,f9400) TRIM(pname), 'LANDUSEF2', TRIM(nf90_strerror(rcode))
            CALL graceful_stop (pname)
          ENDIF
          CALL get_var_3d_int_cdf (cdfid2, 'MOSAIC_CAT_INDEX', dum3d_li, it, rcode)
          IF ( rcode == nf90_noerr ) THEN
           do k=1,nummetlu
            call myinterp(real(dum3d_li(:,:,k)),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
            moscatidx(1:ncols_x,1:nrows_x,k) = int(atmp(1:ncols_x,1:nrows_x))
           enddo 
            WRITE (*,ifmt3) 'MOSAIC_CAT', (moscatidx(lprt_metx,lprt_mety,k),k=1,nummetlu)
          ELSE
            ! Will be filled in getluse.f90, if NOAH Mosaic LSM was used
          ENDIF
        ENDIF
      ELSE  ! land use fractions in GEOGRID file from WPS
        flg = file_geo
        rcode = nf90_open (flg, nf90_nowrite, cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9900) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
        CALL get_var_3d_real_cdf (cdfidg, 'LANDUSEF', dum3d_l, 1, rcode)
        IF ( rcode == nf90_noerr ) THEN
         do k=1,nummetlu
           call myinterp(dum3d_l(:,:,k),met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
           lufrac(1:ncols_x,1:nrows_x,k) = atmp(1:ncols_x,1:nrows_x)
         enddo
          WRITE (*,ifmt2) 'LANDUSEF ', lufrac(lprt_metx,lprt_mety,:)
        ELSE
          WRITE (*,f9400) TRIM(pname), 'LANDUSEF', TRIM(nf90_strerror(rcode))
          CALL graceful_stop (pname)
        ENDIF
        rcode = nf90_close (cdfidg)
        IF ( rcode /= nf90_noerr ) THEN
          WRITE (*,f9950) TRIM(pname)
          CALL graceful_stop (pname)
        ENDIF
      ENDIF
    ENDIF
! No urban canopy in FV3 yet.  :(
    IF ( met_urban_phys >= 1 ) THEN  ! urban canopy model used
      CALL get_var_2d_real_cdf (cdfid2, 'FRC_URB2D', dum2d, it, rcode)
      IF ( rcode == nf90_noerr ) THEN
        frc_urb(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
        frc_urb(met_nx,:) = frc_urb(nxm,:)
        frc_urb(:,met_ny) = frc_urb(:,nym)
        WRITE (*,f6000) 'FRC_URB2D', frc_urb(lprt_metx, lprt_mety), 'fraction'
      ELSE
        WRITE (*,f9400) TRIM(pname), 'FRC_URB2D', TRIM(nf90_strerror(rcode))
      ENDIF
    ENDIF
    IF ( lpv > 0 ) THEN
      CALL get_var_2d_real_cdf (cdfid2, 'F', dum2d, it, rcode)
      IF ( rcode == nf90_noerr ) THEN
        coriolis(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
        coriolis(met_nx,:) = coriolis(nxm,:)
        coriolis(:,met_ny) = coriolis(:,nym)
        WRITE (*,f6000) 'F        ', coriolis(lprt_metx, lprt_mety), 's-1'
      ELSE
        WRITE (*,f9400) TRIM(pname), 'F      ', TRIM(nf90_strerror(rcode))
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF

  IF ( ifznt ) THEN  ! expecting roughness length in file
    CALL get_var_2d_real_cdf (cdfid2, 'sfcr', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
      znt(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
      WRITE (*,f6000) 'sfcr      ', znt(lprt_metx, lprt_mety),    'm'
      gotznt = .TRUE.
    ELSE
      WRITE (*,f9400) TRIM(pname), 'sfcr    ', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    gotznt = .FALSE.
  ENDIF
! FV3 does not have snow cover flag in output, use snow cover fraction instead
  CALL get_var_2d_real_cdf (cdfid2, 'snowc_ave', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    snowcovr(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)*0.01
    WRITE (*,f6000) 'snowc_ave    ', snowcovr(lprt_metx, lprt_mety), 'fraction'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'snowc_ave', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

! FV3 does not have sea ice fraction in output, use sea-ice thickness instead

  CALL get_var_2d_real_cdf (cdfid2, 'icetk', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,1)
    seaice(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    gotseaice = .TRUE.
    WRITE (*,f6000) 'icetk   ', seaice(lprt_metx, lprt_mety), 'XXX'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'icetk', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

! Where seaice thickness in FV3GFS16 > 0 set to 1.0 for grid cell 
! (CMAQ expects seaice fraction, not thickness)
    WHERE ( seaice >  0.0 )
        seaice = 1.0
      ENDWHERE

  CALL get_var_2d_real_cdf (cdfid2, 'snod', dum2d, it, rcode)
  IF ( rcode == nf90_noerr ) THEN
    call myinterp(dum2d,met_nx,met_ny,atmp,xindex,yindex,ncols_x,nrows_x,2)
    snowh(1:ncols_x,1:nrows_x) = atmp(1:ncols_x,1:nrows_x)
    WRITE (*,f6000) 'snod    ', snowh(lprt_metx, lprt_mety), 'm'
  ELSE
    WRITE (*,f9400) TRIM(pname), 'snod', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
! NO PX physics in FV3 for soil hydraulic paramters
  IF ( ifpxwrf41 ) THEN

    CALL get_var_2d_real_cdf (cdfid2, 'WSAT_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      wsat_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      wsat_px(met_nx,:) = wsat_px(nxm,:)
      wsat_px(:,met_ny) = wsat_px(:,nym)
      WRITE (*,f6000) 'WSAT_PX  ', wsat_px(lprt_metx, lprt_mety), 'm3 m-3'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'WSAT_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'WFC_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      wfc_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      wfc_px(met_nx,:) = wfc_px(nxm,:)
      wfc_px(:,met_ny) = wfc_px(:,nym)
      WRITE (*,f6000) 'WFC_PX  ', wfc_px(lprt_metx, lprt_mety), 'm3 m-3'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'WFC_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'WWLT_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      wwlt_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      wwlt_px(met_nx,:) = wwlt_px(nxm,:)
      wwlt_px(:,met_ny) = wwlt_px(:,nym)
      WRITE (*,f6000) 'WWLT_PX  ', wwlt_px(lprt_metx, lprt_mety), 'm3 m-3'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'WWLT_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'CSAND_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      csand_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      csand_px(met_nx,:) = csand_px(nxm,:)
      csand_px(:,met_ny) = csand_px(:,nym)
      WRITE (*,f6000) 'CSAND_PX ', csand_px(lprt_metx, lprt_mety), '1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'CSAND_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'FMSAND_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      fmsand_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      fmsand_px(met_nx,:) = fmsand_px(nxm,:)
      fmsand_px(:,met_ny) = fmsand_px(:,nym)
      WRITE (*,f6000) 'FMSAND_PX', fmsand_px(lprt_metx, lprt_mety), '1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'FMSAND_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'CLAY_PX', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      clay_px(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      clay_px(met_nx,:) = clay_px(nxm,:)
      clay_px(:,met_ny) = clay_px(:,nym)
      WRITE (*,f6000) 'CLAY_PX  ', clay_px(lprt_metx, lprt_mety), '1'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'CLAY_PX', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  ENDIF
! No LU mosaic fractional model in FV3 yet... :(
  IF ( ifmosaic ) THEN

    CALL get_var_3d_real_cdf (cdfid, 'LAI_MOSAIC', dum3d_m, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      lai_mos(1:nxm,   1:nym,   :) = dum3d_m(:,:,:)
      lai_mos(  met_nx, :,      :) = lai_mos(nxm,:,:)
      lai_mos( :,        met_ny,:) = lai_mos(:,nym,:)
      WRITE (*,ifmt4) 'LAI_MOS  ', (lai_mos(lprt_metx,lprt_mety,k),k=1,nummosaic)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'LAI_MOSAIC', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'RS_MOSAIC', dum3d_m, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      rs_mos(1:nxm,   1:nym,   :) = dum3d_m(:,:,:)
      rs_mos(  met_nx, :,      :) = rs_mos(nxm,:,:)
      rs_mos( :,        met_ny,:) = rs_mos(:,nym,:)
      WRITE (*,ifmt4) 'RS_MOS   ', (rs_mos(lprt_metx,lprt_mety,k),k=1,nummosaic)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'RS_MOSAIC', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'TSK_MOSAIC', dum3d_m, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      tsk_mos(1:nxm,   1:nym,   :) = dum3d_m(:,:,:)
      tsk_mos(  met_nx, :,      :) = tsk_mos(nxm,:,:)
      tsk_mos( :,        met_ny,:) = tsk_mos(:,nym,:)
      WRITE (*,ifmt4) 'TSK_MOS  ', (tsk_mos(lprt_metx,lprt_mety,k),k=1,nummosaic)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'TSK_MOSAIC', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'ZNT_MOSAIC', dum3d_m, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      znt_mos(1:nxm,   1:nym,   :) = dum3d_m(:,:,:)
      znt_mos(  met_nx, :,      :) = znt_mos(nxm,:,:)
      znt_mos( :,        met_ny,:) = znt_mos(:,nym,:)
      WRITE (*,ifmt4) 'ZNT_MOS  ', (znt_mos(lprt_metx,lprt_mety,k),k=1,nummosaic)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'ZNT_MOSAIC', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'WSPD', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      wspdsfc(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      wspdsfc(met_nx,:) = wspdsfc(nxm,:)
      wspdsfc(:,met_ny) = wspdsfc(:,nym)
      WRITE (*,f6000) 'WSPDSFC  ', wspdsfc(lprt_metx, lprt_mety), 'm s-1'
    ELSE
      ! Original version was stored in "WSPDSFC"; released WRF code uses "WSPD"
      CALL get_var_2d_real_cdf (cdfid2, 'WSPDSFC', dum2d, it, rcode)
      IF ( rcode == nf90_noerr ) THEN
        wspdsfc(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
        wspdsfc(met_nx,:) = wspdsfc(nxm,:)
        wspdsfc(:,met_ny) = wspdsfc(:,nym)
        WRITE (*,f6000) 'WSPDSFC  ', wspdsfc(lprt_metx, lprt_mety), 'm s-1'
      ELSE
        WRITE (*,f9400) TRIM(pname), 'WSPDSFC', TRIM(nf90_strerror(rcode))
        CALL graceful_stop (pname)
      ENDIF
    ENDIF

    CALL get_var_2d_real_cdf (cdfid2, 'XLAIDYN', dum2d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      xlaidyn(1:nxm,1:nym) = dum2d(:,met_ny:1:-1)
      xlaidyn(met_nx,:) = xlaidyn(nxm,:)
      xlaidyn(:,met_ny) = xlaidyn(:,nym)
      WRITE (*,f6000) 'XLAIDYN  ', xlaidyn(lprt_metx, lprt_mety), 'm2 m-2'
    ELSE
      WRITE (*,f9400) TRIM(pname), 'XLAIDYN', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  ENDIF  ! ifmosaic
!No KF radiative extras in FV3 yet... :(
  IF ( ifkfradextras ) THEN  ! Extra vars from KF scheme w radiative feedbacks

    CALL get_var_3d_real_cdf (cdfid, 'QC_CU', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      qc_cu(1:nxm,   1:nym,   :) = dum3d(:,met_ny:1:-1,met_nz:1:-1)
      qc_cu(  met_nx, :,      :) = qc_cu(nxm,:,:)
      qc_cu( :,        met_ny,:) = qc_cu(:,nym,:)
      WRITE (*,ifmt1a) 'QC_CU    ', (qc_cu(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'QC_CU', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'QI_CU', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      qi_cu(1:nxm,   1:nym,   :) = dum3d(:,met_ny:1:-1,met_nz:1:-1)
      qi_cu(  met_nx, :,      :) = qi_cu(nxm,:,:)
      qi_cu( :,        met_ny,:) = qi_cu(:,nym,:)
      WRITE (*,ifmt1a) 'QI_CU    ', (qi_cu(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'QI_CU', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'CLDFRA_DP', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      cldfra_dp(1:nxm,   1:nym,   :) = dum3d(:,met_ny:1:-1,met_nz:1:-1)
      cldfra_dp(  met_nx, :,      :) = cldfra_dp(nxm,:,:)
      cldfra_dp( :,        met_ny,:) = cldfra_dp(:,nym,:)
      WRITE (*,ifmt1a) 'CLDFRA_DP', (cldfra_dp(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'CLDFRA_DP', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    CALL get_var_3d_real_cdf (cdfid, 'CLDFRA_SH', dum3d, it, rcode)
    IF ( rcode == nf90_noerr ) THEN
      cldfra_sh(1:nxm,   1:nym,   :) = dum3d(:,met_ny:1:-1,met_nz:1:-1)
      cldfra_sh(  met_nx, :,      :) = cldfra_sh(nxm,:,:)
      cldfra_sh( :,        met_ny,:) = cldfra_sh(:,nym,:)
      WRITE (*,ifmt1a) 'CLDFRA_SH', (cldfra_sh(lprt_metx,lprt_mety,k),k=1,met_nz)
    ELSE
      WRITE (*,f9400) TRIM(pname), 'CLDFRA_SH', TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

  ENDIF  ! ifkfradextras

  IF ( met_hybrid >= 0 ) THEN
!FV3 does not have the hybrid coefficients and need to calculate them using
!b parameter, sigmah, an sigmaf
    
    IF ( .NOT. ALLOCATED ( dum1d ) )  &
    ALLOCATE ( dum1d (met_nz + 1 ) )    ! 3D, cross points, full lvls

    rcode = nf90_get_att (cdfid, nf90_global, 'bk', dum1d)
    IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'bk values for Hybrid',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
    ENDIF
    
    b_k(1:met_nz+1) = dum1d(met_nz+1:1:-1)

    c1f(1) = 1.0 
    DO k = 2, met_nz+1 
    c1f(k) = (b_k(k) - b_k(k-1)) / (sigmaf(k) - sigmaf(k-1))
    END DO
    c2f = (1.0-c1f) * (100000.0 - met_ptop)

    c1h(1) = 1.0
    DO k = 2, met_nz
    c1h(k) = (b_k(k) - b_k(k-1)) / (sigmah(k) - sigmah(k-1))
    END DO
    c2h = (1.0-c1h) * (100000.0 - met_ptop)

    DEALLOCATE (dum1d)

  ENDIF

!FV3 output doesn't have soil thicknesses, hardcode to 4 Noah LSM levels!
!Harcoded FV3GFSv16 Noah 4-layers soil thicknesses right now...
  dzs = (/0.1, 0.3, 0.6, 1.0/)

  rcode = nf90_close (cdfid)
  rcode2=nf90_close(cdfid2)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9950) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF


!-------------------------------------------------------------------------------
! If this is the first time in this routine, then determine season.
!-------------------------------------------------------------------------------
! But if global 
  IF ( first ) THEN

    ! These seasons are used in MM5 and WRF for land-use lookup tables.

    startseas = met_startdate(1:4) // "-04-15-00:00:00"
    endseas   = met_startdate(1:4) // "-10-15-00:00:00"

    CALL geth_idts (met_startdate(1:19), startseas,           idts_start)
    CALL geth_idts (endseas,             met_startdate(1:19), idts_end)

    IF ( ( idts_start < 0 ) .OR. ( idts_end < 0 ) ) THEN
      IF ( met_cen_lat >= 0.0 ) THEN  ! Northern Hemisphere
        met_season = 2   ! winter
      ELSE  ! Southern Hemisphere
        met_season = 1   ! summer
      ENDIF
    ELSE
      IF ( met_cen_lat >= 0.0 ) THEN  ! Northern Hemisphere
        met_season = 1   ! summer
      ELSE  ! Southern Hemisphere
        met_season = 2   ! winter
      ENDIF
    ENDIF
!-------------------------------------------------------------------------------
! If roughness length was not available in output, fill it from lookup tables.
! If the urban model was used in WRF, replace roughness length with urban-
! specific arrays.
!-------------------------------------------------------------------------------

    IF ( .NOT. gotznt ) THEN

      IF ( met_season == 1 ) THEN  ! summer

        DO j = 1, nym
          DO i = 1, nxm
            IF ( ( met_lu_src(1:4) == "USGS" ) .AND.  &
                 ( met_lu_water == 16 ) ) THEN
              znt(i,j) = sfz0usgssum(landuse(i,j)) * 0.01  ! cm --> m
            ELSE IF ( ( met_lu_src(1:3) == "OLD" ) .AND.  &
                      ( met_lu_water == 7 ) ) THEN
              znt(i,j) = sfz0oldsum(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:6) == "NLCD50" ) THEN
              znt(i,j) = sfz0nlcd50sum(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:6) == "NLCD40" ) THEN
              znt(i,j) = sfz0nlcd40sum(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:3) == "SIB" ) THEN
              znt(i,j) = sfz0sibsum(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:3) == "MOD" ) THEN
              znt(i,j) = sfz0modsum(landuse(i,j))  * 0.01  ! cm --> m
            ELSE
              WRITE (*,f9700) TRIM(pname), met_lu_src, met_lu_water
              CALL graceful_stop (pname)
            ENDIF
          ENDDO
        ENDDO

      ELSE IF ( met_season == 2 ) THEN  ! winter

        DO j = 1, nym
          DO i = 1, nxm
            IF ( ( met_lu_src(1:4) == "USGS" ) .AND.  &
                 ( met_lu_water == 16 ) ) THEN
              znt(i,j) = sfz0usgswin(landuse(i,j)) * 0.01  ! cm --> m
            ELSE IF ( ( met_lu_src(1:3) == "OLD" ) .AND.  &
                      ( met_lu_water == 7 ) ) THEN
              znt(i,j) = sfz0oldwin(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:6) == "NLCD50" ) THEN
              znt(i,j) = sfz0nlcd50win(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:6) == "NLCD40" ) THEN
              znt(i,j) = sfz0nlcd40win(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:3) == "SIB" ) THEN
              znt(i,j) = sfz0sibwin(landuse(i,j))  * 0.01  ! cm --> m
            ELSE IF ( met_lu_src(1:3) == "MOD" ) THEN
              znt(i,j) = sfz0modwin(landuse(i,j))  * 0.01  ! cm --> m
            ELSE
              WRITE (*,f9700) TRIM(pname), met_lu_src, met_lu_water
              CALL graceful_stop (pname)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

      znt(:,met_ny) = znt(:,nym)
      znt(met_nx,:) = znt(nxm,:)

      IF ( met_urban_phys < 1 ) THEN  ! if UCM, write after urban update
        WRITE (*,f6000) 'ZNT      ', znt   (lprt_metx, lprt_mety), 'm'
      ENDIF

    ENDIF

    first = .FALSE.

  ENDIF

!-------------------------------------------------------------------------------
! If sea ice was not part of the output, set flag to compute it later
! in METVARS2CTM.
!-------------------------------------------------------------------------------

  IF ( .NOT. gotseaice ) THEN
    WRITE (*,f9975) TRIM(pname), 'SEAICE'
    needseaice = .TRUE.
  ELSE
    needseaice = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
! Deallocate arrays.
!-------------------------------------------------------------------------------

! DEALLOCATE ( dum2d )    ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum2d_i )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum2d_u )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum2d_v )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_l )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_li ) ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_m )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_p )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_s )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_t )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_u)  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_v )  ! commented out to avoid memory fragmentation
! DEALLOCATE ( dum3d_w )  ! commented out to avoid memory fragmentation

!-------------------------------------------------------------------------------

END SUBROUTINE rdfv3
