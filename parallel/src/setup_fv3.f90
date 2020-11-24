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

SUBROUTINE setup_fv3 (cdfid, cdfid2, ctmlays)

!-------------------------------------------------------------------------------
! Name:     Set Up the FV3 Domain Attributes
! Purpose:  Establishes bounds for FV3 post-processing.
! Revised:  ?? Jun 2004  Modified from MCIP2.2 for FV3. (S.-B. Kim)
!           26 May 2005  Changed vertical dimension to reflect full-layer
!                        dimension in FV3 header.  Added dynamic calculation
!                        of MET_TAPFRQ.  Converted dimensions to X,Y as opposed
!                        to the (former) convention that aligned with MM5.
!                        Included updates from MCIPv2.3.  Added calculation of
!                        cone factor.  Added logic for moist species, 2-m
!                        temperature, and 10-m winds.  Added definitions for
!                        FV3 base state variables.  Added capability to use all
!                        FV3 layers for MCIP without defining a priori.
!                        Cleaned up code.  (T. Otte)
!           15 Jul 2005  Added debugging on variable retrievals.  Changed check
!                        on 3D mixing ratios from rain to ice.  Corrected RADM
!                        seasons for Southern Hemisphere.  Corrected variable
!                        name for retrieval of surface physics option. (T. Otte)
!           18 Aug 2005  Changed internal variable SIGN to FAC to avoid
!                        confusion with F90 intrinsic function.  (T. Otte)
!           10 Apr 2006  Corrected checking of I/O API variables for Mercator
!                        projection.  (T. Otte)
!           12 May 2006  Corrected setting of I/O API variables for polar
!                        stereographic projection.  Revised defining and
!                        setting projection variables for module METINFO.
!                        Added restriction on using Eta/Ferrier microphysics
!                        scheme where QCLOUD represents total condensate.
!                        (T. Otte)
!           20 Jun 2006  Changed setting of IDTSEC from REAL to INTEGER
!                        value.  (T. Otte)
!           27 Jul 2007  Removed settings for RADMdry variable ISESN and for
!                        MET_INHYD.  Updated read of P_TOP to account for new
!                        method of storing "real" scalars in FV3 I/O API with
!                        FV3.  Added checks for fractional land use, leaf
!                        area index, Monin-Obukhov length, aerodynamic and
!                        stomatal resistances, vegetation fraction, canopy
!                        wetness, and soil moisture, temperature, and type in
!                        FV3 file.  Added read for number of land use
!                        categories...new with FV3V2.2.  Added read for number
!                        of soil layers, MET_RELEASE, MET_FDDA_3DAN and
!                        MET_FDDA_OBS.  Set MET_FDDA_SFAN to 0 for now because
!                        that option is not in FV3 ARW as of V2.2.  Changed
!                        MET_RADIATION into MET_LW_RAD and MET_SW_RAD.
!                        (T. Otte)
!           06 May 2008  Changed criteria for setting NUMMETLU when netCDF
!                        dimension "land_cat_stag" does not exist.  Added
!                        checks to determine if 2-m mixing ratio (Q2) and
!                        turbulent kinetic energy (TKE_MYJ) arrays exist, and
!                        set flags appropriately.  Extract nudging coefficients
!                        from header to use in metadata.  Extract whether or
!                        not the urban canopy model was used.  (T. Otte)
!           27 Oct 2009  Cleaned up file opening and logging in FV3 I/O API to
!                        prevent condition with too many files open for long
!                        simulations.  Added MODIFIED IGBP MODIS NOAH and 
!                        NLCD/MODIS as land-use classification options.
!                        Changed MET_UCMCALL to MET_URBAN_PHYS, and allowed
!                        for variable to be set to be greater than 1.  Chnaged
!                        code to allow for surface analysis nudging option
!                        and coefficients to be defined per FV3.  Define
!                        MET_CEN_LAT, MET_CEN_LON, MET_RICTR_DOT, MET_RJCTR_DOT,
!                        and MET_REF_LAT.  Increased MAX_TIMES to 1000.  Compute
!                        MET_XXCTR and MET_YYCTR.  Corrected setting for
!                        DATE_INIT, and fill variable MET_RESTART.  Read number
!                        of land use categories from FV3 global attributes for
!                        FV3V3.1 and beyond.  Allow output from FV3
!                        Preprocessing System (WPS) routine, GEOGRID, to provide
!                        fractional land use output if it is unavailable in FV3
!                        output.  Fill MET_P_ALP_D and MET_P_BET_D here
!                        rather than in setgriddefs.F for Mercator.  Added
!                        new logical variables IFLUFV3OUT and IFZNT.  (T. Otte)
!           12 Feb 2010  Removed unused variables COMM and SYSDEP_INFO.
!                        (T. Otte)
!           18 Mar 2010  Added CDFID as an input argument, and no longer open
!                        and close FV3 history file here.  Added CDFIDG as an
!                        input argument for subroutine CHKWPSHDR.  (T. Otte)
!           15 Dec 2010  Improved support for long MCIP runs from long FV3
!                        runs by increasing MAX_TIMES to 9999.  Added
!                        MET_RAIN_BUCKET.  (T. Otte)
!           23 Feb 2011  Refined error checking for MET_RAIN_BUCKET.  (T. Otte)
!           11 Aug 2011  Added MET_SHAL_CU to input.  Replaced module PARMS3
!                        with I/O API module M3UTILIO.  (T. Otte)
!           24 Aug 2011  Changed name of module FILE to FILES to avoid conflict
!                        with F90 protected intrinsic.  Updated netCDF commands
!                        to F90, and improved error handling.  Replaced calls
!                        to GET_TIMES_CDF with explicit netCDF functions.
!                        (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           21 Nov 2011  Force 2-m water vapor mixing ratio from FV3 with
!                        YSU PBL to be filled with layer 1 QVAPOR to avoid
!                        occasional Q2 < 0 in wintertime.  (T. Otte)
!           07 Dec 2011  Removed requirement to fill nudging coefficient for
!                        moisture when spectral nudging is used in FV3; as of
!                        FV3, spectral nudging toward moisture is not
!                        released in FV3.  Also added provision to collect
!                        nudging coefficient for geopotential when spectral
!                        nudging is used; was added to FV3 header with FV3v3.2.
!                        (T. Otte)
!           21 Aug 2012  Added MET_PCP_INCR for FV3V3.2+.  (T. Otte)
!           10 Sep 2012  Added handling for 40-category 2006 NLCD-MODIS land
!                        use classification as "NLCD40".  Added alternate name
!                        for 50-category 2001 NLCD-MODIS land use classification
!                        as "NLCD50".  (T. Otte)
!           26 Nov 2014  Added reads of ice, lake, and urban land use indices,
!                        and moved those definitions from getluse.f90 to this
!                        routine.  (T. Spero)
!           10 Apr 2015  Determine if 3D cloud fraction is part of FV3 output
!                        and if it represents resolved clouds.  Fill new logical
!                        variable IFCLD3D appropriately so that if resolved
!                        cloud fraction is available, it will be passed through
!                        in output.  (T. Spero)
!           21 Aug 2015  Added flag to capture whether ACM2 was run so that
!                        Monin-Obukhov length can be recalculated following
!                        the "corrector" part of the predictor-corrector in
!                        FV3/ACM2.  (T. Spero)
!           17 Sep 2015  Changed IFMOLACM to IFMOLPX.  (T. Spero)
!           21 Apr 2017  Added MODIS category 21 as "Lake".  (T. Spero)
!           23 Jun 2017  Added a check for FV3's hybrid vertical coordinate
!                        in FV3v3.9 and beyond.  Currently disabled MCIP when
!                        that coordinate is detected.  To be implemented in
!                        a later release of MCIP.  (T. Spero)
!           09 Feb 2018  Added support for hybrid vertical coordinate in FV3
!                        output.  Added capability to read and process data
!                        from the NOAH Mosaic land-surface model.  (T. Spero)
!           26 Jun 2018  Changed name of module with netCDF IO to broaden its
!                        usage.  Now use netCDF tokens for missing data.
!                        (T. Spero)
!           14 Sep 2018  Removed support for MM5v3 input.  (T. Spero)
!           23 Nov 2018  Modify criteria to determine whether incremental
!                        precipitation is available in FV3 output.  FV3v4.0
!                        allows header variable PREC_ACC_DT to appear even if
!                        the accompanying precipitation fields are not in the
!                        output.  (T. Spero)
!           14 Dec 2018  Added flag (IFRCURB) to determine if fraction of urban
!                        area is obtained from urban canopy model.  (T. Spero)
!           10 May 2019  Removed layer collapsing.  (T. Spero)
!           18 Jun 2019  Added a flag (IFPXFV341) to determine of new surface
!                        variables with PX LSM are available to improve dust
!                        simulation in CCTM.  Added a flag (IFCURADFDBK) to
!                        indicate if the convective scheme included radiative
!                        feedbacks.  Added a flag (IFKFRADEXTRAS) for extra
!                        variables available with KF convective scheme with
!                        radiative feedbacks.  (T. Spero)
!           24 Feb 2020  Adapted for FV3GFSv16 at NOAA-ARL (P. C. Campbell)
!-------------------------------------------------------------------------------

  USE metinfo
  USE date_pack
  USE mcipparm
  USE coord
  USE files
  USE netcdf_io
  USE const, ONLY: pi180
  USE netcdf

  IMPLICIT NONE

  INTEGER     ,       INTENT(IN)    :: cdfid, cdfid2
  INTEGER                           :: cdfidg
  INTEGER                           :: cdfid_vgvf
  INTEGER                           :: cdfid_vlai
  REAL,               INTENT(OUT)   :: ctmlays     ( maxlays )
  REAL                              :: phalf_lays  ( maxlays )
  REAL                              :: pfull_lays  ( maxlays+1 )
  CHARACTER(LEN=32)                 :: date_init
  INTEGER                           :: dimid
  INTEGER                           :: dimids     ( nf90_max_var_dims )
  REAL,               ALLOCATABLE   :: dum1d      ( : )
  REAL,               ALLOCATABLE   :: dum2d      ( : , : )
  INTEGER,            ALLOCATABLE   :: dum2d_i    ( : , : )
  INTEGER                           :: dx
  INTEGER                           :: dy
  REAL                              :: fac
  CHARACTER(LEN=256)                :: fl
  CHARACTER(LEN=256)                :: fl2
  CHARACTER(LEN=256)                :: flg
  CHARACTER(LEN=256)                :: geofile
  CHARACTER(LEN=256)                :: viirsgvf
  CHARACTER(LEN=256)                :: viirslai
  INTEGER                           :: icloud_cu
  INTEGER                           :: id_data
  INTEGER                           :: idtsec
  LOGICAL                           :: ifgeo
  LOGICAL                           :: ifisltyp
  LOGICAL                           :: ifra
  LOGICAL                           :: ifrs
  LOGICAL                           :: ifsmois
  LOGICAL                           :: iftslb
  LOGICAL                           :: ifu10m
  LOGICAL                           :: ifv10m
  INTEGER                           :: it,iv
  INTEGER                           :: ival
  INTEGER                           :: lent
  INTEGER                           :: n_times
  INTEGER                           :: nxm
  INTEGER                           :: nym
  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'SETUP_FV3'
  INTEGER                           :: rcode, rcode2,imax2,jmax2
  REAL                              :: rval
  REAL,  ALLOCATABLE                :: times      ( : )
  INTEGER                           :: varid
  CHARACTER(LEN=80)                 :: fv3version
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f6000 = "(/, 1x, &
    & '- SUBROUTINE SETUP_FV3 - READING FV3 HEADER')"
  CHARACTER(LEN=256), PARAMETER :: f6100 = "(3x, &
    & 'FV3 GRID DIMENSIONS (X,Y,Z) ', i4, 1x, i4, 1x, i3, //)"

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   MISMATCH IN DX AND DY', &
    & /, 1x, '***   DX, DY = ', 2(f7.2), &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNKNOWN LAND USE CLASSIFICATION', &
    & /, 1x, '***   FIRST THREE LETTERS = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9225 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   QCLOUD NOT FOUND IN FV3 OUTPUT...STOPPING', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9250 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ETA/FERRIER SCHEME IS NOT SUPPORTED IN CMAQ', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9275 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   FOUND QCLOUD BUT NOT QRAIN...STOPPING', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   NQSPECIES SET AT 3',&
    & /, 1x, '***   MCIP NEEDS TO BE MODIFIED FOR THIS CASE', &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
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
    & /, 1x, '***   ONLY FOUND ONE FILE WITH ONE TIME PERIOD', &
    & /, 1x, '***   SETTING OUTPUT FREQUENCY TO 1 MINUTE', &
    & /, 1x, 70('*'))" 

  CHARACTER(LEN=256), PARAMETER :: f9550 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   NEED PRECIPITATION ACCUMULATION IN FV3 TO MATCH', &
    & /, 1x, '***   MCIP OUTPUT INTERVAL', &
    & /, 1x, '***   PREC_ACC_DT from FV3: ', i4, &
    & /, 1x, '***   INTVL from MCIP: ', i4, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9600 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING FV3 NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9610 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING VIIRS NETCDF FILE', &
    & /, 1x, '***   WILL NOT USE VIIRS GVF' &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9700 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING FV3 NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9800 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   DID NOT FIND FRACTIONAL LAND USE IN FV3 output', &
    & /, 1x, '***   AND DID NOT FIND GEOGRID FILE' &
    & /, 1x, '***   -- WILL NOT USE FRACTIONAL LAND USE DATA' &
    & /, 1x, 70('*'))"
  
  CHARACTER(LEN=256), PARAMETER :: f9900 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   DID NOT FIND LEAF AREA INDEX IN FV3 output', &
    & /, 1x, '***   AND DID NOT FIND GEOGRID FILE' &
    & /, 1x, '***   -- WILL NOT USE LEAF AREA INDEX DATA' &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9910 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM VIIRS FILE', &
    & /, 1x, '***   WILL NOT USE VIIRS GVF' &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Extract NX, NY, and NZ.
!-------------------------------------------------------------------------------
  WRITE (*,f6000)

  fl = file_mm(1)

  rcode = nf90_get_att (cdfid, nf90_global, 'im', met_nx)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'WEST-EAST_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  
  rcode2 = nf90_get_att (cdfid2, nf90_global, 'im', imax2)
  IF ( rcode2 /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SFC WEST-EAST_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_get_att (cdfid, nf90_global, 'jm',  &
                        met_ny)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SOUTH-NORTH_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode2 = nf90_get_att (cdfid2, nf90_global,'jm',jmax2)
  IF ( rcode2 /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'SFC SOUTH-NORTH_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  
  if(met_nx.ne.imax2.or.met_ny.ne.jmax2) then
   print*,'inconsistent i,j dimension between ATM and SFC files ',&
     met_nx,imax2,met_ny,jmax2
   call graceful_stop(pname)
  endif   

  rcode = nf90_inq_dimid (cdfid, 'phalf', dimid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'ID for phalf',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  rcode = nf90_inquire_dimension (cdfid, dimid, len=ival)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'BOTTOM-TOP_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE

! Set met_nz
  met_nz = ival-1

  ENDIF

  WRITE (*,f6100) met_nx, met_ny, met_nz

  met_rictr_dot = FLOAT(met_nx - 1) / 2.0 + 1.0
  met_rjctr_dot = FLOAT(met_ny - 1) / 2.0 + 1.0

!-------------------------------------------------------------------------------
! Read FV3 Pressure layers.
!-------------------------------------------------------------------------------
  rcode = nf90_inq_varid (cdfid, 'phalf', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'phalf',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE
  ALLOCATE ( dum1d ( met_nz+1 ) ) 
  rcode = nf90_get_var (cdfid, varid, dum1d)

  pfull_lays(1:met_nz+1) = dum1d(met_nz+1:1:-1) 

  DEALLOCATE (dum1d)
  ENDIF
  
  rcode = nf90_inq_varid (cdfid, 'pfull', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'pfull',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ELSE
  ALLOCATE ( dum1d ( met_nz ) )
  rcode = nf90_get_var (cdfid, varid, dum1d)

  phalf_lays(1:met_nz) = dum1d(met_nz:1:-1)

  DEALLOCATE (dum1d)
  ENDIF
!
!---Read FV3 Lat/lon

  allocate(fv3lat(met_ny),fv3lon(met_nx))
  
  CALL get_var_1d_double_cdf (cdfid, 'grid_yt', fv3lat, 1, rcode)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'grid_yt', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF
  CALL get_var_1d_double_cdf (cdfid, 'grid_xt', fv3lon, 1, rcode)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'grid_xt', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

!-------------------------------------------------------------------------------
! Extract domain attributes.
!-------------------------------------------------------------------------------

   rcode = nf90_get_att (cdfid, nf90_global, 'source', fv3version)
   IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'source', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

!FV3 the user sets the input dx/dy resolutions in namelist
   dx=dx_in
   dy=dy_in
!
  IF (dx == dy) THEN
    met_resoln = dx
  ELSE
    WRITE (*,f9000) TRIM(pname), dx, dy
    CALL graceful_stop (pname)
  ENDIF

  met_nxcoarse = met_nx 
  met_nycoarse = met_ny
  met_gratio   = 1
  met_x_11     = 1
  met_y_11     = 1

!     FV3 Gaussian Global Grid
      met_mapproj    = 4                      ! FV3 Gaussian Map Projection      
!      met_proj_clon  = 0.0
!      met_p_alp_d  = 0.0                      ! lat of coord origin [deg]
!      met_p_bet_d  = 0.0                      ! (not used)
!      met_p_gam_d  = met_proj_clon            ! lon of coord origin [deg]
!      met_cone_fac = 0.0                      ! cone factor
!      met_ref_lat  = -999.0                   ! not used
      
!      met_cen_lat  = met_cen_lat_in          ! set from namelist
!      met_cen_lon  = met_cen_lon_in          ! set from namelist                  

!-------------------------------------------------------------------------------
! Extract model run options.
!-------------------------------------------------------------------------------
    met_ns=4  

  ! Determine if NOAH Mosaic was run and created the correct output fields.
  ! Note that this code is temporarily modified later in this subroutine to
  ! toggle IFMOSAIC to FALSE if the fractional land use arrays are also not
  ! available.  That change is temporary until a subroutine is added to
  ! reconstruct the fractional land use rank if it is missing.

  ! FV3 does not have a Noah Mosaic option, so turn off
  ifmosaic = .FALSE.

  ! Define number of land use categories.

  rcode = nf90_inq_varid (cdfid2, 'vtype', varid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'Land Use Type',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
 
  ELSE
  nxm = met_nx
  nym = met_ny 
  ALLOCATE ( dum2d_i ( nxm, nym ) )
  rcode = nf90_get_var (cdfid2, varid, dum2d_i)
  nummetlu=MAXVAL(dum2d_i)
  DEALLOCATE (dum2d_i)
  ENDIF

! FV3 only  has MODIS 20-category ('MOD'):  Harcoded Ice, water, urban, lake         
! FV3 doesn't have lake, category 21, but if using geofile, we add it from MODIS
! input
          met_lu_src = 'MOD'
          met_lu_water = 17 !MODIS IGBP water = 17...
          !But, FV3 sets water category =0, and need extra met_lu_water_fv3
          !variable for calulating PURB and LWMASK later
          met_lu_water_fv3 = 0
          met_lu_ice   = 15
          met_lu_urban = 13
          IF ( file_geo(1:7) == " " ) THEN
           met_lu_lake =  -1
          ELSE
           nummetlu = nummetlu +1 !added lake category from geofile
           met_lu_lake =  21
          ENDIF

   met_lw_rad=0
   met_sw_rad=0
   met_cumulus=0   
   met_expl_moist=0
   met_pbl=0
   met_sfc_lay=0
   met_soil_lsm=0
   met_urban_phys = 0
   ifrcurb = .FALSE.
   met_shal_cu = 0

  met_snow_opt = 1       ! not used for FV3 yet, but set conistent with original MCIP
  met_rain_bucket = -1.0 ! FV3 includes both bucket and prate, and leave bucket turned off
  met_pcp_incr = 1       ! FV3 has time step average prate, so leave pcp_inc turned on

  ! Currently there are no radiative feedbacks in FV3 accompaning the convective 
  ! parameterization scheme, turn off.
  ifcuradfdbk = .FALSE.

!-------------------------------------------------------------------------------
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

  met_startdate =  date_init(13:32) // '.0000'

  ! FV3 always turned off for "met restart" in MCIP
  met_restart = 0
 
  rcode = nf90_inquire_variable (cdfid, varid, dimids=dimids)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9420) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inquire_dimension (cdfid, dimids(1), len=n_times)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9430) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  IF ( ALLOCATED ( times ) ) DEALLOCATE ( times )
  ALLOCATE ( times ( n_times ) )
  rcode = nf90_get_var (cdfid, varid, times)
!                        start=(/1/), count=(/n_times/))
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'time', TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

   met_tapfrq = (times(1))*60  ! convert hrs --> min
   IF ( met_model == 3 ) THEN !IF FV3 and first 000 forecast hour, assume hourly input = 60 min
    IF (met_tapfrq == 0.0) THEN
     met_tapfrq = 60.0 !minutes
    ENDIF
   ENDIF
   
!-------------------------------------------------------------------------------
! If layer structure was not defined in user namelist, use use FV3 pressure
! to calculate ctmlays.
!-------------------------------------------------------------------------------
  IF ( needlayers ) THEN
     nlays=met_nz
     ctmlays(:) = 0.0 !Initialize
     ctmlays = (pfull_lays(1:nlays+1) - pfull_lays(nlays+1)) / ((pfull_lays(1)) - phalf_lays(nlays+1))
  ENDIF

!-------------------------------------------------------------------------------
! Set variables for non-hydrostatic base state.  There is no option for
! hydrostatic run in FV3.  The base state variables are not currently output
! , so fill in "default" values from FV3 namelist.
!
! Note:  In FV3v2.2 NCAR changed the way "real" scalars (e.g., P_TOP) are
!        stored in the FV3 I/O API.
!-------------------------------------------------------------------------------

  met_ptop  = pfull_lays(met_nz+1)*100.0
  met_p00   = 100000.0 ! base state sea-level pressure [Pa]
  met_ts0   =    290.0 ! base state sea-level temperature [K]
  met_tlp   =     50.0 ! base state lapse rate d(T)/d(ln P) from 1000 to 300 mb
  met_tiso  = fillreal ! base state stratospheric isothermal T [K]  ! not used

!-------------------------------------------------------------------------------
! Determine FV3 release.
!-------------------------------------------------------------------------------

  met_release = ' FV3  '

!-------------------------------------------------------------------------------
! Determine FDDA options.
!-------------------------------------------------------------------------------
  met_fdda_3dan = 0
  met_fdda_gv3d = -1.0 
  met_fdda_gt3d = -1.0
  met_fdda_gq3d = -1.0
  met_fdda_gph3d = -1.0
  met_fdda_sfan  =  0  ! sfc analysis nudging not in FV3
  met_fdda_gvsfc = -1.0
  met_fdda_gtsfc = -1.0
  met_fdda_gqsfc = -1.0
  met_fdda_obs = 0       ! not implemented in FV3
  met_fdda_giv = -1.0    ! not implemented in FV3
  met_fdda_git = -1.0    ! not implemented in FV3
  met_fdda_giq = -1.0    ! not implemented in FV3

!-------------------------------------------------------------------------------
! Determine whether or not fractional land use is available in the output.
! Set the flag appropriately.
!-------------------------------------------------------------------------------
  ! If fractional LANDUSEfF becomes available in FV3 or, if user provides a
  ! file_geo with LANDUSEF,it does not stop model if not available (needed for
  ! tiled deposition with STAGE.
  rcode2 = nf90_inq_varid (cdfid2, 'LANDUSEF', varid)
  IF ( rcode2 == nf90_noerr ) THEN
    iflufrc    = .TRUE.   ! fractional land use is available
    ifluwrfout = .TRUE.   ! fractional land use is located in FV3 history file
  ELSE
    ifluwrfout = .FALSE.  ! fractional land use is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9800) TRIM(pname)
      iflufrc = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'LANDUSEF', varid)
      IF ( rcode == nf90_noerr ) THEN
        iflufrc = .TRUE.  ! fractional land use is in the file
      ELSE
        iflufrc = .FALSE. ! fractional land use is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF  

  ! For now, require LANDUSEF2 and MOSAIC_CAT_INDEX to process NOAH Mosaic.
  ! IFMOSAIC is toggled to FALSE here if either field is missing.
  ! Can add a subroutine later to reconstruct those fields from LANDUSEF if
  ! LANDUSEF2 and/or MOSAIC_CAT_INDEX is missing.

  IF ( ifmosaic ) THEN
    rcode = nf90_inq_varid (cdfid, 'LANDUSEF2', varid)
    IF ( rcode == nf90_noerr ) THEN
      rcode = nf90_inq_varid (cdfid, 'MOSAIC_CAT_INDEX', varid)
      IF ( rcode == nf90_noerr ) THEN
        iflu2wrfout = .TRUE.   ! lookup for LANDUSEF2 is in FV3 history
      ELSE
        iflu2wrfout = .FALSE.
        ifmosaic    = .FALSE.
      ENDIF
    ELSE
      iflu2wrfout = .FALSE.  ! frac land use 2 is not available in FV3 history
      ifmosaic    = .FALSE.
    ENDIF
  ELSE
    iflu2wrfout = .FALSE.
  ENDIF
!-------------------------------------------------------------------------------
! Determine whether or not the 2-m temperature, the 2-m mixing ratio, the
! 10-m wind components, and the turbulent kinetic energy are in the output,
! and set the flags appropriately.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid2, 'tmp2m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ift2m = .TRUE.  ! 2-m temperature is in the file
  ELSE
    ift2m = .FALSE. ! 2-m temperature is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'spfh2m', varid)
  IF ( rcode == nf90_noerr ) THEN
    IF ( met_pbl == 1 ) THEN  ! YSU PBL scheme
      ifq2m = .FALSE. ! do not use Q2 from YSU PBL; occasional winter negatives
    ELSE
      ifq2m = .TRUE.  ! 2-m mixing ratio is in the file
    ENDIF
  ELSE
    ifq2m = .FALSE. ! 2-m mixing ratio is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'ugrd10m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifu10m = .TRUE.  ! 10-m u-component wind is in the file
  ELSE
    ifu10m = .FALSE. ! 10-m u-component wind is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'vgrd10m', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifv10m = .TRUE.  ! 10-m v-component wind is in the file
  ELSE
    ifv10m = .FALSE. ! 10-m v-component wind is not in the file
  ENDIF

  IF ( ( ifu10m ) .AND. ( ifv10m ) ) THEN
    ifw10m = .TRUE.
  ELSE
    ifw10m = .FALSE.
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'TKE_MYJ', varid)  !Left as TKE_MYJ, not in FV3 = False
  IF ( rcode == nf90_noerr ) THEN
    IF ( met_pbl == 2 ) THEN  ! Mellor-Yamada-Janjic (Eta)
      iftke  = .TRUE.  ! turbulent kinetic energy is in the file
      iftkef = .FALSE. ! TKE is not on full-levels; it is on half-layers
    ELSE
      iftke  = .FALSE. ! turbulent kinetic energy is not in the file
      iftkef = .FALSE.
    ENDIF
  ELSE
    iftke  = .FALSE. ! turbulent kinetic energy is not in the file
    iftkef = .FALSE.
  ENDIF
!-------------------------------------------------------------------------------
! Determine whether or not some surface variables are in the output, and set
! the flags appropriately.
!-------------------------------------------------------------------------------
 IF ( ( iffengsha_dust ) ) THEN  !User is trying to use Fengsha Windblown Dust in CMAQ

  rcode2 = nf90_inq_varid (cdfid2, 'CLAY_FRAC', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    ifclayf     = .TRUE.   ! clay fraction is in the file
    ifclayfwrfout = .TRUE.   ! clay fraction is not in the file
  ELSE
    ifclayfwrfout = .FALSE.  ! clay fraction is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      ifclayf = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'CLAY_FRAC', varid)
      IF ( rcode == nf90_noerr ) THEN
        ifclayf = .TRUE.  ! clay fraction is in the file
      ELSE
        ifclayf = .FALSE. ! clay fraction is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF

  rcode2 = nf90_inq_varid (cdfid2, 'SAND_FRAC', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    ifsandf     = .TRUE.   ! sand fraction is in the file
    ifsandfwrfout = .TRUE.   ! sand fraction is not in the file
  ELSE
    ifsandfwrfout = .FALSE.  ! sand fraction is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      ifsandf = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'SAND_FRAC', varid)
      IF ( rcode == nf90_noerr ) THEN
        ifsandf = .TRUE.  ! sand fraction is in the file
      ELSE
        ifsandf = .FALSE. ! sand fraction is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF
  
  rcode2 = nf90_inq_varid (cdfid2, 'DRAG_PART', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    ifdrag     = .TRUE.   ! drag partition is in the file
    ifdragwrfout = .TRUE.   ! drag partition is not in the file
  ELSE
    ifdragwrfout = .FALSE.  ! drag partition is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      ifdrag = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'DRAG_PART', varid)
      IF ( rcode == nf90_noerr ) THEN
        ifdrag = .TRUE.  ! drag partition is in the file
      ELSE
        ifdrag = .FALSE. ! drag partition is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF

  rcode2 = nf90_inq_varid (cdfid2, 'SSM', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    ifssm     = .TRUE.   ! sediment supply map is in the file
    ifssmwrfout = .TRUE.   ! sediment supply map is not in the file
  ELSE
    ifssmwrfout = .FALSE.  ! sediment supply map is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      ifssm = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'SSM', varid)
      IF ( rcode == nf90_noerr ) THEN
        ifssm = .TRUE.  ! sediment supply map is in the file
      ELSE
        ifssm = .FALSE. ! sediment supply map is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF

  rcode2 = nf90_inq_varid (cdfid2, 'UTHRES', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    ifuthr     = .TRUE.   ! threshold velocity is in the file
    ifuthrwrfout = .TRUE.   ! threshold velocity is not in the file
  ELSE
    ifuthrwrfout = .FALSE.  ! threshold velocity is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      ifuthr = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'UTHRES', varid)
      IF ( rcode == nf90_noerr ) THEN
        ifuthr = .TRUE.  ! threshold velocity is in the file
      ELSE
        ifuthr = .FALSE. ! threshold velocity is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF
  ENDIF
 ENDIF !Fengsha WB dust variables


 IF ( ( ifviirs_lai ) ) THEN  !Check if user is using VIIRS LAI instead of FV3
   viirslai = TRIM( file_viirs_lai )
   flg = viirslai
   rcode = nf90_open (flg, nf90_nowrite, cdfid_vlai)
   IF ( rcode == nf90_noerr ) THEN !file is present
      rcode = nf90_inq_varid (cdfid_vlai, 'lai', varid)
      IF ( rcode == nf90_noerr ) THEN ! LAI is found in the VIIRS file
        iflai = .TRUE.
        iflai_viirs = .TRUE.
       !---Next Check Latitude dimension and if present read VIIRS Lat
       rcode = nf90_inq_dimid (cdfid_vlai, 'latitude', dimid)
       IF ( rcode == nf90_noerr ) THEN !latitude dimension is there
        rcode = nf90_inquire_dimension (cdfid_vlai, dimid, len=ival)
        met_ny_viirs = ival
        if(.not.allocated(viirslat_lai)) allocate(viirslat_lai(met_ny_viirs))
        CALL get_var_1d_double_cdf (cdfid_vlai, 'latitude', viirslat_lai, 1, rcode)
!        viirslat=viirslat(met_ny_viirs:1:-1)
       ELSE !latitude dimension is not there
         WRITE (*,f9910) TRIM(pname), 'latitude', TRIM(nf90_strerror(rcode))
         iflai = .FALSE.
         iflai_viirs = .FALSE.
         rcode = nf90_inq_varid (cdfid2, 'lai', varid)
         IF ( rcode == nf90_noerr ) THEN
          iflai = .TRUE.  ! LAI is in the file
          iflaiwrfout = .TRUE.   ! leaf area index is in the file
         ELSE
          iflaiwrfout = .FALSE. ! LAI is not in the file
         ENDIF
!         CALL graceful_stop (pname) !NACC will not stop --> default to LAI/LU table calculation.
       ENDIF
       !---Next Check Longitude dimension and if present read VIIRS Lon
       rcode = nf90_inq_dimid (cdfid_vlai, 'longitude', dimid)
       IF ( rcode == nf90_noerr ) THEN !longitude dimension is there
        rcode = nf90_inquire_dimension (cdfid_vlai, dimid, len=ival)
        met_nx_viirs = ival
        if(.not.allocated(viirslon_lai)) allocate(viirslon_lai(met_nx_viirs))
        CALL get_var_1d_double_cdf (cdfid_vlai, 'longitude', viirslon_lai, 1, rcode)
        !conform to fv3 longitude values, which is 0-->360
        do iv=1,met_nx_viirs
         if(viirslon_lai(iv).lt.0) viirslon_lai(iv)=viirslon_lai(iv)+360
        enddo
       ELSE !longitude dimension is not there
         WRITE (*,f9910) TRIM(pname), 'longitude', TRIM(nf90_strerror(rcode))
         iflai = .FALSE.
         iflai_viirs = .FALSE.
         rcode = nf90_inq_varid (cdfid2, 'lai', varid)
         IF ( rcode == nf90_noerr ) THEN
          iflai = .TRUE.  ! LAI is in the file
          iflaiwrfout = .TRUE.   ! leaf area index is in the file
         ELSE
          iflaiwrfout = .FALSE. ! LAI is not in the file
         ENDIF
!         CALL graceful_stop (pname)  !NACC will not stop --> default to LAI/LU table calculation.
       ENDIF

      ELSE !! LAI is not in the VIIRS file-->default to check FV3 LAI if available
        iflai_viirs = .FALSE.
        WRITE (*,f9610) TRIM(pname), TRIM(flg)
        rcode = nf90_inq_varid (cdfid2, 'lai', varid)
        IF ( rcode == nf90_noerr ) THEN
         iflai = .TRUE.  ! LAI is in the file
         iflaiwrfout = .TRUE.   ! leaf area index is in the file
        ELSE
         iflaiwrfout= .FALSE. ! LAI is not in the file
        ENDIF
!        CALL graceful_stop (pname)   !NACC will not stop --> default to LAI/LU table calculation.
      ENDIF
   ELSE !Can't find/open VIIRS LAI file, default to check FV3 LAI if available
      iflai_viirs = .FALSE.
      WRITE (*,f9610) TRIM(pname), TRIM(flg)
      rcode = nf90_inq_varid (cdfid2, 'lai', varid)
      IF ( rcode == nf90_noerr ) THEN
       iflai = .TRUE.  ! LAI is in the file
       iflaiwrfout = .TRUE.   ! leaf area index is in the file
      ELSE
       iflaiwrfout = .FALSE. ! LAI is not in the file
      ENDIF
!      CALL graceful_stop (pname)  !NACC will not stop --> default to LAI/LU table calculation.
   ENDIF
   rcode = nf90_close (cdfid_vlai)

 ELSE !User is not using VIIRS LAI (Default)-->Check LAI is in either FV3 or Geofile
   iflai_viirs = .FALSE.
   rcode2 = nf90_inq_varid (cdfid2, 'LAI', varid) !not in FV3GFSv16
  IF ( rcode2 == nf90_noerr ) THEN
    iflai    = .TRUE.   ! leaf area index is in the file
    iflaiwrfout = .TRUE.   ! leaf area index is in the file
  ELSE
    iflaiwrfout = .FALSE.  ! leaf area index is not available in FV3 history
    geofile = TRIM( file_geo )
    INQUIRE ( FILE=geofile, EXIST=ifgeo )
    IF ( .NOT. ifgeo ) THEN
      WRITE (*,f9900) TRIM(pname)
      iflai = .FALSE.
    ELSE
      flg = file_geo
       rcode = nf90_open (flg, nf90_nowrite, cdfidg)

      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9600) TRIM(pname), TRIM(flg)
        CALL graceful_stop (pname)
      ENDIF
      rcode = nf90_inq_varid (cdfidg, 'LAI', varid)
      IF ( rcode == nf90_noerr ) THEN
        iflai = .TRUE.  ! leaf area index is in the file
      ELSE
        iflai = .FALSE. ! leaf area index is not in the file
      ENDIF
      rcode = nf90_close (cdfidg)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9700)  TRIM(pname),TRIM(flg)
!        CALL graceful_stop (pname)  !NACC will not stop --> default to LAI/LU table calculation.
      ENDIF
    ENDIF
  ENDIF

 ENDIF !End VIIRS/FV3/Geofile LAI Checks

  rcode = nf90_inq_varid (cdfid2, 'RMOL', varid) !not in FV3GFSv16
  IF ( rcode == nf90_noerr ) THEN
    ifmol = .TRUE.  ! (inverse) Monin-Obukhov length is in the file
  ELSE
    ifmol = .FALSE. ! (inverse) Monin-Obukhov length is not in the file
  ENDIF

    ifmolpx   = .FALSE.
    ifpxwrf41 = .FALSE.

  rcode = nf90_inq_varid (cdfid2, 'acond', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifra = .TRUE.  ! aerodynamic resistance is in the file
  ELSE
    ifra = .FALSE. ! aerodynamic resistance is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'RS', varid)  !not in FV3GFSv16??
  IF ( rcode == nf90_noerr ) THEN
    ifrs = .TRUE.  ! stomatal resistance is in the file
  ELSE
    ifrs = .FALSE. ! stomatal resistance is not in the file
  ENDIF

  IF ( ( ifra ) .AND. ( ifrs ) ) THEN
    ifresist = .TRUE.
  ELSE
    ifresist = .FALSE.
  ENDIF

  IF ( ( ifviirs_gvf ) ) THEN  !User is using VIIRS GVF for vegetation fraction instead of FV3
   viirsgvf = TRIM( file_viirs_gvf )
   flg = viirsgvf
   rcode = nf90_open (flg, nf90_nowrite, cdfid_vgvf)
   IF ( rcode == nf90_noerr ) THEN !file is present
      rcode = nf90_inq_varid (cdfid_vgvf, 'VEG_surface', varid)
      IF ( rcode == nf90_noerr ) THEN ! vegetation fraction is in the VIIRS file
        ifveg_viirs = .TRUE.  
       !---Next Check Latitude dimension and if present read VIIRS Lat
       rcode = nf90_inq_dimid (cdfid_vgvf, 'latitude', dimid)
       IF ( rcode == nf90_noerr ) THEN !latitude dimension is there
        rcode = nf90_inquire_dimension (cdfid_vgvf, dimid, len=ival)
        met_ny_viirs = ival
        if(.not.allocated(viirslat_gvf)) allocate(viirslat_gvf(met_ny_viirs))
        CALL get_var_1d_double_cdf (cdfid_vgvf, 'latitude', viirslat_gvf, 1, rcode)
!        viirslat=viirslat(met_ny_viirs:1:-1)
       ELSE !latitude dimension is not there
         WRITE (*,f9910) TRIM(pname), 'latitude', TRIM(nf90_strerror(rcode))
         ifveg_viirs = .FALSE.
         rcode = nf90_inq_varid (cdfid2, 'veg', varid)
         IF ( rcode == nf90_noerr ) THEN
          ifveg = .TRUE.  ! vegetation fraction is in the file
         ELSE
          ifveg = .FALSE. ! vegetation fraction is not in the file
         ENDIF
!         CALL graceful_stop (pname)
       ENDIF
       !---Next Check Longitude dimension and if present read VIIRS Lon
       rcode = nf90_inq_dimid (cdfid_vgvf, 'longitude', dimid)
       IF ( rcode == nf90_noerr ) THEN !longitude dimension is there
        rcode = nf90_inquire_dimension (cdfid_vgvf, dimid, len=ival)
        met_nx_viirs = ival
        if(.not.allocated(viirslon_gvf)) allocate(viirslon_gvf(met_nx_viirs))
        CALL get_var_1d_double_cdf (cdfid_vgvf, 'longitude', viirslon_gvf, 1, rcode)
        !conform to fv3 longitude values, which is 0-->360
        do iv=1,met_nx_viirs
         if(viirslon_gvf(iv).lt.0) viirslon_gvf(iv)=viirslon_gvf(iv)+360
        enddo
       ELSE !longitude dimension is not there
         WRITE (*,f9910) TRIM(pname), 'longitude', TRIM(nf90_strerror(rcode))
         ifveg_viirs = .FALSE.
         rcode = nf90_inq_varid (cdfid2, 'veg', varid)
         IF ( rcode == nf90_noerr ) THEN
          ifveg = .TRUE.  ! vegetation fraction is in the file
         ELSE
          ifveg = .FALSE. ! vegetation fraction is not in the file
         ENDIF
!         CALL graceful_stop (pname)
       ENDIF

      ELSE !! vegetation fraction is not in the VIIRS file-->default to check FV3 veg if available
        ifveg_viirs = .FALSE.
        WRITE (*,f9610) TRIM(pname), TRIM(flg)
        rcode = nf90_inq_varid (cdfid2, 'veg', varid)
        IF ( rcode == nf90_noerr ) THEN
         ifveg = .TRUE.  ! vegetation fraction is in the file
        ELSE
         ifveg = .FALSE. ! vegetation fraction is not in the file
        ENDIF
!        CALL graceful_stop (pname)   !doesnt stop NACC
      ENDIF
   ELSE !Can't find/open VIIRS GVF file, default to check FV3 veg if available
      ifveg_viirs = .FALSE.
      WRITE (*,f9610) TRIM(pname), TRIM(flg)
      rcode = nf90_inq_varid (cdfid2, 'veg', varid)
      IF ( rcode == nf90_noerr ) THEN
       ifveg = .TRUE.  ! vegetation fraction is in the file
      ELSE
       ifveg = .FALSE. ! vegetation fraction is not in the file
      ENDIF
!      CALL graceful_stop (pname)  !doesnt stop NACC
   ENDIF
   rcode = nf90_close (cdfid_vgvf)

  ELSE !User is not using VIIRS (Default)-->check to make sure vegetation fraction is in FV3
   ifveg_viirs = .FALSE.
   rcode = nf90_inq_varid (cdfid2, 'veg', varid) 
   IF ( rcode == nf90_noerr ) THEN
     ifveg = .TRUE.  ! vegetation fraction is in the file
   ELSE
     ifveg = .FALSE. ! vegetation fraction is not in the file
   ENDIF

  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'cnwat', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifwr = .TRUE.  ! canopy wetness is in the file
  ELSE
    ifwr = .FALSE. ! canopy wetness is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'soill1', varid) 
  IF ( rcode == nf90_noerr ) THEN
    ifsmois = .TRUE.  ! soil moisture is in the file
  ELSE
    ifsmois = .FALSE. ! soil moisture is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'soilt1', varid)
  IF ( rcode == nf90_noerr ) THEN
    iftslb = .TRUE.  ! soil temperature is in the file
  ELSE
    iftslb = .FALSE. ! soil temperature is not in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'sotyp', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifisltyp = .TRUE.  ! soil type is in the file
  ELSE
    ifisltyp = .FALSE. ! soil type is not in the file
  ENDIF

  If ( ( ifsmois ) .AND. ( iftslb ) .AND. ( ifisltyp ) ) THEN
    ifsoil = .TRUE.
  ELSE
    ifsoil = .FALSE.
  ENDIF

  rcode = nf90_inq_varid (cdfid2, 'sfcr', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifznt = .TRUE.  ! roughness length is in the file
  ELSE
    ifznt = .FALSE. ! roughness length is not in the file
  ENDIF

!KF-Radiative Feedback is not in FV3

    ifkfradextras = .FALSE.

!-------------------------------------------------------------------------------
! Determine the number of 3D cloud moisture species.  Assume that cloud water
! mixing ratio and rain water mixing ratio will occur together.  Also assume
! that cloud ice mixing ratio and cloud snow mixing ratio will occur together,
! but check for availability.  Check for graupel, as well.
! Note:  In FV3v2.1.2 and prior, the Eta/Ferrier microphysics scheme only
! outputs QCLOUD which represents total condensate, not cloud water mixing
! ratio.  CMAQv4.6 and prior cannot handle this field, so MCIP will stop in
! this case.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid, 'clwmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = 1  ! QCLOUD is in the file
  ELSE  ! need hydrometeor fields for CMAQ
    WRITE (*,f9225) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'rwmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QRAIN is in the file
  ELSE
    IF ( met_expl_moist == 5 ) THEN  ! Eta/Ferrier scheme
      WRITE (*,f9250) TRIM(pname)
      CALL graceful_stop (pname)
    ELSE
      WRITE (*,f9275) TRIM(pname)
      CALL graceful_stop (pname)
    ENDIF
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'icmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QICE is in the file
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'snmr', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QSNOW is in the file
  ENDIF

  IF ( nqspecies == 3 ) THEN  ! not set up for QI w/o QS or vice versa
    WRITE (*,f9300) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_varid (cdfid, 'grle', varid)
  IF ( rcode == nf90_noerr ) THEN
    nqspecies = nqspecies + 1  ! QGRAUP is in the file
  ENDIF

  IF ( nqspecies == 3 ) THEN  ! not set up for QG without QI and QS
    WRITE (*,f9300) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

!-------------------------------------------------------------------------------
! Determine whether 3D resolved cloud fraction is part of FV3 output.  If
! Kain-Fritsch scheme with radiative feedbacks to subgrid clouds is used (new
! in FV3v3.6) or if MSKF is used (new in FV3v3.7) in FV3, then the 3D cloud
! fraction includes both resolved and subgrid clouds.
!-------------------------------------------------------------------------------

  rcode = nf90_inq_varid (cdfid, 'cld_amt', varid)
  IF ( rcode == nf90_noerr ) THEN
    ifcld3d = .TRUE.  ! 3D resolved cloud fraction is in the file
  ELSE
    ifcld3d = .FALSE. ! 3D cloud fraction is not if the file
  ENDIF
!-------------------------------------------------------------------------------
! Determine if the hybrid vertical coordinate has been used in FV3.  It is
! available as of FV3v3.9.
!-------------------------------------------------------------------------------
   met_hybrid = 2   ! FV3 employs a hybrid vertical coordinate as default = 2 (mimic WRF hybrid flag)

END SUBROUTINE setup_fv3
