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

SUBROUTINE setup (ctmlays, itimestep)

!-------------------------------------------------------------------------------
! Name:     Set Up the Input Meteorology Domain Attributes
! Purpose:  Establishes bounds for MM5 or WRF post-processing.
! Revised:  10 Sep 2001  Original version.  (T. Otte)
!           07 Jan 2002  Changed file name to explicit file rather than
!                        Fortran unit to improve portability.  (S. Howard
!                        and T. Otte)
!           09 Jan 2002  Changed calls to "abort" to calls to "m3exit" for
!                        graceful shut-down of I/O API files.  (T. Otte)
!           26 May 2005  Added WRF capability.  Changed routine name from
!                        SETUPMM5 to SETUP to make code more general.  (T. Otte)
!           09 Apr 2007  Removed option to handle MM5v2-formatted data.
!                        (T. Otte)
!           22 Apr 2008  Set WRF DYN_OPT to 2 (mass core) for WRFv3 and
!                        beyond because support for other cores within WRF-ARW
!                        was discontinued in WRFv3.  (T. Otte)
!           17 Mar 2010  Changed all calls to netCDF routines to use the
!                        Fortran interface rather than the C interface.
!                        Rearranged subroutine to improve efficiency.  Removed
!                        dependency on module WRF_NETCDF.  Improved clarity
!                        in some error-handling messages.  Added CDFID to the
!                        argument list for subroutine SETUP_WRFEM.  (T. Otte)
!           31 Aug 2011  Changed name of module FILE to FILES to avoid conflict
!                        with F90 protected intrinsic.  Updated netCDF commands
!                        to F90, and improved error handling.  Changed F77
!                        character declarations to F90 standard.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           14 Sep 2018  Removed support for MM5v3 input.  (T. Spero)
!           15 Nov 2018  Allow WRFv4.0 input to be used.  (T. Spero)
!           18 Nov 2019  Modified for FV3GFS Capability. (P. C. Campbell)
!-------------------------------------------------------------------------------

  USE mcipparm
  USE metinfo
  USE files
  USE netcdf

  IMPLICIT NONE

  INTEGER                           :: cdfid, cdfid2
  REAL,               INTENT(INOUT) :: ctmlays   ( maxlays )
  CHARACTER(LEN=19)                 :: gridtype
  INTEGER                           :: istat
  CHARACTER(LEN=16),  PARAMETER     :: pname     = 'SETUP'
  INTEGER                           :: rcode, rcode2
  CHARACTER(LEN=80)                 :: wrfversion
  character                 :: ctmp3*3
  integer, intent(IN)       :: itimestep
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR OPENING WRF or FV3 NETCDF FILE', &
    & /, 1x, '***   FILE = ', a, &
    & /, 1x, '***   NCF:  ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9200 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNKNOWN WRF or FV3 OUTPUT VERSION', &
    & /, 1x, '***   IVERSION = ', i3, &
    & /, 1x, '***   GRIDTYPE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9300 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM WRF or FV3 FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   UNKNOWN OR UNSUPPORTED WRF or FV3 OUTPUT VERSION', &
    & /, 1x, '***   VERSION = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9500 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR CLOSING WRF or FV3 FILE', &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Try to determine if input meteorology file is in NetCDF format or not.
! If NetCDF format, it is probably WRF or FV3.
!-------------------------------------------------------------------------------
  met_model=inmetmodel

    !---------------------------------------------------------------------------
    ! If WRF, determine whether or not the Advanced Research WRF, ARW, formerly
    ! known as Eulerian mass, EM) version was used.
    !---------------------------------------------------------------------------

   IF ( met_model == 2 ) THEN

    rcode = nf90_open (file_mm(1), nf90_nowrite, cdfid)
    IF ( rcode .ne. nf90_noerr) then
      WRITE (*,f9000) TRIM(pname), TRIM(file_mm(1)), TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF
    
    rcode = nf90_get_att (cdfid, nf90_global, 'DYN_OPT', met_iversion)
    IF ( rcode /= nf90_noerr ) THEN
      rcode = nf90_get_att (cdfid, nf90_global, 'TITLE', wrfversion)
      IF ( rcode /= nf90_noerr ) THEN
        WRITE (*,f9300) TRIM(pname), 'TITLE', TRIM(nf90_strerror(rcode))
        CALL graceful_stop (pname)
      ENDIF
      IF ( wrfversion(18:19) >= "V3" ) THEN
        met_iversion = 2  ! NCAR only supports mass core in WRFv3 and beyond
      ELSE
        WRITE (*,f9400) TRIM(pname), TRIM(wrfversion)
        CALL graceful_stop (pname)
      ENDIF
    ENDIF

    rcode = nf90_get_att (cdfid, nf90_global, 'GRIDTYPE', gridtype)
    IF ( rcode /= nf90_noerr ) THEN
      WRITE (*,f9300) TRIM(pname), 'GRIDTYPE', rcode
      CALL graceful_stop (pname)
    ENDIF

    IF ( ( met_iversion == 2 ) .AND. ( gridtype(1:1) == "C" ) ) THEN
      CALL setup_wrfem (cdfid, ctmlays)
    ELSE
      WRITE (*,f9200) TRIM(pname), met_iversion, gridtype
      CALL graceful_stop (pname)
    ENDIF

    rcode = nf90_close (cdfid)
    IF ( rcode /= nf90_noerr ) THEN
      WRITE (*,f9500) TRIM(pname), TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

   ELSE  ! FV3
    write(ctmp3,'(i3.3)')itimestep
    rcode = nf90_open (trim(file_mm(1))//ctmp3//trim(file_mm(2)),nf90_nowrite, cdfid)
    rcode2 = nf90_open (trim(file_sfc(1))//ctmp3//trim(file_sfc(2)),nf90_nowrite, cdfid2)
 
    IF ( rcode.ne.nf90_noerr .or. rcode2.ne.nf90_noerr ) then
     WRITE (*,f9000) TRIM(pname),trim(file_mm(1))//ctmp3//trim(file_mm(2)), TRIM(nf90_strerror(rcode))
     WRITE (*,f9000) TRIM(pname),trim(file_sfc(1))//ctmp3//trim(file_sfc(2)), TRIM(nf90_strerror(rcode2))
     CALL graceful_stop (pname)
    endif 
     
    rcode = nf90_get_att (cdfid, nf90_global, 'source', fv3_version)

    IF ( rcode /= nf90_noerr ) THEN
      WRITE (*,f9300) TRIM(pname), TRIM(nf90_strerror(rcode))
      CALL graceful_stop (pname)
    ENDIF

    rcode = nf90_get_att (cdfid, nf90_global, 'grid', gridtype)
    IF ( rcode /= nf90_noerr ) THEN
      WRITE (*,f9300) TRIM(pname), 'grid', rcode
      CALL graceful_stop (pname)
    ENDIF

    IF ( ( fv3_version == "FV3GFS" ) .AND. ( gridtype(1:8) == "gaussian" ) ) THEN
      CALL setup_fv3 (cdfid, cdfid2, ctmlays)
    ELSE
      WRITE (*,f9200) TRIM(pname), fv3_version, gridtype
      CALL graceful_stop (pname)
    ENDIF
    rcode = nf90_close (cdfid)
    rcode = nf90_close (cdfid2)

   END IF

END SUBROUTINE setup
