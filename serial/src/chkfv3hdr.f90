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

SUBROUTINE chkfv3hdr (fl, cdfid)

!-------------------------------------------------------------------------------
! Name:     Check WRF Header
! Purpose:  Check WRF header variables from one WRF output file against the
!           "base" WRF output file used for this MCIP run to ensure that the
!           WRF output files are from the same simulation.
! Notes:    This routine is not thorough, but it should be enough to spot-check
!           key variables that would indicate a different WRF simulation.
!           This routine assumes that FL (input argument) is already opened.
! Revised:  15 May 2008  Original version.  (T. Otte)
!           11 May 2009  Correct bug in checking surface layer scheme in
!                        subsequent files.  (T. Otte)
!           25 Sep 2009  Removed netCDF file opening to prevent condition with
!                        too many open files for long simulations.  Added
!                        check on urban physics option and surface analysis
!                        nudging options.  Changed code to allow for GRID_FDDA
!                        to be greater than 1.  Corrected bug in checking to
!                        ensure that data are from same simulation (i.e.,
!                        restarted rather than reinitialized).  Corrected bug
!                        in checking observation nudging coefficient for
!                        temperature.  (T. Otte)
!           12 Feb 2010  Removed unused variables CDFID, DX, DY, N_TIMES, and
!                        VARID, and removed unused format statements 9600 and
!                        9700.  Changed RTOL to 1.0e-4.  (T. Otte)
!           18 Mar 2010  Added CDFID as an input argument.  Removed dependency
!                        on module WRF_NETCDF.  (T. Otte)
!           23 Aug 2011  Changed name of module FILE to FILES to avoid conflict
!                        with F90 protected intrinsic.  Updated netCDF commands
!                        to F90, and improved error handling.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           24 Feb 2020  Adapted for FV3GFSv16 at NOAA-ARL (P. C. Campbell)
!-------------------------------------------------------------------------------

  USE files
  USE metinfo
  USE netcdf
  USE mcipparm

  IMPLICIT NONE

  INTEGER,            INTENT(IN)    :: cdfid
  CHARACTER(LEN=80)                 :: cval
  CHARACTER(LEN=8)                  :: cval8
  INTEGER                           :: dimid
  CHARACTER(LEN=256), INTENT(IN)    :: fl
  CHARACTER(LEN=256)                :: fl1
  INTEGER                           :: ival
  CHARACTER(LEN=16),  PARAMETER     :: pname     = 'CHKFV3HDR'
  INTEGER                           :: rcode
  REAL,               PARAMETER     :: rtol      = 1.0e-4
  REAL                              :: rval

!-------------------------------------------------------------------------------
! Error-handling section.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   FV3 FILES DO NOT SEEM TO BE FROM SAME DOMAIN', &
    & /, 1x, '***   VARIABLE = ', a)"

  CHARACTER(LEN=256), PARAMETER :: f9100 = "( &
    & /, 1x, '***   FIRST FILE = ', a, &
    & /, 1x, '***   VALUE IN FIRST FILE = ', i4, &
    & /, 1x, '***   NEW FILE = ', a, &
    & /, 1x, '***   VALUE IN NEW FILE = ', i4, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9200 = "( &
    & /, 1x, '***   FIRST FILE = ', a, &
    & /, 1x, '***   VALUE IN FIRST FILE = ', f13.3, &
    & /, 1x, '***   NEW FILE = ', a, &
    & /, 1x, '***   VALUE IN NEW FILE = ', f13.3, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9300 = "( &
    & /, 1x, '***   FIRST FILE = ', a, &
    & /, 1x, '***   VALUE IN FIRST FILE = ', a, &
    & /, 1x, '***   NEW FILE = ', a, &
    & /, 1x, '***   VALUE IN NEW FILE = ', a, &
    & /, 1x, 70('*'))"

  CHARACTER(LEN=256), PARAMETER :: f9400 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   ERROR RETRIEVING VARIABLE FROM WRF FILE', &
    & /, 1x, '***   VARIABLE = ', a, &
    & /, 1x, '***   NCF: ', a, &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Check NX, NY, and NZ.
!-------------------------------------------------------------------------------

  fl1 = file_mm(1)

  rcode = nf90_get_att (cdfid, nf90_global, 'im', ival)
  IF ( rcode == nf90_noerr ) THEN
    IF ( ival /= met_nx ) THEN
      WRITE (*,f9000) TRIM(pname), 'WEST-EAST_GRID_DIMENSION'
      WRITE (*,f9100) TRIM(fl1), met_nx, TRIM(fl), ival
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    WRITE (*,f9400) TRIM(pname), 'WEST-EAST_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_get_att (cdfid, nf90_global, 'jm', ival)
  IF ( rcode == nf90_noerr ) THEN
    IF ( ival /= met_ny ) THEN
      WRITE (*,f9000) TRIM(pname), 'SOUTH-NORTH_GRID_DIMENSION'
      WRITE (*,f9100) TRIM(fl1), met_ny, TRIM(fl), ival
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    WRITE (*,f9400) TRIM(pname), 'SOUTH-NORTH_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inq_dimid (cdfid, 'phalf', dimid)
  IF ( rcode /= nf90_noerr ) THEN
    WRITE (*,f9400) TRIM(pname), 'ID for phalf',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

  rcode = nf90_inquire_dimension (cdfid, dimid, len=ival)
  IF ( rcode == nf90_noerr ) THEN
    IF ( ival-1 /= met_nz ) THEN
      WRITE (*,f9000) TRIM(pname), 'BOTTOM-TOP_GRID_DIMENSION'
      WRITE (*,f9100) TRIM(fl1), met_nz, TRIM(fl), MIN(maxlays,ival)
      CALL graceful_stop (pname)
    ENDIF
  ELSE
    WRITE (*,f9400) TRIM(pname), 'BOTTOM-TOP_GRID_DIMENSION',  &
                    TRIM(nf90_strerror(rcode))
    CALL graceful_stop (pname)
  ENDIF

END SUBROUTINE chkfv3hdr
