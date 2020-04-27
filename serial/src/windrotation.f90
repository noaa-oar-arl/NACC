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

SUBROUTINE windrotation (ua,va,londot,ix,jy,reflon,reflat)

!-------------------------------------------------------------------------------
! Name:     Wind Rotation
! Purpose:  Rotates U and V components with respect to true north 
!           (meteorological standard)
! Notes:       
! Revised:  25 Feb 2020 Adapted for FV3GFSv16 at NOAA-ARL (Y. Tang and 
!           P. C. Campbell)
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  REAL,               INTENT(INOUT) :: ua         ( : , : )
  REAL,               INTENT(INOUT) :: va         ( : , : )
  REAL,               INTENT(IN)    :: londot     ( : , : )
  INTEGER,            INTENT(IN)    :: ix, jy
  REAL,               INTENT(IN)    :: reflon, reflat
  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'WINDROTATION'
  REAL,               ALLOCATABLE   :: angle2     ( : , : )
  REAL,               ALLOCATABLE   :: sinx2      ( : , : )
  REAL,               ALLOCATABLE   :: cosx2      ( : , : )
  REAL,               ALLOCATABLE   :: ut         ( : , : )
  REAL,               ALLOCATABLE   :: vt         ( : , : )
  INTEGER                           :: i, j
  REAL                              :: d2r, rotcon_p
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   INPUT ARRAY SIZES DO NOT MATCH', &
    & /, 1x, 70('*'))"

!-------------------------------------------------------------------------------
! Check array sizes.
!-------------------------------------------------------------------------------

  IF ( ( SIZE(ua,1) /= SIZE(londot,1)  ) .AND.  &
       ( SIZE(ua,2) /= SIZE(londot,2) ) ) THEN
    WRITE (*,f9000) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

  IF ( ( SIZE(va,1) /= SIZE(londot,1)  ) .AND.  &
       ( SIZE(va,2) /= SIZE(londot,2) ) ) THEN
    WRITE (*,f9000) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

!-------------------------------------------------------------------------------
! Allocate necessary arrays.
!-------------------------------------------------------------------------------

  ALLOCATE ( angle2 ( ix, jy ) )
  ALLOCATE ( sinx2  ( ix, jy ) )
  ALLOCATE ( cosx2  ( ix, jy ) )
  ALLOCATE ( ut     ( ix, jy ) )
  ALLOCATE ( vt     ( ix, jy ) )

!-------------------------------------------------------------------------------
! Compute wind rotation.
! from https://www.mcs.anl.gov/~emconsta/wind_conversion.txt
! or   https://github.com/matplotlib/basemap/issues/269
!-------------------------------------------------------------------------------
      d2r=atan(1.)/45.
      rotcon_p=sin(reflat*d2r)
        do j=1,jy
        do i=1,ix
           angle2(i,j) = rotcon_p*(londot(i,j)-reflon)*d2r
           sinx2(i,j) = sin(angle2(i,j))
           cosx2(i,j) = cos(angle2(i,j))
           ut(i,j) = ua(i,j)
           vt(i,j) = va(i,j)
           ua(i,j) = cosx2(i,j)*ut(i,j)-sinx2(i,j)*vt(i,j)
           va(i,j) = sinx2(i,j)*ut(i,j)+cosx2(i,j)*vt(i,j)
         end do
        end do
          
!-------------------------------------------------------------------------------
! Deallocate arrays.
!-------------------------------------------------------------------------------

  DEALLOCATE ( angle2 )
  DEALLOCATE ( sinx2  )
  DEALLOCATE ( cosx2  )
  DEALLOCATE ( ut     )
  DEALLOCATE ( vt     )

END SUBROUTINE windrotation
