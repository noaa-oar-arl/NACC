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

SUBROUTINE myinterp (ain,met_nx,met_ny,aout,xindex,yindex,iout,jout,iflag)

!-------------------------------------------------------------------------------
! Name:     Horizontally Interpolate Grid To LCC Output for CMAQ 
! Purpose:  Interpolates input grid to LCC output for CMAQ using namelist 
!           specified grid definitions  
! Notes:       
! Revised:  26 Feb 2020 Adapted for FV3GFSv16 at NOAA-ARL (Y. Tang and 
!           P. C. Campbell)
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  REAL,               INTENT(IN)    :: ain        ( : , : )
  REAL,               INTENT(IN)    :: xindex     ( : , : )
  REAL,               INTENT(IN)    :: yindex     ( : , : )
  INTEGER,            INTENT(IN)    :: met_nx, met_ny
  INTEGER,            INTENT(IN)    :: iout, jout
  INTEGER,            INTENT(IN)    :: iflag
  REAL,               INTENT(OUT)   :: aout       ( : , : )

  INTEGER                           :: i, j
  
  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'MYINTERP'
  REAL,               ALLOCATABLE   :: x          ( : , : )
  REAL,               ALLOCATABLE   :: y          ( : , : )
  REAL,               ALLOCATABLE   :: xratio     ( : , : )
  REAL,               ALLOCATABLE   :: yratio     ( : , : )
!-------------------------------------------------------------------------------
! Error, warning, and informational messages.
!-------------------------------------------------------------------------------

  CHARACTER(LEN=256), PARAMETER :: f9000 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   INPUT ARRAY SIZES DO NOT MATCH', &
    & /, 1x, 70('*'))"

   CHARACTER(LEN=256), PARAMETER :: f9100 = "(/, 1x, 70('*'), &
    & /, 1x, '*** SUBROUTINE: ', a, &
    & /, 1x, '***   WRONG IFLAG (OPTIONS = 1 or 2)', &
    & /, 1x, 70('*'))"

   
!-------------------------------------------------------------------------------
! Check array sizes.
!-------------------------------------------------------------------------------

  IF ( ( SIZE(xindex,1) /= SIZE(yindex,1)  ) .AND.  &
       ( SIZE(xindex,2) /= SIZE(yindex,2) ) ) THEN
    WRITE (*,f9000) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF


!-------------------------------------------------------------------------------
! Check available interpolation options.
!-------------------------------------------------------------------------------
  IF ( iflag /= 1 .AND. iflag /= 2   ) THEN
    WRITE (*,f9100) TRIM(pname)
    CALL graceful_stop (pname)
  ENDIF

!-------------------------------------------------------------------------------
! Allocate necessary arrays.
!-------------------------------------------------------------------------------

  ALLOCATE ( x      ( iout, jout ) )
  ALLOCATE ( y      ( iout, jout ) )
  ALLOCATE ( xratio ( iout, jout ) )
  ALLOCATE ( yratio ( iout, jout ) )

!-------------------------------------------------------------------------------
! Compute horizontal interploation.
! 1 for nearest neighbor, 2 for bilinear
!-------------------------------------------------------------------------------
  do i=1,iout
   do j=1,jout
    if(iflag.eq.1) then ! nearest neighbor
     aout(i,j)=ain(nint(xindex(i,j)),nint(yindex(i,j)))
     else if(iflag.eq.2) then ! binear
     x(i,j)=xindex(i,j)
     y(i,j)=yindex(i,j)
     xratio(i,j)=x(i,j)-int(x(i,j))
     yratio(i,j)=y(i,j)-int(y(i,j))
     aout(i,j)=(1-yratio(i,j))*(ain(int(x(i,j)),int(y(i,j)))*         &
          (1-xratio(i,j))+ain(int(x(i,j))+1,int(y(i,j)))*xratio(i,j))+     &
          yratio(i,j)*(ain(int(x(i,j)),int(y(i,j))+1)*(1-xratio(i,j))+     &
          ain(int(x(i,j))+1,int(y(i,j))+1)*xratio(i,j))
    endif
   enddo
  enddo
 
!-------------------------------------------------------------------------------
! Deallocate arrays.
!-------------------------------------------------------------------------------

  DEALLOCATE ( x       )
  DEALLOCATE ( y       )
  DEALLOCATE ( xratio  )
  DEALLOCATE ( yratio  )

END SUBROUTINE myinterp
