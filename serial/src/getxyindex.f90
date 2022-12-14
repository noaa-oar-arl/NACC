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

SUBROUTINE getxyindex (xlat,xlon,xi,yj,tlat,tlon,ix,jy)

!-------------------------------------------------------------------------------
! Name:     Get XY Index in a target domain 
! Purpose:  Used for interpolates input grid to LCC output for CMAQ using namelist 
!           specified grid definitions  
! Notes:       
! Revised:  26 Feb 2020 Adapted for FV3GFSv16 at NOAA-ARL (Y. Tang and 
!           P. C. Campbell)
!           14 Dec 2022 Modified to account for polar stereographic projections
!           (Y. Tang and P. C. Campbell)
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  REAL,               INTENT(IN)    :: xlat        
  REAL,               INTENT(IN)    :: xlon        
  REAL,               INTENT(IN)    :: tlat ( : )       
  REAL,               INTENT(IN)    :: tlon ( : )       
  INTEGER,            INTENT(IN)    :: ix, jy
  REAL,               INTENT(OUT)   :: xi
  REAL,               INTENT(OUT)   :: yj

  CHARACTER(LEN=16),  PARAMETER     :: pname      = 'GETXYINDEX'
  INTEGER                           :: i, j
  REAL                              :: xlontmp    

!-------------------------------------------------------------------------------
! Get XY Index in a target domain.
! 
!-------------------------------------------------------------------------------
 xlontmp=xlon
 if(xlontmp.lt.0) xlontmp=xlontmp+360

 if(xlontmp.gt.tlon(ix).and.xlontmp.le.tlon(1)+360) then
  xi=ix+(xlontmp-tlon(ix))/(tlon(1)-tlon(ix)+360)
 else if (xlontmp.le.tlon(ix)) then 
  do i=1,ix-1
   if(xlontmp.ge.tlon(i).and.xlontmp.le.tlon(i+1)) then
    xi=i+(xlontmp-tlon(i))/(tlon(i+1)-tlon(i))
    exit
   endif
  enddo
 else
  print*,'wrong xlontmp ',xlontmp,tlon(ix),tlon(1)
 endif

  yj=0.
  do j=1,jy-1
   if(xlat.le.amax1(tlat(j),tlat(j+1)).and.xlat.ge.amin1(tlat(j),tlat(j+1))) then
! fv3 is from north to south
    yj=j+(xlat-tlat(j))/(tlat(j+1)-tlat(j))
    exit
   endif
  enddo
  if(yj.le.1e-32) then
     if(tlat(1).gt.tlat(2)) then ! north to south
        if(xlat.ge.tlat(1)) then
          yj=1.
        else if (xlat.le.tlat(jy)) then
          yj=jy
        endif
      else  ! south to north
        if(xlat.le.tlat(1)) then
           yj=1.
        else if (xlat.ge.tlat(jy)) then
           yj=jy
        endif
      endif
   if(yj.le.1e-32) then
      print*,'wrong xlat ',xlat
      stop
   endif
 endif
 
END SUBROUTINE getxyindex
