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

SUBROUTINE xy2ll_ps (xx, yy, phi1, phi2, lambda0, phi, lambda)

!-------------------------------------------------------------------------------
! Name:     (X,Y) to Latitude-Longitude for Polar Secant Stereographic Projection
! Purpose:  Calcluates latitude-longitude for a given (X,Y) pair from origin
!           and Polar Secant Stereographic projection information.
! Notes:    This subroutine converts from Polar Stereographic (X,Y) coordinates  
!           to geodetic latitude and longitude for the polar regions. The
!           equations are from Snyder, J. P., 1982,  Map Projections Used by the 
!           U.S. Geological Survey, Geological Survey Bulletin 1532, U.S. 
!           Government Printing Office.  See JPL Technical Memorandum 
!           3349-85-101 for further details.  
!           http://www.ccpo.odu.edu/~msd/GOVARS/mapxy.m                                                  
!           https://www.mathworks.com/matlabcentral/fileexchange/32907-polar-stereographic-coordinate-transformation-map-to-lat-lon
! Revised:  12 Dec 2022  Original version.  (P.C. Campbell)
!-------------------------------------------------------------------------------

  USE const, ONLY: rearth

  IMPLICIT NONE

  REAL(8)                      :: deg2rad ! convert degrees to radians
  REAL(8)                      :: rad2deg ! convert radians to degrees
  REAL(8)                      :: drearth ! earth radius [m]
  REAL(8)                      :: hemi    ! +/-1 for input Northern/Southern Hemis 
  REAL,          INTENT(OUT)   :: lambda  ! longitude [deg]
  REAL,          INTENT(IN)    :: lambda0 ! central meridian [deg]
  REAL(8)                      :: lambda02rad ! central meridian [rad]
  REAL,          INTENT(OUT)   :: phi     ! latitude [deg]
  REAL,          INTENT(IN)    :: phi1    ! 1.0 for North Polar, -1.0 for South Polar
  REAL,          INTENT(IN)    :: phi2    ! Secant latitude (latitude of true scale)
  REAL(8)                      :: phi2rad ! true latitude 1 [rad]
  REAL(8)                      :: pi
  REAL(8)                      :: piover4 ! pi/4
  REAL(8)                      :: rho
  REAL(8)                      :: tt, cm
  REAL(8)                      :: ee, e
  REAL(8)                      :: chi
  REAL,          INTENT(IN)    :: xx      ! X-coordinate from origin
  REAL,          INTENT(IN)    :: yy      ! Y-coordinate from origin
  REAL(8)                      :: xxd, yyd
!-------------------------------------------------------------------------------
! Compute constants.
!-------------------------------------------------------------------------------

    piover4 = DATAN(1.0d0)
    pi      = 4.0d0 * piover4
    deg2rad = pi / 1.8d2
    rad2deg = 1.8d2 / pi
    ee = 6.694379852e-3           !Eccentricity squared
    e  = sqrt(ee)
    drearth = DBLE(rearth)
    xxd = DBLE(xx)     
    yyd = DBLE(yy)     
    hemi = phi1

    phi2rad = DBLE(phi2) * deg2rad  ! convert true latitude from degrees to radians
    lambda02rad = DBLE(lambda0) * deg2rad  ! convert central meridian from degrees to radians

    !if the standard parallel is in S.Hemi., switch signs.
    if (phi2rad.lt.0) then
      hemi=-1        !plus or minus, north lat. or south
      phi2rad=-1*phi2rad
      lambda02rad=-1*lambda02rad
      xxd=-1*xxd
      yyd=-1*yyd
    else
      hemi=1
    end if

!this is not commented very well. See Snyder for details.
    tt=tan(pi/4 - phi2rad/2)/((1-e*sin(phi2rad))/(1+e*sin(phi2rad)))**(e/2)
    cm=cos(phi2rad)/sqrt(1-e**2*(sin(phi2rad))**2)
    rho=sqrt(xxd**2+yyd**2)
    tt=rho*tt/(drearth*cm)

!find phi (lat) with a series instead of iterating.
    chi=pi/2 - 2 * atan(tt);
    phi=chi+(e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360)*sin(2*chi) &
    + (7*e**4/48 + 29*e**6/240 + 811*e**8/11520)*sin(4*chi) &
    + (7*e**6/120+81*e**8/1120)*sin(6*chi) &
    + (4279*e**8/161280)*sin(8*chi)

!find lambda (lon)
    lambda=lambda02rad + atan2(xxd,-yyd)

!correct the signs and phasing
    phi=hemi*phi
    lambda=hemi*lambda
    lambda=mod(lambda+pi,2*pi)-pi !want longitude in the range -pi to pi

!convert back to degrees
    phi=rad2deg*phi
    lambda=rad2deg*lambda

END SUBROUTINE xy2ll_ps
