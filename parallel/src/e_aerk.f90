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

FUNCTION e_aerk (tempc)

!-------------------------------------------------------------------------------
! Name:     Saturation Vapor Pressure
! Purpose:  Returns saturation vapor pressure [Pa] as a funtion of temperature
!           in degrees Celsius.
! Revised:  ?? ??? ????  Original version as a statement function in MCIP
!                        routines bcldprc_ak.f90 and getpblht.f90.
!           23 Feb 2011  Converted to independent routine.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!-------------------------------------------------------------------------------

  USE const

  IMPLICIT NONE

  REAL                         :: e_aerk   ! saturation vapor pressure [Pa]
  REAL,          INTENT(IN)    :: tempc    ! temperature [deg C]

  e_aerk = vp0 * EXP( 17.625 * tempc / ( 243.04 + tempc ) )

END FUNCTION e_aerk
