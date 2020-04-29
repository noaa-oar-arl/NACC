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

SUBROUTINE alloc_met

!-------------------------------------------------------------------------------
! Name:     Allocate Meteorology Variables
! Purpose:  Allocate arrays for input meteorology variables.
! Revised:  10 Sep 2001  Original version.  (T. Otte)
!           29 May 2003  Added SNOWCOVR.  (D. Schwede)
!           09 Aug 2004  Added QGA, VEGOLD, and T2.  (D. Schwede and T. Otte)
!           29 Nov 2004  Added LUFRAC.  (T. Otte)
!           04 Apr 2005  Changed array dimensions from I,J to X,Y to make code
!                        more general.  Now all meteorology arrays will have
!                        the east-west index first, as opposed to the standard
!                        MM5 convention.  For MM5 input, arrays are transposed
!                        to the new convention in either RDMM5V2 or RDMM5V3.
!                        Removed unused variables REGIME and MAVAIL.  Added PH,
!                        PHB, PB, MU, and MUB for WRF.  Added U10 and V10.
!                        (T. Otte and S.-B. Kim)
!           11 Aug 2005  Removed unused variable FSOIL.  (T. Otte)
!           25 Jul 2007  Removed internal variables for emissivity and net
!                        radiation.  Eliminated logical variable "PX" to make
!                        code more general.  (T. Otte)
!           05 May 2008  Added 2-m mixing ratio (Q2) and turbulent kinetic
!                        energy (TKE) arrays.  Added urban fraction (FRC_URB)
!                        and urban roughness length (Z0C_URB2D) for
!                        MET_UCMCALL=1.  (T. Otte)
!           29 Oct 2009  Changed MET_UCMCALL to MET_URBAN_PHYS, and allowed
!                        for variable to be set to be greater than 1.  Added
!                        THETA and CORIOLIS for when potential vorticity is
!                        needed.  Added LATU, LONU, MAPU, LATV, LONV, and
!                        and MAPV.  Removed Z0C_URB2D.  (T. Otte)
!           15 Dec 2010  Added sea ice.  Added tipping buckets for convective
!                        and non-convective precipitation.  (T. Otte)
!           07 Sep 2011  Updated disclaimer.  (T. Otte)
!           11 Sep 2012  Added LANDMASK to be read from WRF.  (T. Otte)
!           10 Apr 2015  Added new array CLDFRA to pass 3D resolved cloud
!                        fraction to output.  (T. Spero)
!           21 Aug 2015  Changed latent heat flux from QFX to LH.  Fill THETA
!                        and add moisture flux (QFX) for IFMOLACM.  (T. Spero)
!           17 Sep 2015  Changed IFMOLACM to IFMOLPX.  (T. Spero)
!           16 Mar 2018  Added SNOWH to output.  Added C1H, C2H, C1F, and C2F to
!                        support hybrid vertical coordinate in WRF.  Added
!                        LUFRAC2, MOSCATIDX, ZNT_MOS, TSK_MOS, RA_MOS, RS_MOS,
!                        and LAI_MOS for NOAH Mosaic land-surface model.
!                        Added DZS, SOIT3D, and SOIM3D.  Added WSPDSFC and
!                        XLAIDYN for Noah.  (T. Spero)
!           27 Jun 2018  Removed local aliases for dimensions of input
!                        meteorological fields.  (T. Spero)
!           14 Sep 2018  Changed condition to enable hybrid vertical coordinate
!                        from WRF.  Removed support for MM5v3 input.  (T. Spero)
!           18 Jun 2019  Added new surface variables with PX LSM that can
!                        improve dust simulation in CCTM.  Added optional
!                        variables from KF convective scheme with radiative
!                        feedbacks.  (T. Spero)
!           24 Feb 2020  Adapted for FV3GFSv16 at NOAA-ARL (P. C. Campbell)
!-------------------------------------------------------------------------------

  USE metinfo
  USE metvars
  USE mcipparm

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Allocate time-invariant fields.
!-------------------------------------------------------------------------------
  integer ix, jy
  if(met_model.eq.2) then
   ix=met_nx
   jy=met_ny
  else if (met_model.eq.3) then ! fv3
   ix=ncols_x+1
   jy=nrows_x+1
!   allocate(dum3d(met_nx,met_ny,met_nz))
  else
   print*,'wrong met_model ',met_model
   stop
  endif
       
  ALLOCATE ( albedo   (ix, jy) )   ! time varying in P-X LSM
  ALLOCATE ( landmask (ix, jy) )   ! time varying in NOAH LSM
  ALLOCATE ( landuse  (ix, jy) )
  ALLOCATE ( latcrs   (ix, jy) )
  ALLOCATE ( latdot   (ix, jy) )
  ALLOCATE ( latu     (ix, jy) )
  ALLOCATE ( latv     (ix, jy) )
  ALLOCATE ( loncrs   (ix, jy) )
  ALLOCATE ( londot   (ix, jy) )
  ALLOCATE ( lonu     (ix, jy) )
  ALLOCATE ( lonv     (ix, jy) )
  ALLOCATE ( mapcrs   (ix, jy) )
  ALLOCATE ( mapdot   (ix, jy) )
  ALLOCATE ( mapu     (ix, jy) )
  ALLOCATE ( mapv     (ix, jy) )
  ALLOCATE ( sigmaf                   (met_nz+1) )
  ALLOCATE ( sigmah                   (met_nz) )
  IF ( met_model == 3 ) THEN  ! FV3
   ALLOCATE ( pfull(met_nz+1) )
   ALLOCATE ( phalf(met_nz) )
   ALLOCATE ( dpres(ix, jy,met_nz) )
   ALLOCATE ( delz (ix, jy, met_nz) )
   ALLOCATE ( b_k  (met_nz+1) )
  ENDIF
  ALLOCATE ( terrain  (ix, jy) )
  ALLOCATE ( znt      (ix, jy) )

  IF ( iflufrc ) THEN
    ALLOCATE ( lufrac (ix, jy, nummetlu) )
    IF ( ifmosaic ) THEN
      ALLOCATE ( lufrac2   (ix, jy, nummetlu) )
      ALLOCATE ( moscatidx (ix, jy, nummetlu) )
    ENDIF
  ENDIF

  IF ( lpv > 0 ) THEN  ! potential vorticity; get Coriolis
    ALLOCATE ( coriolis (ix, jy) )
  ENDIF

   IF ( met_hybrid >= 0 ) THEN
    ALLOCATE ( c1f (met_nz+1) )
    ALLOCATE ( c1h (met_nz)   )
    ALLOCATE ( c2f (met_nz+1) )
    ALLOCATE ( c2h (met_nz)   )
   ENDIF

  IF ( met_ns > 0 ) THEN
    ALLOCATE ( dzs (met_ns) )
  ENDIF

!-------------------------------------------------------------------------------
! Allocate time-varying fields.
!-------------------------------------------------------------------------------

  ALLOCATE ( glw     (ix, jy) )
  ALLOCATE ( groundt (ix, jy) )
  ALLOCATE ( hfx     (ix, jy) )
  ALLOCATE ( i_rainc (ix, jy) )
  ALLOCATE ( i_rainnc(ix, jy) )
  ALLOCATE ( ircold  (ix, jy) )
  ALLOCATE ( irnold  (ix, jy) )
  ALLOCATE ( lh      (ix, jy) )
  ALLOCATE ( pp      (ix, jy, met_nz) )
  ALLOCATE ( psa     (ix, jy) )
  ALLOCATE ( qca     (ix, jy, met_nz) )
  ALLOCATE ( qga     (ix, jy, met_nz) )
  ALLOCATE ( qia     (ix, jy, met_nz) )
  ALLOCATE ( qra     (ix, jy, met_nz) )
  ALLOCATE ( qsa     (ix, jy, met_nz) )
  ALLOCATE ( qva     (ix, jy, met_nz) )
  ALLOCATE ( raincon (ix, jy) )
  ALLOCATE ( rainnon (ix, jy) )
  ALLOCATE ( rcold   (ix, jy) )   ! save this variable on each call
  ALLOCATE ( rgrnd   (ix, jy) )
  ALLOCATE ( rnold   (ix, jy) )   ! save this variable on each call
  ALLOCATE ( seaice  (ix, jy) )
  ALLOCATE ( snowcovr(ix, jy) )
  ALLOCATE ( snowh   (ix, jy) )
  ALLOCATE ( ta      (ix, jy, met_nz) )
  ALLOCATE ( ua      (ix, jy, met_nz) )
  ALLOCATE ( ust     (ix, jy) )
  ALLOCATE ( va      (ix, jy, met_nz) )
  ALLOCATE ( wa      (ix, jy, met_nz) )
  ALLOCATE ( zpbl    (ix, jy) )

  IF ( ift2m ) THEN  ! 2-m temperature available
    ALLOCATE ( t2    (ix, jy) )
  ENDIF

  IF ( ifq2m ) THEN  ! 2-m mixing ratio available
    ALLOCATE ( q2    (ix, jy) )
  ENDIF

  IF ( ifw10m ) THEN  ! 10-m wind components available
    ALLOCATE ( u10   (ix, jy) )
    ALLOCATE ( v10   (ix, jy) )
  ENDIF

  IF ( met_model == 2 ) THEN  ! WRF
    ALLOCATE ( mu    (ix, jy) )
    ALLOCATE ( mub   (ix, jy) )
    ALLOCATE ( pb    (ix, jy, met_nz)   )
    ALLOCATE ( ph    (ix, jy, met_nz+1) )
    ALLOCATE ( phb   (ix, jy, met_nz+1) )
  ELSE 
    !test
  ENDIF

  IF ( iflai ) THEN  ! leaf area index available
    ALLOCATE ( lai    (ix, jy) )
  ENDIF

  IF ( ifmol ) THEN  ! Monin-Obukhov length available
    ALLOCATE ( mol    (ix, jy) )
  ENDIF

  IF ( ifresist ) THEN  ! aerodynamic and stomatal resistances available
    ALLOCATE ( ra     (ix, jy) )
    ALLOCATE ( rstom  (ix, jy) )
  ENDIF

  IF ( ifveg ) THEN  ! vegetation fraction available
    ALLOCATE ( veg    (ix, jy) )
  ENDIF

  IF ( ifwr ) THEN  ! canopy wetness available
    ALLOCATE ( wr     (ix, jy) )
  ENDIF

  IF ( ifsoil ) THEN  ! soil moisture, temperature, and type available
    ALLOCATE ( isltyp (ix, jy) )
    ALLOCATE ( soilt1 (ix, jy) )
    ALLOCATE ( soilt2 (ix, jy) )
    ALLOCATE ( w2     (ix, jy) )
    ALLOCATE ( wg     (ix, jy) )
    ALLOCATE ( soim3d (ix, jy, met_ns) )
    ALLOCATE ( soit3d (ix, jy, met_ns) )
  ENDIF
  
  IF ( iftke ) THEN  ! turbulent kinetic energy available
    IF ( iftkef ) THEN  ! TKE on full-levels
      ALLOCATE ( tke   (ix, jy, met_nz+1) )
    ELSE  ! TKE on half-levels
      ALLOCATE ( tke   (ix, jy, met_nz) )
    ENDIF
  ENDIF

  IF ( lpv > 0 .OR. ifmolpx ) THEN  ! need potential temperature
    ALLOCATE ( theta (ix, jy, met_nz) )
  ENDIF

  IF ( ifmolpx ) THEN  ! recalculate Monin-Obukhov length for WRF-ACM2
    ALLOCATE ( qfx   (ix, jy) )
  ENDIF

  IF ( met_urban_phys >= 1 ) THEN  ! urban canopy model in WRF
    ALLOCATE ( frc_urb   (ix, jy) )
  ENDIF

  IF ( ifcld3d ) THEN
    ALLOCATE ( cldfra (ix, jy, met_nz) )
  ENDIF

  IF ( ifmosaic ) THEN
    ALLOCATE ( lai_mos (ix, jy, nummosaic) )
    ALLOCATE ( ra_mos  (ix, jy, nummosaic) )
    ALLOCATE ( rs_mos  (ix, jy, nummosaic) )
    ALLOCATE ( tsk_mos (ix, jy, nummosaic) )
    ALLOCATE ( znt_mos (ix, jy, nummosaic) )
    ALLOCATE ( wspdsfc (ix, jy) )
    ALLOCATE ( xlaidyn (ix, jy) )
  ENDIF

  IF ( ifpxwrf41 ) THEN
    ALLOCATE ( lai_px    (ix, jy) )
    ALLOCATE ( wsat_px   (ix, jy) )
    ALLOCATE ( wfc_px    (ix, jy) )
    ALLOCATE ( wwlt_px   (ix, jy) )
    ALLOCATE ( csand_px  (ix, jy) )
    ALLOCATE ( fmsand_px (ix, jy) )
    ALLOCATE ( clay_px   (ix, jy) )
  ENDIF

  IF ( ifkfradextras ) THEN
    ALLOCATE ( qc_cu     (ix, jy, met_nz) )
    ALLOCATE ( qi_cu     (ix, jy, met_nz) )
    ALLOCATE ( cldfra_dp (ix, jy, met_nz) )
    ALLOCATE ( cldfra_sh (ix, jy, met_nz) )
  ENDIF

END SUBROUTINE alloc_met
