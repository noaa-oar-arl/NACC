
PROGRAM BDYSCALE

    !!***********************************************************************
    !! Version "$Id: BDYSCALE.f90 1 2017-06-10 18:05:20Z coats $"
    !!   EDSS/Models-3 M3TOOLS.
    !!   Copyright (C) 2017 UNC Institute for the Environment.
    !!   Distributed under the GNU GENERAL PUBLIC LICENSE version 2
    !!   See file "GPL.txt" for conditions of use.
    !!.........................................................................
    !!  program body starts at line  140
    !!
    !!  DESCRIPTION:
    !!      For each time step, reads all variables, rescale the 
    !!      north and west edges by a specified scalar, , and write them
    !!       to the specified output file.
    !!
    !!  REVISION  HISTORY:
    !!       Prototype 8/2017 by Carlie J. Coats, Jr., UNC IE
    !!***********************************************************************
    USE M3UTILIO


    IMPLICIT NONE

    !!...........   PARAMETERS and their descriptions:

    CHARACTER*16, PARAMETER :: PNAME = 'BDYSCALE'
    CHARACTER*16, PARAMETER :: BLANK = ' '
    CHARACTER*72, PARAMETER :: BAR   =  &
    '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'


    !!...........   LOCAL VARIABLES and their descriptions:

    CHARACTER*16    FNAME   !  input data  file logical name
    LOGICAL         EFLAG
    CHARACTER*256   MESG
    INTEGER         LDEV        !  log-device
    INTEGER         ISTAT       !  allocation-status

    INTEGER         C, R, L, V, N

    INTEGER         SDATE, STIME, NRECS, JDATE, JTIME, TSTEP

    INTEGER         WEST0, WEST1, NORTH0, NORTH1, PERIM

    INTEGER         NCOLS1      ! number of grid columns
    INTEGER         NROWS1      ! number of grid rows
    INTEGER         NLAYS1      ! number of layers
    INTEGER         GDTYP1      ! grid type:  1=LAT-LON, 2=UTM, ...
    REAL*8          P_ALP1      ! first, second, third map
    REAL*8          P_BET1      ! projection descriptive
    REAL*8          P_GAM1      ! parameters.
    REAL*8          XCENT1      ! lon for coord-system X=0
    REAL*8          YCENT1      ! lat for coord-system Y=0
    REAL*8          XORIG1      ! X-coordinate origin of grid (map units)
    REAL*8          YORIG1      ! Y-coordinate origin of grid
    REAL*8          XCELL1      ! X-coordinate cell dimension
    REAL*8          YCELL1      ! Y-coordinate cell dimension

    REAL            NWFAC

    REAL, ALLOCATABLE :: BUF1( :,: )


    !!***********************************************************************
    !!   begin body of program BDYSCALE

    LDEV  = INIT3()
    EFLAG = .FALSE.

  WRITE( *, '( 5X, A )' ) BLANK, BAR, BLANK,                                &
'Program BDYSCALE to read all time steps of all variables from a boundary', &
'conditions file, re-scale the north and west edges by a specified',        &
'scalar, and write the result to the specified output file.',               &
'',                                                                         &
'PRECONDITIONS REQUIRED:',                                                  &
'',                                                                         &
'    setenv BDY_IN      <path-name>',                                       &
'    setenv BDY_OUT     <path-name>',                                       &
'    setenv BDY_FAC     <scalar>',                                          &
'',                                                                         &
'    ${BDY_IN} is an I/O API boundary-condition file.',                     &
'',                                                                         &
'Program copyright (C) 1992-2002 MCNC, (C) 1995-2013 Carlie J. Coats, Jr.', &
'(C) 2003-2010 Baron Advanced Meteorological Systems, LLC., and',           &
'(C) 2014-2016 UNC Institute for the Environment.',                         &
'Released under Version 2 of the GNU General Public License. See',          &
'enclosed GPL.txt, or URL',                                                 &
''  ,                                                                       &
'    https://www.gnu.org/licenses/old-licenses/gpl-2.0.html',               &
''  ,                                                                       &
'Comments and questions are welcome and can be sent to'  ,                  &
'',                                                                         &
'    Carlie J. Coats, Jr.    carlie@jyarborough.com',                       &
'or',                                                                       &
'    UNC Institute for the Environment',                                    &
'    100 Europa Dr., Suite 490 Rm 405',                                     &
'    Campus Box 1105',                                                      &
'    Chapel Hill, NC 27599-1105',                                           &
'',                                                                         &
'Program version: ',                                                        &
'$Id: BDYSCALE.f90 1 2017-06-10 18:05:20Z coats $',&
' '

    IF ( .NOT. GETVAL( 'Continue with program?', .TRUE. ) ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Program terminated at user request', 2 )
    END IF


    !!...............  Open and get description for input data file

    IF (      .NOT.OPEN3( 'BDY_IN', FSREAD3, PNAME ) ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Could not open "BDY_IN"', 2 )
    ELSE IF ( .NOT.DESC3( 'BDY_IN' ) ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Could not DESC3( "BDY_IN" )', 2 )
    ELSE
        NCOLS1 = NCOLS3D
        NROWS1 = NROWS3D
        NLAYS1 = NLAYS3D
        GDTYP1 = GDTYP3D
        P_ALP1 = P_ALP3D
        P_BET1 = P_BET3D
        P_GAM1 = P_GAM3D
        XCENT1 = XCENT3D
        YCENT1 = YCENT3D
        XORIG1 = XORIG3D
        YORIG1 = YORIG3D
        XCELL1 = XCELL3D
        YCELL1 = YCELL3D
        TSTEP  = TSTEP3D
        SDATE  = SDATE3D
        STIME  = STIME3D
        NRECS  = MXREC3D
        NORTH0 =   NCOLS1   + NROWS1 + 3
        NORTH1 = 2*NCOLS1 +   NROWS1 + 3      !!  note that north1+1 = west0
        WEST0  = 2*NCOLS1 +   NROWS1 + 4
        WEST1  = 2*NCOLS1 + 2*NROWS1 + 4
        PERIM  = 2*NCOLS1 + 2*NROWS1 + 4
    END IF

    !!...............  Get scale factor

    NWFAC = ENVGET( 'BDY_FAC', 'Scale factor for north, west boundary-edges', 1.0, ISTAT )
    IF ( ISTAT .GT. 0 ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Bad environment variable "BDY_FAC"', 2 )
    END IF

    !!...............  Open output data file (borrowing file description...)

    IF (      .NOT.OPEN3( 'BDY_OUT', FSUNKN3, PNAME ) ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Could not open "BDY_OUT"', 2 )
    END IF


    !!...............  Allocate buffers; compute re-gridding matrix

    ALLOCATE( BUF1( PERIM, NLAYS1 ),  STAT = ISTAT )
    IF ( ISTAT .NE. 0 ) THEN
        WRITE( MESG, '( A, I10)' ) 'Buffer allocation failed:  STAT=', ISTAT
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF


    !!...............  Process output time step sequence

    JDATE = SDATE
    JTIME = STIME

    DO  N = 1, NRECS

        WRITE( MESG, '( A, I7.7, A, I6.6 )' ) 'Processing  ', JDATE, ':', JTIME
        CALL M3MSG2( ' ' )
        CALL M3MSG2( MESG )

        DO  V = 1, NVARS3D

            IF ( .NOT. READ3( 'BDY_IN', VNAME3D( V ), ALLAYS3, JDATE, JTIME, BUF1 ) ) THEN
                EFLAG = .TRUE.
                CALL M3MESG( 'ERROR:  read-failure' )
                CYCLE
            END IF

            BUF1( NORTH0:WEST1, : ) = NWFAC * BUF1( NORTH0:WEST1, : ) 

            IF ( .NOT.WRITE3( 'BDY_OUT', VNAME3D( V ), JDATE, JTIME, BUF1 ) ) THEN
                EFLAG = .TRUE.
                CALL M3MESG( 'ERROR:  write-failure' )
            END IF

        END DO      !  end loop on variables

        CALL NEXTIME( JDATE, JTIME, TSTEP )

    END DO          !  end loop on output time steps


    IF ( EFLAG ) THEN
        MESG  = 'Failure in program'
        ISTAT = 2
    ELSE
        MESG  = 'Success in program'
        ISTAT = 0
    END IF

    CALL M3EXIT( PNAME, 0, 0, MESG, ISTAT )


END PROGRAM BDYSCALE





