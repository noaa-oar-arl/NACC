
PROGRAM  VERSTAT

    !!***********************************************************************
    !! Version "$Id: VERSTAT.f90 22 2017-09-14 16:31:19Z coats $"
    !! EDSS/Models-3 M3TOOLS.
    !! Copyright (C) 2017 Carlie J. Coats, Jr, and
    !! (C) 2002-2010 Baron Advanced Meteorological Systems, LLC.,
    !! UNC Institute for the Environment
    !! Distributed under the GNU GENERAL PUBLIC LICENSE version 2
    !! See file "GPL.txt" for conditions of use.
    !!.........................................................................
    !!  program body starts at line  103
    !!
    !!  FUNCTION:
    !!
    !!  PRECONDITIONS REQUIRED:
    !!
    !!  REVISION  HISTORY:
    !!      Prototype 10/2017 by CJC
    !!***********************************************************************

    USE M3UTILIO
    IMPLICIT NONE

    !!...........   PARAMETERS and their descriptions:

    CHARACTER*16, PARAMETER :: PNAME = 'VERSTAT'
    CHARACTER*1,  PARAMETER :: BLANK = ' '
    CHARACTER*72, PARAMETER :: BAR   =  &
    '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'

    !!...........   LOCAL VARIABLES and their descriptions:

    INTEGER         LDEV, RDEV, ISTAT
    INTEGER         COL0, COL1, ROW0, ROW1
    INTEGER         C, R, L, N, V

    LOGICAL         EFLAG

    CHARACTER*16    VNAME
    CHARACTER*16    RNAME
    CHARACTER*256   MESG

    !!     GRIDDESC name, parameters for grid

    CHARACTER*16    GDNAM
    INTEGER         FTYPE
    INTEGER         GDTYP
    INTEGER         NCOLS
    INTEGER         NROWS
    INTEGER         NLAYS
    INTEGER         NTHIK
    INTEGER         NVARS
    REAL*8          P_ALP
    REAL*8          P_BET
    REAL*8          P_GAM
    REAL*8          XCENT
    REAL*8          YCENT
    REAL*8          XORIG
    REAL*8          YORIG
    REAL*8          XCELL
    REAL*8          YCELL
    INTEGER         SDATE   ! starting date
    INTEGER         STIME   ! starting time
    INTEGER         EDATE   ! ending date
    INTEGER         ETIME   ! ending time
    INTEGER         JDATE   !  starting date, from user
    INTEGER         JTIME   !  starting time, from user
    INTEGER         TSTEP   !  time step, from INAME header
    INTEGER         NRECS   !  duration in TSTEPs

    REAL   , ALLOCATABLE :: ABUF( :,:,: )
    REAL                 :: AVAL
    REAL                 :: AMAX
    REAL                 :: AMIN
    REAL*8               :: ABAR, DIV
    REAL*8               :: ASIG
    INTEGER              :: RMAX
    INTEGER              :: RMIN
    INTEGER              :: CMAX
    INTEGER              :: CMIN
    INTEGER              :: NMAX
    INTEGER              :: NMIN

    REAL   , ALLOCATABLE :: FAMAX( : )
    REAL   , ALLOCATABLE :: FAMIN( : )
    REAL*8 , ALLOCATABLE :: FABAR( : )
    REAL*8 , ALLOCATABLE :: FASIG( : )
    INTEGER, ALLOCATABLE :: FRMAX( : )
    INTEGER, ALLOCATABLE :: FRMIN( : )
    INTEGER, ALLOCATABLE :: FCMAX( : )
    INTEGER, ALLOCATABLE :: FCMIN( : )
    INTEGER, ALLOCATABLE :: FNMAX( : )
    INTEGER, ALLOCATABLE :: FNMIN( : )

    !!.........................................................................
    !!   begin body of program  VERSTAT

    LDEV = INIT3()
    EFLAG = .FALSE.

    WRITE( *, '( 5X, A )' ) BAR,                                                &
'Program "verstat" to compute per-layer statistics for a "window" into a',      &
'specified variable for the specified time period, and write the result',       &
' to the program-log or to a REPORT file.',                                     &
'',                                                                             &
'',                                                                             &
'PRECONDITIONS REQUIRED:',                                                      &
'    setenv  INFILE    <path name for  input file>',                            &
'    setenv  REPORT    <path name> or "LOG"',                                   &
'',                                                                             &
'    ${INFILE} is of type GRIDDED.',                                            &
'',                                                                             &
'    The requested variable must be of type REAL',                              &
'',                                                                             &
'THE PROGRAM WILL PROMPT YOU for the variable-name, the starting and',          &
'ending date&time for the report period, and the starting and ending',          &
'values for the column and row range being processed.',                         &
'',                                                                             &
'Program copyright (C) 2017 Carlie J. Coats, Jr.',                              &
'Released under Version 2 of the GNU General Public License.',                  &
'See enclosed GPL.txt, or URL',                                                 &
''  ,                                                                           &
'    https://www.gnu.org/licenses/old-licenses/gpl-2.0.html',                   &
BLANK, BAR, ''

    IF ( .NOT. GETVAL( 'Continue with program?', .TRUE. ) ) THEN
        MESG = 'Program terminated at user request'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF


    !!...............  Open files:

    IF ( .NOT.OPEN3( 'INFILE', FSREAD3, PNAME ) ) THEN
        EFLAG = .TRUE.
        MESG  = 'ERROR:  Could not open "INFILE"'
        CALL M3MESG( MESG )
    ELSE IF ( .NOT.DESC3( 'INFILE' ) ) THEN
        EFLAG = .TRUE.
        MESG  = 'ERROR:  Could not get description for "INFILE"'
        CALL M3MESG( MESG )
    ELSE IF ( FTYPE3D .NE. GRDDED3 ) THEN
        EFLAG = .TRUE.
        MESG  = 'ERROR:  unsupported file type for "INFILE"'
        CALL M3MESG( MESG )
    ELSE

        GDNAM = GDNAM3D
        FTYPE = FTYPE3D
        GDTYP = GDTYP3D
        NCOLS = NCOLS3D
        NROWS = NROWS3D
        NLAYS = NLAYS3D
        NTHIK = NTHIK3D
        NVARS = NVARS3D
        GDTYP = GDTYP3D
        P_ALP = P_ALP3D
        P_BET = P_BET3D
        P_GAM = P_GAM3D
        XCENT = XCENT3D
        YCENT = YCENT3D
        XORIG = XORIG3D
        YORIG = YORIG3D
        XCELL = XCELL3D
        YCELL = YCELL3D

    END IF              !  if not.open3(INFILE...); else...

    CALL ENVSTR( 'REPORT', 'GRIDDESC-name for REPORT file', 'LOG', RNAME, ISTAT )
    IF ( ISTAT .GT. 0 ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Bad environment variable "REPORT"', 2 )
    ELSE IF ( RNAME .EQ. 'LOG' ) THEN
        RDEV = LDEV
    ELSE
        RDEV = GETEFILE( 'REPORT', .FALSE., .TRUE., PNAME )
        IF ( RDEV .LT. 0 ) THEN
            CALL M3EXIT( PNAME, 0, 0, 'Failure opening REPORT file', 2 )
        END IF
    END IF

    IF ( EFLAG ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Fatal file-related error(s)', 2 )
    END IF

    ALLOCATE( ABUF( NCOLS, NROWS,NLAYS ),    &
              FAMAX( NLAYS ),                &
              FAMIN( NLAYS ),                &
              FABAR( NLAYS ),                &
              FASIG( NLAYS ),                &
              FRMAX( NLAYS ),                &
              FRMIN( NLAYS ),                &
              FCMAX( NLAYS ),                &
              FCMIN( NLAYS ),                &
              FNMAX( NLAYS ),                &
              FNMIN( NLAYS ),   STAT = ISTAT )
    IF ( ISTAT .NE. 0 ) THEN
        WRITE( MESG, '( A, I10 )' ) 'Buffer allocation failed:  STAT=', ISTAT
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF


    !!...............  Get date&time, run specs:

    CALL RUNSPEC( 'INFILE', .FALSE., SDATE, STIME, TSTEP, NRECS )

    COL0 = GETNUM( 1,    NCOLS,     1, 'Enter first column for window to analyze' )
    COL0 = GETNUM( COL0, NCOLS, NCOLS, 'Enter  last column for window to analyze' )
    ROW0 = GETNUM( 1,    NROWS,     1, 'Enter first  row   for window to analyze' )
    ROW0 = GETNUM( ROW0, NROWS, NROWS, 'Enter  last  row   for window to analyze' )

    CALL M3MESG( 'The list of REAL variables is:' )
    DO V = 1, NVARS
        IF ( VTYPE3D( V ) .EQ. M3REAL ) THEN
            WRITE( MESG, '( I4, ":  ", A, " (", A, "): ", A )' )     &
                V, VNAME3D(V), UNITS3D(V), VDESC3D(V)
        END IF
    END DO
    V = GETNUM( 1, NVARS, 1, 'Enter number for variable to be analyzed' )
    IF ( VTYPE3D( V ) .NE. M3REAL ) THEN
        CALL M3EXIT( PNAME, 0, 0, 'Variable not of type REAL', 2 )
    END IF
    VNAME = VNAME3D( V )

    !!...............   Process time step sequence

    WRITE( RDEV, '( A4, 2( 2X, A14 ), 2( 2X, A14, A8, 1X, A8 ) )' )  &
        'LAY', 'MEAN', 'SIGMA', 'MAX', ' @ COL', 'ROW', 'MIN', ' @ COL', 'ROW'

    
    DIV = 1.0D0 / DBLE( ( ROW1 - ROW0 + 1 )*( COL1 - COL0 + 1 ) )       !!  = 1/window-size

    FAMAX( : ) =  BADVAL3
    FAMIN( : ) = -BADVAL3
    FABAR( : ) =  0.0D0
    FASIG( : ) =  0.0D0
    FRMAX( : ) = -1
    FRMIN( : ) = -1
    FCMAX( : ) = -1
    FCMIN( : ) = -1
    FNMAX( : ) = -1
    FNMIN( : ) = -1

    JDATE = SDATE
    JTIME = STIME
    CALL NEXTIME( JDATE, JTIME, -TSTEP )
    
    DO N = 1, NRECS
    
        CALL NEXTIME( JDATE, JTIME, TSTEP )
        IF ( .NOT.READ3( 'INFILE', VNAME, ALLAYS3, JDATE, JTIME, ABUF ) ) THEN
            EFLAG = .TRUE.
            CYCLE
        END IF

        WRITE( RDEV, '(  1X,  I9.7, ";", I6.6 )' ) JDATE, JTIME

!$OMP   PARALLEL DO                                                     &
!$OMP&      DEFAULT( NONE ),                                            &
!$OMP&       SHARED( NLAYS, ROW0, ROW1, COL0, COL1, ABUF,               &
!$OMP&               FABAR, FASIG, FAMAX, FAMIN, FRMAX, FRMIN,          &
!$OMP&               FCMAX, FCMIN, FNMAX, FNMIN, DIV ),                 &
!$OMP&      PRIVATE( L, R, C, AVAL, ABAR, ASIG, AMAX, AMIN,             &
!$OMP&               RMAX, RMIN, CMAX, CMIN )
        
        DO L = 1, NLAYS

            ABAR = 0.0D0
            ASIG = 0.0D0
            AMAX = ABUF( COL0,ROW0,L )
            AMIN = ABUF( COL0,ROW0,L )
            CMAX = COL0
            CMIN = COL0
            RMAX = ROW0
            RMIN = ROW0

            DO R = ROW0, ROW1
            DO C = COL0, COL1
                AVAL = ABUF( C,R,L )
                ABAR = ABAR + AVAL
                ASIG = ABAR + AVAL**2
                IF ( AVAL .GT. AMAX ) THEN
                    AMAX = AVAL
                    RMAX = R
                    CMAX = C
                END IF
                IF ( AVAL .LT. AMIN ) THEN
                    AMIN = AVAL
                    RMIN = R
                    CMIN = C
                END IF
            END DO
            END DO
            
            FABAR( L ) = FABAR( L ) + ABAR
            FASIG( L ) = FASIG( L ) + ASIG

            IF ( AMAX .GT. FAMAX( L ) ) THEN
                FAMAX( L ) = AMAX
                FCMAX( L ) = CMAX
                FRMAX( L ) = RMAX
                FNMAX( L ) = N
            END IF
            
            IF ( AMIN .LT. FAMIN( L ) ) THEN
                FAMIN( L ) = AMIN
                FCMIN( L ) = CMIN
                FRMIN( L ) = RMIN
                FNMIN( L ) = N
            END IF
            
            ABAR = DIV * ABAR
            ASIG = SQRT( MAX( ASIG * DIV - ABAR**2, 0.0D0 ) )
            WRITE( RDEV, '( I3, ":", 2( 2X, 1PD14.6 ), 2(  2X, 1PE14.6, " @ (", I4, ",", I4, ")" ) )' )     &
                L, ABAR, ASIG, AMAX, CMAX, RMAX, AMIN, CMIN, RMIN

        END DO
        
    END DO

    IF ( EFLAG ) THEN
        MESG  = 'Failure(s) in program'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    DIV = DIV / DBLE( NRECS )
    
    WRITE( RDEV, '(A)' ) BLANK
    DO L = 1, NLAYS
        FABAR( L ) = DIV * FABAR( L )
        FASIG( L ) = SQRT( MAX( DIV * FASIG( L ) - FABAR( L )**2 , 0.0D0 ) )
        WRITE( RDEV, '( A18, 2( 2X, 1PD14.6 ), 2( 2X, 1PE14.6, " @ (", I4, ",", I4, ",", I4, ")" ) )' )                     &
            'File Stats:', FABAR(L), FASIG(L),          &
            FAMAX(L), FCMAX(L), FRMAX(L), FNMAX( L ),   &
            FAMIN(L), FCMIN(L), FRMIN(L), FNMIN( L )
    END DO

    CALL M3EXIT( PNAME, 0, 0, 'Success in program', 0 )
    
END PROGRAM  VERSTAT
