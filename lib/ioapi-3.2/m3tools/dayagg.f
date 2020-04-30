
        PROGRAM DAYAGG

C***********************************************************************
C Version "$Id: dayagg.f 1 2017-06-10 18:05:20Z coats $"
C EDSS/Models-3 M3TOOLS.
C Copyright (C) 1992-2002 MCNC, (C) 1995-2002,2005-2013 Carlie J. Coats, Jr.,
C and (C) 2002-2010 Baron Advanced Meteorological Systems. LLC.
C Distributed under the GNU GENERAL PUBLIC LICENSE version 2
C See file "GPL.txt" for conditions of use.
C.........................................................................
C  program body starts at line  101
C
C  DESCRIPTION:
C       Program DAYAGG to construct "seasonal standard day" weighted
C       averages from a set of input files and write the result to a
C       "standard day" output file.
C
C  PRECONDITIONS REQUIRED:
C       All files have the same grid structure and contain all the
C       requested output variables on the requested time step structure.
C       Requested variables are all of type REAL
C       setenv STTIME   <starting time for the "seasonal standard day"
C                       in format HHMMSS>
C       setenv TSTEP    <time step for the "seasonal standard day"
C                       in format HHMMSS>
C       setenv RUNLEN   <duration for the "seasonal standard day"
C                       in format HHMMSS>
C       setenv FILELIST <comma-delimited list of input files>
C       setenv WGTSLIST <comma-delimited list of weighting factors
C                       for the input-files>
C       setenv VBLELIST <comma-delimited list of output variables>
C       setenv DATELIST <comma-delimited list of starting dates>
C       setenv AGGFILE  <physical (path) name>
C       For each of the input files,
C               setenv <logical name> <physical (path) name>
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 12/6/2000 by Carlie J. Coats, Jr., NCSC, adapted
C       from "m3tproc'
C
C       Version   11/2001 by CJC for I/O API Version 2.1
C
C       Version   11/2005 by CJC:  eliminate unused vbles
C
C       Version   06/2008 by CJC:  output starting date is
C         1 = 1000*YEAR + DAY
C
C       Version 02/2010 by CJC for I/O API v3.1:  Fortran-90 only;
C       USE M3UTILIO, and related changes.
C
C       Version 10/2014 by CJC:  Fix handling of ENV*() calls.
C***********************************************************************

      USE M3UTILIO

      IMPLICIT NONE

C...........   PARAMETERS and their descriptions:

        CHARACTER*16, PARAMETER :: PNAME   = 'DAYAGG'

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER     STTIME
        INTEGER     TSTEP
        INTEGER     RUNLEN
        INTEGER     NFILES, NSTEPS, NVARS
        INTEGER     JDATE, JTIME
        LOGICAL     EFLAG

        CHARACTER*16  VNAME( MXVARS3 )
        CHARACTER*16  UNITS( MXVARS3 )
        CHARACTER*80  VDESC( MXVARS3 )

        CHARACTER*16  FNAME( MXFILE3 )
        INTEGER      SDATES( MXFILE3 )
        INTEGER      JDATES( MXFILE3 )
        INTEGER      JTIMES( MXFILE3 )
        REAL        WEIGHTS( MXFILE3 )

        REAL,    ALLOCATABLE::   ABUF( :, :, : )
        REAL,    ALLOCATABLE::   TBUF( :, :, : )

        REAL            WTFAC

        INTEGER         ESTAT
        INTEGER         I, C, R, L, V, N, T

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

        CHARACTER*256   MESG


C***********************************************************************
C   begin body of program DAYAGG

        EFLAG = .FALSE.
        I = INIT3()
        WRITE( *,'( 5X, A )' )
     &  ' ',
     &  'Program DAYAGG to construct a "standard day" file for a set',
     &  'of selected variables, by constructing the weighted average',
     &  'of the variables from a set of input files, and write the',
     &  'result to a "standard day" output file.',
     &  ' ',
     &  'PRECONDITIONS REQUIRED: ',
     &  '     All files have the same grid structure and contain all',
     &  '     the requested output variables on the requested time',
     &  '     step structure.',
     &  '     Requested variables are all of type REAL',
     &  '     setenv STTIME   <starting time (HHMMSS) for the output>',
     &  '     setenv TSTEP    <time step (HHMMSS) for the output>',
     &  '     setenv RUNLEN   <duration (HHMMSS) for the output>',
     &  '     setenv AGGFILE  <output-file physical (path) name>',
     &  '     setenv FILELIST <(comma-delimited) list of input file',
     &  '                      logical names>',
     &  '     setenv WGTSLIST <list of weighting factors>',
     &  '     setenv VBLELIST <list of output variables>',
     &  '     setenv DATELIST <list of starting dates>',
     &  ' For each of the input files,',
     &  '     setenv <logical name> <physical name>',
     &' ',
     &'See URL',
     &'https://www.cmascenter.org/ioapi/documentation/3.1/html#tools',
     &' ',
     &'Program copyright (C) 1992-2002 MCNC, (C) 1995-2013',
     &'Carlie J. Coats, Jr., and (C) 2002-2010 Baron Advanced',
     &'Meteorological Systems, LLC.  Released under Version 2',
     &'of the GNU General Public License. See enclosed GPL.txt, or',
     &'URL http://www.gnu.org/copyleft/gpl.html',
     &' ',
     &'See URL  http://www.baronams.com/products/ioapi/AA.html#tools',
     &'Comments and questions are welcome and can be sent to',
     &' ',
     &'    Carlie J. Coats, Jr.    carlie@jyarborough.com',
     &'or',
     &'    UNC Institute for the Environment',
     &'    137 E. Franklin St. Suite 602 Room 613-C',
     &'    Campus Box 1105',
     &'    Chapel Hill, NC 27599-1105',
     &' ',
     &'Program version: ',
     &'$Id:: dayagg.f 1 2017-06-10 18:05:20Z coats                   $',
     &' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( PNAME, 0, 0,
     &                   'Program terminated at user request', 2 )
        END IF

        STTIME = ENVINT( 'STTIME',
     &                   'Starting TIME (HHMMSS)', 0, ESTAT )
        IF ( ESTAT .GT. 0 ) THEN
            MESG  = 'Invalid/missing environment variable "STTIME"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        TSTEP = ENVINT( 'TSTEP',
     &                  'TIME-STEP (HHMMSS)', 10000, ESTAT )
        IF ( ESTAT .NE. 0 ) THEN
            MESG  = 'Invalid/missing environment variable "TSTEP"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        RUNLEN = ENVINT( 'RUNLEN',
     &                  'Run duration (HHMMSS)', 0, ESTAT )
        IF ( ESTAT .NE. 0 ) THEN
            MESG  = 'Invalid/missing environment variable "RUNLEN"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        IF ( .NOT. STRLIST( 'VBLELIST',
     &                      'List of extracted output variables',
     &                      MXVARS3, NVARS, VNAME ) ) THEN
            MESG  = 'Bad environment vble "VBLELIST"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        IF ( .NOT. STRLIST( 'FILELIST',
     &                      'List of input-file logical names',
     &                      MXFILE3, NFILES, FNAME ) ) THEN
            MESG  = 'Bad environment vble "FILELIST"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        IF ( .NOT. INTLIST( 'DATELIST',
     &                      'List of input-file starting dates',
     &                      MXFILE3, N, SDATES ) ) THEN
            MESG  = 'Bad environment vble "DATELIST"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        ELSE IF ( N .NE. NFILES ) THEN
            WRITE ( MESG, '( A, I4, 2X, A, I4 )' )
     &         'Inconsistency:  # of input files:', NFILES,
     &         '# of DATES:', N
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        IF ( .NOT. REALIST( 'WGTSLIST',
     &                      'List of input-file weighting factors',
     &                      MXFILE3, N, WEIGHTS ) ) THEN
            MESG  = 'Bad environment vble "WGTSLIST"'
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        ELSE IF ( N .NE. NFILES ) THEN
            WRITE ( MESG, '( A, I4, 2X, A, I4 )' )
     &         'Inconsistency:  # of input files:', NFILES,
     &         '# of WEIGHTS:', N
            EFLAG = .TRUE.
            CALL M3MESG( MESG )
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( PNAME, 0,0, 'Bad environment', 2 )
        END IF


C...............  Open input files

        N = 1
        IF ( .NOT. OPEN3( FNAME( N ), FSREAD3, PNAME ) ) THEN
            MESG = 'Could not open ' // FNAME( N )
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ELSE IF ( .NOT. DESC3( FNAME( N ) ) ) THEN
            MESG = 'Could not get description for ' // FNAME( N )
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

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
        DO  V = 1, NVARS
            I = INDEX1( VNAME( V ), NVARS3D, VNAME3D )
            IF ( 0 .GE. I ) THEN
               MESG = 'Variable "' // TRIM( VNAME( V ) )  //
     &                '" not available in file ' // FNAME( N )
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            ELSE IF ( VTYPE3D( I ) .NE. M3REAL ) THEN
                MESG = 'Variable "' // TRIM( VNAME( V ) )  //
     &                '" not of type REAL'
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF
            UNITS( V ) = UNITS3D( I )
            VDESC( V ) = VDESC3D( I )

        END DO

        DO  N = 2, NFILES       !  open/check the rest of the input files

            IF ( .NOT. OPEN3( FNAME(N), FSREAD3, PNAME ) ) THEN
                MESG = 'Could not open ' // FNAME(N)
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            ELSE IF ( .NOT. DESC3( FNAME(N) ) ) THEN
                MESG = 'Could not get description for ' // FNAME(N)
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            ELSE IF ( .NOT.FILCHK3( FNAME(N),  GRDDED3,
     &                      NCOLS1, NROWS1, NLAYS1, NTHIK3D ) ) THEN
                MESG = 'Inconsistent dimensions  for ' // FNAME(N)
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            ELSE IF ( .NOT.GRDCHK3( FNAME(N),
     &                      P_ALP1, P_BET1, P_GAM1, XCENT1, YCENT1,
     &                      XORIG1, YORIG1, XCELL1, YCELL1,
     &                      NLAYS1, VGTYP3D, VGTOP3D, VGLVS3D ) ) THEN
                MESG = 'Inconsistent coord/grid  for ' // FNAME(N)
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

            DO  V = 1, NVARS
                IF ( 0 .GE. INDEX1( VNAME( V ),
     &                              NVARS3D, VNAME3D ) ) THEN
                   MESG = 'Variable "' // TRIM( VNAME( V ) )  //
     &                '" not available in file ' // FNAME( N )
                    CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
                END IF
            END DO

        END DO          !  end loop opening the rest of the input files


C...............  Open the output file

        SDATE3D = 1
        STIME3D = STTIME
        TSTEP3D = TSTEP
        NVARS3D = NVARS

        DO  V = 1, NVARS
            VNAME3D( V ) = VNAME( V )
            UNITS3D( V ) =  UNITS( V )
            VDESC3D( V ) =  VDESC( V )
            VTYPE3D( V ) =  M3REAL
        END DO

        IF ( .NOT. OPEN3( 'AGGFILE', FSUNKN3, PNAME ) ) THEN
            MESG = 'Could not open AGGFILE'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF


C...............  File-related setup:

        WTFAC = 0.0             !  for normalizing the weights
        DO  N = 1, NFILES
            WTFAC = WTFAC + WEIGHTS( N )
        END DO
        WTFAC  = 1.0 / WTFAC
        DO  N = 1, NFILES
            WEIGHTS( N ) = WTFAC * WEIGHTS( N )
        END DO

        JDATE = 1
        JTIME = STTIME
        DO  N = 1, NFILES
            JDATES( N ) = SDATES( N )
            JTIMES( N ) = STTIME
        END DO

        NSTEPS = TIME2SEC( RUNLEN ) / TIME2SEC( TSTEP )

        ALLOCATE( ABUF( NCOLS1, NROWS1, NLAYS1 ),
     &            TBUF( NCOLS1, NROWS1, NLAYS1 ), STAT = ESTAT )
        IF ( ESTAT .NE. 0 ) THEN
            WRITE( MESG, '( A, I10)' )
     &               'Buffer allocation failed:  STAT=', ESTAT
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF


C...............  Process the requested time period and variables:

        DO  T = 1, NSTEPS

            DO  V = 1, NVARS

                !!  initialize with first input file:

                N = 1
                IF ( .NOT. READ3( FNAME(N), VNAME(V), ALLAYS3,
     &                            JDATES(N), JTIMES(N), TBUF ) ) THEN
                    MESG = 'Could not read "' // TRIM( VNAME(V) )//
     &                     '" from file ' // FNAME(N)
                    CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
                END IF

                DO  L = 1, NLAYS1
                DO  R = 1, NROWS1
                DO  C = 1, NCOLS1
                    TBUF( C,R,L ) = TBUF( C,R,L ) * WEIGHTS( N )
                END DO
                END DO
                END DO

                DO  N = 2, NFILES       !  rest of the input files

                    IF ( .NOT. READ3( FNAME(N), VNAME(V), ALLAYS3,
     &                                JDATES(N), JTIMES(N),
     &                                ABUF ) ) THEN
                        MESG = 'Could not read "'//TRIM( VNAME(V) )//
     &                         '" from file ' // FNAME(N)
                        CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                    DO  L = 1, NLAYS1
                    DO  R = 1, NROWS1
                    DO  C = 1, NCOLS1
                        TBUF( C,R,L )  =  TBUF( C,R,L )
     &                                  + ABUF( C,R,L ) * WEIGHTS( N )
                    END DO
                    END DO
                    END DO

                END DO  !  end loop on the rest of the input files

                IF ( .NOT. WRITE3( 'AGGFILE', VNAME(V),
     &                             JDATE, JTIME, TBUF )) THEN
                    MESG = 'Could not write "' // TRIM( VNAME(V) )//
     &                     '" to file AGGFILE'
                    CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
                END IF

            END DO      !  end loop on variables

            CALL NEXTIME( JDATE, JTIME, TSTEP )
            DO  N = 1, NFILES
                CALL NEXTIME( JDATES( N ), JTIMES( N ), TSTEP )
            END DO

        END DO          !  end loop on time steps

        CALL M3EXIT( PNAME, 0, 0,
     &               'Successful completion of program DAYAGG', 0 )

        END PROGRAM DAYAGG

