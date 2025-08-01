
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/tmpbeis.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE TMPBEIS312( JDATE, JTIME, TSTEP, N_SPC_BIOG, BEIS_MAP, EMISL )

C-----------------------------------------------------------------------
 
C Description:
C   Computes hourly time stepped gridded biogenic emissions using 
C   normalized gridded emissions from Normbeis3 (3.12) and postprocessed MM5
C   meteorology.
 
C Preconditions:
C   Postprocessed MM5 meteorology that contains temperature, 
C   solar radiation, and pressure data. 
C   Normalized gridded emissions B3GRD from Normbeis3 (3.12) 
 
C Subroutines and Functions Called:
C   BEIS, HRNO, CZANGLE, 
C   ENVYN, CURRSTEP, HHMMSS, NAMEVAL, UPCASE, CHECKMEM, M3EXIT, M3MSG2,
C   WRITE3, OPEN3, CLOSE3, NEXTIME, TIME2SEC
 
C Revision History:
C   3/01: Prototype by Jeff Vukovich
C         Tested only on 36km Lambert domain 
C         Summer/winter switch file option not tested
C   8/04: Updated for BEIS v3.12
C  10/06: J.Young - CMAQ inline
C                 - assumes first call is scenario start date/time
C                 - optional integral average mass diagnostic emissions file
C  10/09: J.Young - intel compiler problem; move last EMISS calculation as
C                   separate assignment
C  01/10: J.Young - move soil NO data and required input data processing to hrno
C  01/10: D. Wong - Eliminate potential race condition at the MYPE = 0
C                   to open the diagnostic file. Eliminate the BARRIER
C                   and the CLOSE3 and OPEN3.

C-----------------------------------------------------------------------
C Modified from:
 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id: tmpbeis.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C smoke@unc.edu
C Pathname: $Source: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/tmpbeis.F,v $
C Last updated: $Date: 2011/02/22 01:10:45 $ 
C-----------------------------------------------------------------------

      USE HGRD_DEFN             ! horizontal domain specifications
      USE BIOG_EMIS             ! beis

      IMPLICIT NONE

C Includes:
      INCLUDE SUBST_IOPARMS    ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC    ! file header data structure
#include      SUBST_IODECL     ! I/O definitions and declarations
        
C Arguments:
      INTEGER, INTENT(  IN ) :: JDATE   ! current simulation date (YYYYDDD)
      INTEGER, INTENT(  IN ) :: JTIME   ! current simulation time (HHMMSS)
      INTEGER, INTENT(  IN ) :: TSTEP( 2 )      ! timestep vector
      INTEGER, INTENT(  IN ) :: N_SPC_BIOG      ! model biogenic species
      INTEGER, INTENT(  IN ) :: BEIS_MAP( : )   ! mapping from beis3 to model species
      REAL,    INTENT( OUT ) :: EMISL( :,:,: )  ! emissions in moles/sec

C External Functions
      INTEGER,         EXTERNAL :: TIME2SEC 
      INTEGER,         EXTERNAL :: ENVINT 
      LOGICAL,         EXTERNAL :: ENVYN
      LOGICAL,         EXTERNAL :: CHKGRID
      CHARACTER( 16 ), EXTERNAL :: PROMPTMFILE
      CHARACTER( 10 ), EXTERNAL :: HHMMSS

C Parameters:
      REAL,    PARAMETER :: HR2SEC = 1.0 / 3600.0

C Latitude and longitude for zenith angle calculation:
      REAL,    ALLOCATABLE, SAVE :: LAT( :,: ) ! grid lat (deg) -90 <= LAT <= 90
      REAL,    ALLOCATABLE, SAVE :: LON( :,: ) ! grid lon (deg) -180 <= LON <= 180 

C Gridded meteorology data:
      REAL,    ALLOCATABLE, SAVE :: TASFC ( :,: )     ! level-1 air temperature (K)
      REAL,    ALLOCATABLE, SAVE :: SSOLAR( :,: )     ! Surface radiation [W/m**2]
      REAL,    ALLOCATABLE, SAVE :: COSZEN( :,: )     ! cosine of zenith angle
      REAL,    ALLOCATABLE, SAVE :: PRES  ( :,: )     ! surface pressure [Pa]
      INTEGER, ALLOCATABLE, SAVE :: SWITCH( :,: )     ! Seasonal switch

C Gridded beis emissions
      REAL,    ALLOCATABLE, SAVE :: AVGEMIS( :,:,:,: ) ! avg emissions (3.12)
      REAL,    ALLOCATABLE, SAVE :: AVGLAI ( :,:,:,: ) ! avg leaf index

C Gridded normalized emissions:
      REAL,    ALLOCATABLE, SAVE :: SEMIS( :,:,: )    ! Normailized Beis emissions
      REAL,    ALLOCATABLE, SAVE :: SLAI ( :,:,: )    ! Normailized Beis LAI

C Mole and mass factors:
      REAL MLFC, MSFC, EMPL

C Diagnostic file avg factor
      REAL DIVFAC

C BEIS3 internal, output species:
      REAL,    ALLOCATABLE, SAVE :: EMPOL( :,:,: )    ! emissions of biogenic categories
C EMPOL( :,:,1:NSEF-1 ) = beis3 emissions
C EMPOL( :,:,NSEF )     = hrno  emissions
      REAL,    ALLOCATABLE, SAVE :: EMISS( :,:,: )    ! emissions in g/s

C Logical names and unit numbers:
            
      CHARACTER( 16 ), SAVE :: SNAME   ! logical name for diagnostic emis output (mass)
      CHARACTER( 16 ), SAVE :: NNAME   ! logical name for normalized-emissions input
      CHARACTER( 16 ), SAVE :: GNAME   ! logical name for GRID_CRO_2D
      CHARACTER( 16 ), SAVE :: BNAME   ! logical name for frost switch input
      CHARACTER( 16 ), SAVE :: MNAME   ! logical name for MET_CRO_2D

C Other variables:
      INTEGER, SAVE :: IHR       ! current simulation hour
      INTEGER          IOS       ! IO or memory allocation status
      INTEGER, SAVE :: LDATE     ! previous simulation date
      INTEGER, SAVE :: MDATE     ! test date to update beis diag avg
      INTEGER, SAVE :: MTIME     ! test time to update beis diag avg
      INTEGER, SAVE :: SDATE     ! scenario start date
      INTEGER, SAVE :: STIME     ! scenario start time
      INTEGER, SAVE :: NSTEPS    ! run duration
      INTEGER, SAVE :: NTICS = 0 ! no. of substeps within an output tstep
      INTEGER          PARTYPE   ! method number to calculate PAR
      INTEGER, SAVE :: LOGDEV
      INTEGER          C, R, I, J, K, S  ! loop counters and subscripts

      LOGICAL          OK        ! check flag
      LOGICAL, SAVE :: SWITCH_FILE   = .TRUE. ! use frost switch file
      LOGICAL       :: ASSUME_SUMMER = .TRUE. ! use summer normalized emissions
      LOGICAL, SAVE :: BEMIS_DIAG             ! true: write diagnostic emiss file
      LOGICAL, SAVE :: LASTTIC                ! true: last sync step this output tstep
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      CHARACTER( 256 ) :: MESG       ! message buffer
      CHARACTER( 256 ) :: EQNAME     ! equivalent filename
      CHARACTER(  16 ) :: VAR        ! variable name
      CHARACTER(  16 ) :: PNAME = 'TMPBEIS312'   ! prodedure name

      INTEGER      GXOFF, GYOFF           ! global origin offset from file
C for INTERPX
      INTEGER       :: STRTCOLNB3, ENDCOLNB3, STRTROWNB3, ENDROWNB3  ! Normalized BEIS3
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2  ! GRID_DRO_2D
      INTEGER, SAVE :: STRTCOLBS,  ENDCOLBS,  STRTROWBS,  ENDROWBS   ! BIOSEASON
      INTEGER, SAVE :: STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2  ! MET_CRO_2D

D     REAL,    ALLOCATABLE, SAVE :: EMISX( : )
D     INTEGER, ALLOCATABLE, SAVE :: CX( : ), RX( : )

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.

         LOGDEV = INIT3()

C Determine last timestamp
         SDATE = ENVINT( 'CTM_STDATE', 'Scenario Start (YYYYJJJ)', 0, IOS )
         STIME = ENVINT( 'CTM_STTIME', 'Scenario Start (HHMMSS)', 0, IOS )
         NSTEPS = ENVINT( 'CTM_RUNLEN', 'Scenario Duration (HHMMSS)', 24, IOS )
         NSTEPS = NSTEPS / TSTEP( 1 )           ! initscen guarantees divisibility

C Check to see if frost date switch file to be used
         MESG = 'Using a frost date switch file?'
         SWITCH_FILE = ENVYN ( 'BIOSW_YN', MESG, .TRUE., IOS )

C Get normalized emissions file, BGRD
         NNAME = PROMPTMFILE( 
     &           'Enter name for Normalized Emissions input file',
     &           FSREAD3, 'B3GRD', PNAME )

C Read description of normalized emissions file
         IF ( .NOT. DESC3( NNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( NNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
         END IF

C Initialize grid definition 
         OK = CHKGRID( NNAME )

C Get domain decomp info for the normalized emissions file
         CALL SUBHFILE ( NNAME, GXOFF, GYOFF,
     &                   STRTCOLNB3, ENDCOLNB3, STRTROWNB3, ENDROWNB3 )

C Open and check bioseason file if using
         IF ( SWITCH_FILE ) THEN

            BNAME = PROMPTMFILE( 
     &              'Enter name for season switch input file',
     &              FSREAD3, 'BIOSEASON', PNAME )
            
C Read description of switch file
            IF ( .NOT. DESC3( BNAME ) ) THEN
               MESG = 'Could not get description of file "' //
     &                 TRIM( BNAME ) // '"'
               CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
            END IF

C Check grid definition 
            IF ( .NOT. CHKGRID( BNAME ) ) THEN
               MESG = 'Grid in file "' // TRIM( BNAME ) //
     &                '" does not match previously set grid.'
               CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
            END IF

C Get domain decomp info for the bioseasons switch file
            CALL SUBHFILE ( BNAME, GXOFF, GYOFF,
     &                      STRTCOLBS, ENDCOLBS, STRTROWBS, ENDROWBS )

            ALLOCATE( SWITCH( NCOLS,NROWS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SWITCH', PNAME )
            SWITCH = 0   ! array

         ELSE

            MESG = 'Use summer normalized emissions?'
            ASSUME_SUMMER = ENVYN ( 'SUMMER_YN', MESG, .TRUE., IOS )

         END IF

C Open met file
         MNAME = PROMPTMFILE( 
     &           'Enter name for gridded met input file',
     &           FSREAD3, 'MET_CRO_2D', PNAME )

C Get description of met file 
         IF ( .NOT. DESC3( MNAME ) ) THEN
            MESG = 'Could not get description of file "'
     &           // TRIM( MNAME ) // '"'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
         END IF

C Check that grid description matches BGRD file
         IF ( .NOT. CHKGRID( MNAME ) ) THEN
            MESG = 'Grid in file "' // TRIM( MNAME )
     &           // '" does not match previously set grid.'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
         END IF

C Get domain decomp info for the met file
         CALL SUBHFILE ( MNAME, GXOFF, GYOFF,
     &                   STRTCOLMC2, ENDCOLMC2, STRTROWMC2, ENDROWMC2 )

C Write emissions diagnostic file?
         BEMIS_DIAG = ENVYN( 'B3GTS_DIAG', 'Write emissions diagnostic file?',
     &                       .FALSE., IOS )
           
         IF ( BEMIS_DIAG ) THEN
C Build description for, and open beis diagnostic file
C (all but variables-table and horizontal domain in description is borrowed from MNAME)
            MDATE = SDATE; MTIME = STIME
            CALL NEXTIME( MDATE, MTIME, TSTEP( 1 ) )
            SDATE3D = MDATE
            STIME3D = MTIME
            TSTEP3D = TSTEP( 1 )
            XORIG3D = XORIG_GD
            YORIG3D = YORIG_GD
            NCOLS3D = GL_NCOLS
            NROWS3D = GL_NROWS
            NLAYS3D = 1
!           NVARS3D = MSPCS
            NVARS3D = N_SPC_BIOG

            J = 0
            DO I = 1, MSPCS
               IF ( BEIS_MAP( I ) .GT. 0 ) THEN
                  J = J + 1
                  VNAME3D( J ) = EMSPC( I )
                  VDESC3D( J ) = 'biogenic emissions of the indicated species'
                  VTYPE3D( J ) = M3REAL
                  UNITS3D( J ) = 'gm/s'
               END IF
            END DO

            FDESC3D = ' '   ! array

            FDESC3D( 1 ) = 'Gridded biogenic emissions from CMAQ-BEIS3'
            FDESC3D( 2 ) = '/from/ ' // PNAME
            FDESC3D( 3 ) = '/Version/ CMAQ'

C Open mass output file (tons/hour)

            SNAME = 'B3GTS_S'
            IF ( MYPE .EQ. 0 ) THEN
               IF ( .NOT. OPEN3( SNAME, FSNEW3, PNAME ) ) THEN
                  MESG = 'Could not open "' // TRIM( SNAME ) // '" file'
                  CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT1 )
               END IF
            END IF

         END IF

C Build name table for variables in normalized emissions file
         ALLOCATE( AVGEMIS( NCOLS,NROWS,NSEF-1,NSEASONS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'AVGEMIS', PNAME )

         ALLOCATE( AVGLAI( NCOLS,NROWS,NLAI,NSEASONS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'AVGLAI', PNAME )

         AVGEMIS = 0.0 ! array
         AVGLAI  = 0.0 ! array

C Open 2-D grid parameters file to get LAT and LON
         GNAME = PROMPTMFILE( 
     &           'Enter name for 2D GRID PARAMETERS input file',
     &           FSREAD3, 'GRID_CRO_2D', PNAME )

         IF ( .NOT. DESC3( GNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( GNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
         END IF

C Check grid description against BGRD File 
         IF ( .NOT. CHKGRID( GNAME ) ) THEN
            MESG = 'Grid in file "' // TRIM( GNAME ) //
     &             '" does not match previously set grid.'
            CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
         END IF

C Get domain decomp info from the normalized emissions file (assume met files have
C same grid dimensions)
         CALL SUBHFILE ( NNAME, GXOFF, GYOFF,
     &                   STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2 )

C Allocate memory for data and read
         ALLOCATE( LAT( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'LAT', PNAME )

         ALLOCATE( LON( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'LON', PNAME )

         ALLOCATE( COSZEN( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'COSZEN', PNAME )

         IF ( .NOT. INTERPX( GNAME, 'LAT', PNAME,
     &                       STRTCOLGC2,ENDCOLGC2, STRTROWGC2,ENDROWGC2, 1,1,
     &                       0, 0, LAT ) ) THEN
            MESG = 'Could not read LAT from file "' //
     &              TRIM( GNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
         END IF

         IF ( .NOT. INTERPX( GNAME, 'LON', PNAME,
     &                       STRTCOLGC2,ENDCOLGC2, STRTROWGC2,ENDROWGC2, 1,1,
     &                       0, 0, LON ) ) THEN
            MESG = 'Could not read LON from file "' //
     &              TRIM( GNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
         END IF

C Read the various categories of normalized emissions
         DO I = 1, NSEASONS 
             
            DO J = 1, NSEF-1
               VAR = 'AVG_' // TRIM( BIOTYPES( J ) ) // SEASON( I )
         
               IF ( .NOT. XTRACT3( NNAME, VAR,
     &                             1,1, STRTROWNB3,ENDROWNB3, STRTCOLNB3,ENDCOLNB3,
     &                             0, 0, AVGEMIS( 1,1,J,I ) ) ) THEN
                  MESG = 'Could not read "' // TRIM( VAR ) //
     &                   '" from file "' // TRIM( NNAME ) // '"'
                  CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
               END IF
            END DO

            DO K = 1, NLAI
               VAR = 'LAI_' // TRIM( LAITYPES( K ) ) // SEASON( I )
 
               IF ( .NOT. XTRACT3( NNAME, VAR,
     &                             1,1, STRTROWNB3,ENDROWNB3, STRTCOLNB3,ENDCOLNB3,
     &                             0, 0, AVGLAI( 1,1,K,I ) ) ) THEN
                  MESG = 'Could not read "' // TRIM( VAR ) //
     &                   '" from file "' // TRIM( NNAME ) // '"'
                  CALL M3EXIT( PNAME, 0, 0, MESG, XSTAT2 )
               END IF
            END DO

         END DO   ! end loop over seasons

C Allocate memory for arrays
         ALLOCATE( TASFC( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'TASFC', PNAME )
         
         ALLOCATE( SSOLAR( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'SSOLAR', PNAME )

         ALLOCATE( EMPOL( NCOLS,NROWS,NSEF ), STAT=IOS )
         CALL CHECKMEM( IOS, 'EMPOL', PNAME )

         ALLOCATE( PRES( NCOLS,NROWS ), STAT=IOS )
         CALL CHECKMEM( IOS, 'PRES', PNAME )

         IF ( BEMIS_DIAG ) THEN
            ALLOCATE( EMISS( NCOLS,NROWS,N_SPC_BIOG ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMISS', PNAME )
         END IF

         ALLOCATE( SEMIS( NCOLS,NROWS,NSEF-1 ), STAT=IOS )
         CALL CHECKMEM( IOS, 'SEMIS', PNAME )

         ALLOCATE( SLAI( NCOLS,NROWS,NLAI ), STAT=IOS )
         CALL CHECKMEM( IOS, 'SLAI', PNAME )

C Initialize normalized emissons to be used 
         IF ( ASSUME_SUMMER ) THEN
            SEMIS = AVGEMIS( 1:MY_NCOLS,1:MY_NROWS,1:NSEF-1,NSUMMER )
            SLAI  = AVGLAI ( 1:MY_NCOLS,1:MY_NROWS,1:NLAI,  NSUMMER )
         ELSE
            SEMIS = AVGEMIS( 1:MY_NCOLS,1:MY_NROWS,1:NSEF-1,NWINTER )
            SLAI  = AVGLAI ( 1:MY_NCOLS,1:MY_NROWS,1:NLAI,  NWINTER )
         END IF

         LDATE = 0

D        ALLOCATE( EMISX( MSPCS ), STAT=IOS )
D        CALL CHECKMEM( IOS, 'EMISX', PNAME )
D        ALLOCATE( CX( MSPCS ), STAT=IOS )
D        CALL CHECKMEM( IOS, 'CX', PNAME )
D        ALLOCATE( RX( MSPCS ), STAT=IOS )
D        CALL CHECKMEM( IOS, 'RX', PNAME )

         IF ( BEMIS_DIAG ) EMISS = 0.0   ! array

      END IF   ! FIRSTIME

C Loop thru the number of time steps (hourly)

      EMISL = 0.0   ! array
      EMPOL = 0.0   ! array
D     EMISX = 0.0   ! array

C Compute zenith angle
      CALL CZANGLE( JDATE, JTIME, MY_NCOLS, MY_NROWS, LAT, LON, COSZEN ) 

      IF ( JDATE .NE. LDATE ) THEN

         CALL WRDAYMSG( JDATE, MESG )               

C If new date, read season switch 
         IF ( SWITCH_FILE ) THEN
            MESG = 'Reading gridded season switch data...'
            CALL M3MSG2( MESG ) 
             
            IF ( .NOT. XTRACT3( BNAME, 'SEASON',
     &                          1,1, STRTROWBS,ENDROWBS, STRTCOLBS,ENDCOLBS,
     &                          JDATE, 0, SWITCH ) ) THEN
               MESG = 'Could not read SEASON from file "' //
     &                TRIM( BNAME ) // '"'
               CALL M3EXIT( PNAME, JDATE, 0, MESG, XSTAT2 )
            END IF

            MESG = 'Applying gridded season switch data...' 
            CALL M3MSG2( MESG )

            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  IF ( SWITCH( C,R ) .EQ. 0 ) THEN   ! use winter normalized emissions
                     SEMIS( C,R,1:NSEF-1 ) = AVGEMIS( C,R,1:NSEF-1,NWINTER )
                     SLAI ( C,R,1:NLAI   ) =  AVGLAI( C,R,1:NLAI,  NWINTER )
                  ELSE
                     SEMIS( C,R,1:NSEF-1 ) = AVGEMIS( C,R,1:NSEF-1,NSUMMER )
                     SLAI ( C,R,1:NLAI   ) =  AVGLAI( C,R,1:NLAI,  NSUMMER )
                  END IF                      
               END DO
            END DO
 
         END IF   ! if using switch file

      END IF   ! if new day

      WRITE( MESG,94030 ) HHMMSS( JTIME )
      CALL M3MSG2( MESG )

C Read temperature data
!     IF ( .NOT. INTERPX( MNAME, 'TEMP10', PNAME,
      IF ( .NOT. INTERPX( MNAME, 'TEMP2', PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    JDATE, JTIME, TASFC ) ) THEN
!        MESG = 'Could not read "' // 'TEMP10' // 
         MESG = 'Could not read "' // 'TEMP2' // 
     &          '" from file "' // TRIM( MNAME ) // '"'
         CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
      END IF

C Read surface radiation data
      IF ( .NOT. INTERPX( MNAME, 'RGRND', PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    JDATE, JTIME, SSOLAR ) ) THEN
         MESG = 'Could not read "' // 'RGRND' // 
     &          '" from file "' // TRIM( MNAME ) // '"' 
         CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
      END IF

C Read surface pressure data 
      IF ( .NOT. INTERPX( MNAME, 'PRSFC', PNAME,
     &                    STRTCOLMC2,ENDCOLMC2, STRTROWMC2,ENDROWMC2, 1,1,
     &                    JDATE, JTIME, PRES ) ) THEN
         MESG = 'Could not read "' // 'PRSFC' // 
     &           '" from file "' // TRIM( MNAME ) // '"'
         CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
      END IF

C Calculate temporal non-speciated beis emissions
      CALL BEIS3( JDATE, JTIME, MY_NCOLS, MY_NROWS, MSPCS, COSZEN, SEMIS,
     &            SLAI, TASFC, SSOLAR, PRES, EMPOL )

C Calculate temporal non-speciated soil NO emissions
      CALL HRNO( JDATE, JTIME, TSTEP, TASFC, EMPOL )

C Speciate emissions

D     EMISX = 0.0; CX = 0; RX = 0    ! arrays
      DO K = 1, NSEF
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               EMPL = EMPOL( C,R,K ) * HR2SEC
               IF ( EMPL .LE. 0.0 ) CYCLE
               I = 0
               DO S = 1, MSPCS
                  IF ( BEIS_MAP( S ) .GT. 0 ) THEN
                     I = I + 1
                     MLFC = MLFAC( S,K )
                     EMISL( I,C,R ) = EMISL( I,C,R ) + EMPL * MLFC
D                    IF ( EMISL( I,C,R ) .GT. EMISX( I ) ) THEN
D                        EMISX( I ) = EMISL( I,C,R )
D                        CX( I ) = C
D                        RX( I ) = R
D                    END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

D     WRITE( LOGDEV,* ) ' '
D     WRITE( LOGDEV,* ) '       from TMPBEIS312'
D     WRITE( LOGDEV,* ) '    Spc   EMISL  MaxC MaxR'
D     DO S = 1, N_SPC_BIOG
D        WRITE( LOGDEV,'( 4X, I3, F10.5, 2I4 )' ) S, EMISX( S ), CX( S ), RX( S )
D     END DO

D     EMISX = 0.0    ! array
D     WRITE( LOGDEV,* ) ' '
D     WRITE( LOGDEV,* ) '                   from TMPBEIS312'
D     WRITE( LOGDEV,* ) '    Spc SpcName  Pol  MLFAC     EMPOL     EMISL( Spc,32,7 )'
D     C = 32; R = 7
D     DO K = 1, NSEF
D        EMPL = EMPOL( C,R,K ) * HR2SEC
D        IF ( EMPL .LE. 0.0 ) CYCLE
D        I = 0
D        DO S = 1, MSPCS
D           IF ( BEIS_MAP( S ) .GT. 0 ) THEN
D              I = I + 1
D              MLFC = MLFAC( S,K )
D              EMISX( I ) = EMISX( I ) + EMPL * MLFC
D              IF ( MLFC .NE. 0.0 )
D    &            WRITE( LOGDEV,'( 4X, I3, 2X, A6, 2X, I3, 1X, 3( 1PE10.3 ) )' )
D    &                              S, EMSPC( I ), K,  MLFC, EMPL, EMISX( I )
D           END IF
D        END DO
D     END DO

      NTICS = NTICS + 1
      LASTTIC = NTICS .EQ. TIME2SEC( TSTEP( 1 ) ) / TIME2SEC( TSTEP( 2 ) )

      IF ( BEMIS_DIAG ) THEN
         IF ( NTICS .EQ. 1 ) THEN
            DO K = 1, NSEF
               I = 0
               DO S = 1, MSPCS
                  IF ( BEIS_MAP( S ) .GT. 0 ) THEN
                     I = I + 1
                     MSFC = HR2SEC * MSFAC( S,K )
                     EMISS( :,:,I ) = EMISS( :,:,I ) + EMPOL( :,:,K ) * MSFC
                  END IF
               END DO
            END DO
         ELSE IF ( .NOT. LASTTIC ) THEN
            DO K = 1, NSEF
               I = 0
               DO S = 1, MSPCS
                  IF ( BEIS_MAP( S ) .GT. 0 ) THEN
                     I = I + 1
                     MSFC = HR2SEC * MSFAC( S,K )
                     EMISS( :,:,I ) = EMISS( :,:,I ) + 2.0 * EMPOL( :,:,K ) * MSFC
                  END IF
               END DO
            END DO
         ELSE   ! LASTTIC
            DIVFAC = 0.5 / FLOAT( NTICS )
            DO K = 1, NSEF
               I = 0
               DO S = 1, MSPCS
                  IF ( BEIS_MAP( S ) .GT. 0 ) THEN
                     I = I + 1
                     MSFC = HR2SEC * MSFAC( S,K )
                     EMISS( :,:,I ) = EMISS( :,:,I ) + EMPOL( :,:,K ) * MSFC
                  END IF
               END DO
            END DO
            EMISS = DIVFAC * EMISS   ! array assignment
            IF ( .NOT. WRITE3( SNAME, 'ALL', MDATE, MTIME, EMISS ) ) THEN
               MESG = 'Could not write to output file "' // TRIM( SNAME ) // '"'
               CALL M3EXIT( PNAME, JDATE, JTIME, MESG, XSTAT2 )
            END IF
            WRITE( LOGDEV,94040 )
     &            'Timestep written to', SNAME,
     &            'for date and time', MDATE, MTIME
            EMISS = 0.0   ! array
            CALL NEXTIME( MDATE, MTIME, TSTEP( 1 ) )
         END IF
      END IF

C Save date
      LDATE = JDATE

      IF ( IHR .LT. NSTEPS ) THEN
         IF ( LASTTIC ) THEN
            NTICS = 0
         END IF
         RETURN
      ELSE IF ( NTICS .LT. TIME2SEC( TSTEP( 1 ) ) / TIME2SEC( TSTEP( 2 ) ) ) THEN
         RETURN
      END IF

C------------------  FORMAT  STATEMENTS   ------------------------------

94010 FORMAT( 10( A, ':', I8, ':', 1X ) )
94020 FORMAT( 1X, I7, ":", I6.7 )
94030 FORMAT( 5X, 'Temporal BEIS at time ', A8 )
94040 FORMAT( /5X, 3( A, :, 1X ), I8, ":", I6.6 )

      END SUBROUTINE TMPBEIS312  

