
C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/init/yamo/load_cgrid.F,v 1.2 2011/08/11 09:56:01 jaroslav Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE LOAD_CGRID ( FNAME, JDATE, JTIME, SPC_CAT, CMIN, CGRID )

C-----------------------------------------------------------------------
C Function:
C   Initialize the model CGRID array from file data

C Revision history:
C   Jeff - Dec 00 - split out from initscen.F
C                 - move CGRID_MAP into f90 module
C   Jeff - Feb 01 - module LOCALFILE to set CCTM IOAPI file defn's
C                 - add BUFSIZE for dfio option
C                 - assumed shape arrays
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX;
C   30 Oct 01 J.Young: fix ICBC_FAC
C    4 Sep 03 J.Young: fix SPC/NDX bug if ASO4J IC's are unavailable
C   20 Nov 03 J.Young: enable loading RHOJ
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_ICBC     ! gas chem ic/bc surrogate names and map table
      INCLUDE SUBST_AE_ICBC     ! aerosol ic/bc surrogate names and map table
      INCLUDE SUBST_NR_ICBC     ! non-react ic/bc surrogate names and map table
      INCLUDE SUBST_TR_ICBC     ! tracer ic/bc surrogate names and map table
      INCLUDE SUBST_GC_G2AE     ! gas chem aerosol species and map table
      INCLUDE SUBST_CONST       ! constants
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: FNAME
      INTEGER      JDATE
      INTEGER      JTIME
      CHARACTER( 2 ) ::  SPC_CAT
      REAL         CMIN
!     REAL      :: CGRID( :,:,:,: )  ! for initial CONC
      REAL, POINTER :: CGRID( :,:,:,: )  ! for initial CONC

C Parameters:

C minimum aerosol sulfate concentration [ ug/m**3 ]
      REAL, PARAMETER :: AEROCONCMIN = 0.001

C The following two factors assume that sulfate density is 1.8e3 [ kg/m**3 ]
C and that the geometric mean diameter and geometric standard deviations
C for the Aitken mode are 0.01e-6 [ m ] and 1.7 respectively
C and are 0.07e-6 and 2.0 respectively for the accumulation mode.

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Aitken mode [ ug ].
      REAL, PARAMETER :: NUMFACT_I = 2.988524E11

C factor to calculate aerosol number concentration from aerosol sulfate mass
C concentration in the Accumulation mode [ ug ].
      REAL, PARAMETER :: NUMFACT_J = 3.560191E08

C fraction of sulfuric acid vapor taken as aerosol for first time step
      REAL, PARAMETER :: SO4VAPTOAER = 0.999
C initial fraction of total aerosol sulfate in the Aitken mode
      REAL, PARAMETER :: IFRACATKN = 0.04

      INTEGER, PARAMETER :: MXSPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC

C File variables:

      REAL         DENS( NCOLS,NROWS,NLAYS )       ! air density (kg/m^3)
      REAL         RHOJ( NCOLS,NROWS,NLAYS ) ! air density X Jacobian (kg/m^2)

C External Functions (not already declared by IODECL3.EXT):

      INTEGER, EXTERNAL :: FINDEX       !  looks up number in table.
      INTEGER, EXTERNAL :: INDEX1       !  looks up name in table.
      INTEGER, EXTERNAL :: TRIMLEN      !  string length, excl. trailing blanks

C Local Variables

      REAL         MWH2SO4                           ! H2SO4 molec. wt.
      REAL         H2SO4CONV                         ! ppm -> ug/m**3
      INTEGER      LSULF                             ! Gas chem CGRID index
      INTEGER      ISO4AJ, ISO4AI, INUMATKN, INUMACC ! CGRID aerosol indices

      INTEGER      GXOFF, GYOFF               ! global origin offset from file

C for INTERPX
!     INTEGER :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C for XTRACT3
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER          STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI
      REAL         DBUFF( NCOLS,NROWS,NLAYS )

      INTEGER      SPC_STRT
      INTEGER      N_SPCS
      INTEGER      INDX( MXSPC ), ISUR        ! Variable and surrogate indices
      INTEGER      NDX                        ! loop copy of INDX
      REAL         ICBC_FAC( MXSPC )          ! Factor to be applied to ICs
      INTEGER      C, R, L, SPC, V            ! loop counters
      INTEGER      ASPC                       ! CGRID RHOJ pointer
      INTEGER      NCOLSDENS, NROWSDENS       ! local for DENS

      CHARACTER( 16 ) :: PNAME = 'LOAD_CGRID'
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 16 ) :: CONCMIN
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 40 ) :: CHWARN = 'Domain extents different from model for '

      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER, SAVE :: LOGDEV

C-----------------------------------------------------------------------
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )
C IOFDESC common now loaded with MET_CRO_3D header
         END IF
      WRITE( CONCMIN,'(1PE8.2)' ) CMIN
      IF ( .NOT. OPEN3( FNAME, FSREAD3, PNAME ) ) THEN
         XMSG = 'Could not open ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      IF ( .NOT. DESC3( FNAME ) ) THEN
         XMSG = 'Could not get ' // FNAME // ' file description'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
 
      IF ( GL_NCOLS .NE. NCOLS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NCOLS: ', GL_NCOLS, '  NCOLS3D: ', NCOLS3D
         END IF
 
      IF ( GL_NROWS .NE. NROWS3D ) THEN
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    WARNING: ' // CHWARN // FNAME
         WRITE( LOGDEV,* ) '>>  GL_NROWS: ', GL_NROWS, '  NROWS3D: ', NROWS3D
         END IF
 
      IF ( NLAYS .NE. NLAYS3D ) THEN
         XMSG = 'Wrong number of layers in ' // FNAME // ' file'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get CGRID offsets

      CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C Get INDX
      IF ( SPC_CAT .EQ. 'GC' ) THEN
         SPC_STRT = GC_STRT
         N_SPCS = N_GC_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( GC_SPC( SPC ), NVARS3D, VNAME3D )
            ICBC_FAC( SPC ) = 1.0
            IF ( NDX .EQ. 0 ) THEN
               ISUR = FINDEX ( SPC , N_GC_ICBC, GC_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  NDX = INDEX1( GC_ICBC ( ISUR ), NVARS3D, VNAME3D )
                  ICBC_FAC( SPC ) = GC_ICBC_FAC( ISUR )
                  END IF
               END IF
            INDX( SPC ) = NDX
            IF ( NDX .EQ. 0 ) THEN
               XMSG = 'No IC found for species ' //
     &                GC_SPC( SPC )( 1:TRIMLEN( GC_SPC( SPC ) ) )  //
     &                ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                ' Set to ' // CONCMIN( 1:TRIMLEN( CONCMIN ) ) 
               CALL M3MESG( XMSG )
               END IF
            END DO
         ELSE IF ( SPC_CAT .EQ. 'AE' ) THEN
         SPC_STRT = AE_STRT
         N_SPCS = N_AE_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( AE_SPC( SPC ), NVARS3D, VNAME3D )
            ICBC_FAC( SPC ) = 1.0
            IF ( NDX .EQ. 0 ) THEN
               ISUR = FINDEX ( SPC , N_AE_ICBC, AE_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  NDX = INDEX1( AE_ICBC ( ISUR ), NVARS3D, VNAME3D )
                  ICBC_FAC( SPC ) = AE_ICBC_FAC( ISUR )
                  END IF
               END IF
            INDX( SPC ) = NDX
            IF ( NDX .EQ. 0 ) THEN
               XMSG = 'No IC found for species ' //
     &                AE_SPC( SPC )( 1:TRIMLEN( AE_SPC( SPC ) ) )  //
     &                ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                ' Set to ' // CONCMIN( 1:TRIMLEN( CONCMIN ) ) 
               CALL M3MESG( XMSG )
               END IF
            END DO
         ELSE IF ( SPC_CAT .EQ. 'NR' ) THEN
         SPC_STRT = NR_STRT
         N_SPCS = N_NR_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( NR_SPC( SPC ), NVARS3D, VNAME3D )
            ICBC_FAC( SPC ) = 1.0
            IF ( NDX .EQ. 0 ) THEN
               ISUR = FINDEX ( SPC , N_NR_ICBC, NR_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  NDX = INDEX1( NR_ICBC ( ISUR ), NVARS3D, VNAME3D )
                  ICBC_FAC( SPC ) = NR_ICBC_FAC( ISUR )
                  END IF
               END IF
            INDX( SPC ) = NDX
            IF ( NDX .EQ. 0 ) THEN
               XMSG = 'No IC found for species ' //
     &                NR_SPC( SPC )( 1:TRIMLEN( NR_SPC( SPC ) ) )  //
     &                ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                ' Set to ' // CONCMIN( 1:TRIMLEN( CONCMIN ) ) 
               CALL M3MESG( XMSG )
               END IF
            END DO
         ELSE IF ( SPC_CAT .EQ. 'TR' ) THEN
         SPC_STRT = TR_STRT
         N_SPCS = N_TR_SPC
         DO SPC = 1, N_SPCS
            NDX = INDEX1( TR_SPC( SPC ), NVARS3D, VNAME3D )
            ICBC_FAC( SPC ) = 1.0
            IF ( NDX .EQ. 0 ) THEN
               ISUR = FINDEX ( SPC , N_TR_ICBC, TR_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  NDX = INDEX1( TR_ICBC ( ISUR ), NVARS3D, VNAME3D )
                  ICBC_FAC( SPC ) = TR_ICBC_FAC( ISUR )
                  END IF
               END IF
            INDX( SPC ) = NDX
            IF ( NDX .EQ. 0 ) THEN
               XMSG = 'No IC found for species ' //
     &                TR_SPC( SPC )( 1:TRIMLEN( TR_SPC( SPC ) ) )  //
     &                ' in ' // FNAME( 1:TRIMLEN( FNAME ) ) //
     &                ' Set to ' // CONCMIN( 1:TRIMLEN( CONCMIN ) ) 
               CALL M3MESG( XMSG )
               END IF
            END DO
         ELSE IF ( SPC_CAT .EQ. 'RJ' ) THEN
         N_SPCS = 0
!        WRITE( LOGDEV,'(/ 5X, A)' ) 'loading Density*Jacobian into CGRID'
         ELSE
         XMSG = 'Species categories incorrect for CGRID '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF
        
C Read into CGRID

      CALL SUBHFILE ( FNAME, GXOFF, GYOFF,
     &                STRTCOLINI, ENDCOLINI, STRTROWINI, ENDROWINI )
C IOFDESC common now loaded with FNAME header

      DO SPC = 1, N_SPCS
         V = SPC_STRT - 1 + SPC
         NDX = INDX( SPC )

         IF ( NDX .GT. 0 ) THEN
            IF ( .NOT. XTRACT3( FNAME, VNAME3D ( NDX ),
     &           1,NLAYS, STRTROWINI,ENDROWINI, STRTCOLINI,ENDCOLINI,
     &           JDATE, JTIME, DBUFF ) ) THEN
               XMSG = 'Could not read ' //
     &                VNAME3D( NDX )( 1:TRIMLEN( VNAME3D( NDX ) ) ) //
     &                ' from ' // FNAME
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,V ) = ICBC_FAC( SPC ) * DBUFF( C,R,L )
                     END DO
                  END DO
               END DO

            END IF   ! INDX .GT. 0

         END DO

      IF ( N_SPCS .NE. 0 ) WRITE( LOGDEV,'(/ 5X, A)' )
     &                            SPC_CAT // ' loaded into CGRID'

      IF ( SPC_CAT .EQ. 'RJ' ) THEN

C Load RHOJ for transport and mixing ratio advection adjustment

         VNAME = 'DENSA_J'
         IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME,
     &              1,NLAYS, STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
     &              JDATE, JTIME, RHOJ ) ) THEN
             XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         ASPC = GC_STRT - 1 + N_GC_SPCD
!        DO R = 1, MY_NROWS
!           DO C = 1, MY_NCOLS
!              DO L = 1, NLAYS
!                 CGRID( ASPC,L,C,R ) = RHOJ( C,R,L )
         DO L = 1, NLAYS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  CGRID( C,R,L,ASPC ) = RHOJ( C,R,L )
                  END DO
               END DO
            END DO

         WRITE( LOGDEV,'(/ 5X, A)' ) 'Density*Jacobian loaded into CGRID'

         END IF

      IF ( SPC_CAT .EQ. 'AE' ) THEN

C are ASO4J IC's available on the file?

         VNAME = 'ASO4J'
         NDX = INDEX1( VNAME, NVARS3D, VNAME3D )
     
         IF ( NDX .EQ. 0 ) THEN  ! ASO4J not on file

C Set pointers for gas (vapor) phase sulfur species

            NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
            IF ( NDX .NE. 0 ) THEN
               ISO4AJ = AE_STRT - 1 + NDX
               ELSE
               XMSG = 'Could not find ' // VNAME // 'in aerosol table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

            VNAME = 'SULF'
            NDX = INDEX1( VNAME, N_GC_G2AE, GC_G2AE )
            IF ( NDX .NE. 0 ) THEN
               LSULF   = GC_STRT - 1 + GC_G2AE_MAP( NDX )
               MWH2SO4 = GC_MOLWT( GC_G2AE_MAP( NDX ) )
               ELSE
               XMSG = 'Could not find ' // VNAME // 'in gas chem aerosol table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

            VNAME = 'ASO4I'
            NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
            IF ( NDX .NE. 0 ) THEN
               ISO4AI = AE_STRT - 1 + NDX
               ELSE
               XMSG = 'Could not find ' // VNAME // 'in aerosol table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            VNAME = 'NUMATKN'
            NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
            IF ( NDX .NE. 0 ) THEN
               INUMATKN = AE_STRT - 1 + NDX
               ELSE
               XMSG = 'Could not find ' // VNAME // 'in aerosol table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            VNAME = 'NUMACC'
            NDX = INDEX1( VNAME, N_AE_SPC, AE_SPC )
            IF ( NDX .NE. 0 ) THEN
               INUMACC = AE_STRT - 1 + NDX
               ELSE
               XMSG = 'Could not find ' // VNAME // 'in aerosol table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

            NCOLSDENS = ENDCOLMC3 - STRTCOLMC3 + 1
            NROWSDENS = ENDROWMC3 - STRTROWMC3 + 1

            IF ( NCOLSDENS .NE. MY_NCOLS .OR.
     &           NROWSDENS .NE. MY_NROWS ) THEN
               WRITE( XMSG,'( A, 4I8 )' ) 'Local Columns or Rows incorrect',
     &         NCOLSDENS, MY_NCOLS, NROWSDENS, MY_NROWS
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

C IOFDESC now loaded with MET_CRO_3D header

            VNAME = 'DENS'
!           IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
!    &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
!    &                          1,NLAYS,
!    &                          JDATE, JTIME, DENS ) ) THEN
            IF ( .NOT. XTRACT3( MET_CRO_3D, VNAME, 1,NLAYS,
     &                          STRTROWMC3,ENDROWMC3, STRTCOLMC3,ENDCOLMC3,
     &                          JDATE, JTIME, DENS ) ) THEN
               XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

C Partition the aerosol sulfate arrays with a fraction of the initial SO4 

            H2SO4CONV = 1.0E3 * MWH2SO4 / MWAIR * SO4VAPTOAER

            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS

C total accumulation mode sulfate:

                     CGRID( C,R,L,ISO4AJ )   = MAX ( AEROCONCMIN,
     &                                         ( 1.0 - IFRACATKN )
     &                                       * H2SO4CONV
     &                                       * DENS ( C,R,L )
     &                                       * CGRID( C,R,L,LSULF ) )

C Accumulation mode number:
    
                     CGRID( C,R,L,INUMACC )  = NUMFACT_J
     &                                       * CGRID( C,R,L,ISO4AJ )

C Aitken mode sulfate:
    
                     CGRID( C,R,L,ISO4AI )   = MAX ( AEROCONCMIN,
     &                                         IFRACATKN
     &                                       * H2SO4CONV
     &                                       * DENS ( C,R,L )
     &                                       * CGRID( C,R,L,LSULF ) )
    
C Aitken mode number:
    
                     CGRID( C,R,L,INUMATKN ) = NUMFACT_I
     &                                       * CGRID( C,R,L,ISO4AI )
    
C correct sulfate vapor concentration for part removed:
    
                     CGRID( C,R,L,LSULF )    = ( 1.0 - SO4VAPTOAER )
     &                                       * CGRID( C,R,L,LSULF)
    
                     END DO
                  END DO
               END DO

            XMSG = 'No IC''s found for aerosol sulfate. ' //
     &             'Gas Chem sulfate used for partitioning.'
            CALL M3MESG( XMSG )

            END IF  ! NDX .EQ. 0

         END IF  !  SPC_CAT .EQ. 'AE'

C Close the file

!      IF ( .NOT. CLOSE3( FNAME ) ) THEN
!         XMSG = 'Could not close ' // FNAME // ' file'
!         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

      RETURN
      END
