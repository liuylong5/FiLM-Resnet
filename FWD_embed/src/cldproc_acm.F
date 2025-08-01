
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
C $Header: /project/work/rep/CCTM/src/cloud/cloud_acm_ae5/cldproc_acm.F,v 1.5 2009/12/09 17:40:03 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CLDPROC ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C
C  FUNCTION:  RADM/ACM and Resolved cloud process driver
C
C  PRECONDITIONS REQUIRED:
C       Dates and times represented YYYYDDD:HHMMSS.
C
C  REVISION  HISTORY:
C      Date   Who             What
C    -------- ---             -----------------------------------------
C     11/21/00 J. Young       enable appending timesteps to wetdep1 and wetdep2
C     3/01/98 S.Roselle       modified to output 2 wet deposition files
C     8/12/97 S.Roselle       added conversion of cgrid units to mol/mol units
C     6/14/94 Dongming Hwang  configuration management prototype
C     Dec 00  J. Young        move CGRID_MAP into f90 module
C     Sep 01  J. Young        Dyn Alloc - Use HGRD_DEFN
C     Jan 05  J. Young        dyn alloc - establish both horizontal & vertical
C                             domain specifications in one module
C     May 05  J. Pleim        Replaced RADMcld with RADMacmcld
C     6/08/05 S.Roselle       added new cloud diagnostic variables
C     7/21/09 D. Wong         fixed parallel bug in opening WET_DEP_2
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications
#ifndef Incl
      USE CGRID_SPCS          ! CGRID mechanism species
#else
      USE CGRID_SPCS          ! CGRID species number and offsets
#endif

      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex
      USE CHK_FILES_DEFS
      USE CHK_FILES

      IMPLICIT NONE

C...........INCLUDES

      INCLUDE SUBST_IOPARMS          ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC          ! file header data structure
#include      SUBST_IODECL           ! I/O definitions and declarations
#ifdef Incl
      INCLUDE SUBST_GC_SPC           ! gas chemistry species table
      INCLUDE SUBST_GC_WDEP          ! wet deposition table for gases
      INCLUDE SUBST_AE_SPC           ! aerosol species table
      INCLUDE SUBST_AE_WDEP          ! wet deposition table for aerosols
      INCLUDE SUBST_NR_SPC           ! non-reactive species table
      INCLUDE SUBST_NR_WDEP          ! wet deposition table for non-reactives
      INCLUDE SUBST_TR_SPC           ! tracer species table
      INCLUDE SUBST_TR_WDEP          ! wet deposition table for tracers
#endif

      INCLUDE SUBST_CONST            ! constants
      INCLUDE SUBST_FILES_ID         ! file name parameters

      CHARACTER( 120 ) :: XMSG = ' ' ! exit status message string

C...........PARAMETERS

C # of wet deposition species
#ifndef Incl
      INTEGER, SAVE :: N_SPC_WDEP
#else
      INTEGER, PARAMETER :: N_SPC_WDEP = N_GC_WDEP
     &                                 + N_AE_WDEP
     &                                 + N_NR_WDEP
     &                                 + N_TR_WDEP
      INTEGER, PARAMETER :: N_SPC_WDEPD = N_SPC_WDEP + 1
#endif

      REAL, PARAMETER :: CNV1 = MWAIR * 1.0E-9
      REAL, PARAMETER :: CNV1I = 1.0 / CNV1
      REAL, PARAMETER :: CNV2 = MWAIR * 1.0E-3
      REAL, PARAMETER :: CNV2I = 1.0 / CNV2
!     REAL, PARAMETER :: CNV3 = MWAIR * 1.0E+3 / AVO  ! -> ppmV
!     REAL, PARAMETER :: CNV3 = CNV2 / AVO            ! -> mol/mol
      REAL, PARAMETER :: CNV3 = CNV2                  ! -> #/mol
      REAL, PARAMETER :: CNV3I = 1.0 / CNV3

C...........ARGUMENTS

!     REAL          CGRID( NCOLS, NROWS, NLAYS, * )  ! concentrations
!      REAL       :: CGRID( :,:,:,: )                 ! concentrations
      REAL, POINTER :: CGRID( :,:,:,: )                 ! concentrations
      INTEGER       JDATE            ! current model date, coded YYYYDDD
      INTEGER       JTIME            ! current model time, coded HHMMSS
      INTEGER       TSTEP( 2 )       ! model time step, coded HHMMSS

C...........Local Variables

      LOGICAL, SAVE :: FIRSTIME = .TRUE. ! flag for first pass thru

      LOGICAL, SAVE :: CLD_DIAG          ! flag to output cloud diagnostic files

      CHARACTER( 16 ), SAVE :: PNAME = 'CLDPROC' ! driver program name
      CHARACTER( 16 ) :: VNAME            ! input variable name list

      INTEGER       COL              ! column subscript indices
      INTEGER       FINI             ! ending position
      INTEGER       ICLDTYPE         ! 1: raining, 2: either CNP or PFW
      INTEGER       L                ! loop counter
      INTEGER       LAY              ! layer subscript indices
      INTEGER       MDATE            ! middle of this time step
      INTEGER       MTIME            ! middle of this time step
      INTEGER       NDATE            ! middle of this time step
      INTEGER       NTIME            ! middle of this time step
      INTEGER, SAVE :: NNAE          ! number of #/m3 species
      INTEGER, SAVE :: NQAE          ! number of ug/m3 species
      INTEGER, SAVE :: NSAE          ! number of m2/m3 species
      INTEGER       ROW              ! row subscript indices
      INTEGER       SPC              ! species subscript indices
      INTEGER       STRT             ! starting position
      INTEGER       STATUS           !  ENVINT status
      INTEGER       VAR              ! variable subscript indices
      INTEGER       ALLOCSTAT
#ifndef Incl
      INTEGER, ALLOCATABLE, SAVE :: WDEP_MAP( : ) ! wet deposition map to CGRID
      INTEGER, ALLOCATABLE, SAVE :: QAE( : ) ! CGRID pointer to ug/m3 species
      INTEGER, ALLOCATABLE, SAVE :: NAE( : ) ! CGRID pointer to #/m3 species
      INTEGER, ALLOCATABLE, SAVE :: SAE( : ) ! CGRID pointer to m2/m3 species
#else
      INTEGER, SAVE :: WDEP_MAP( N_SPC_WDEPD ) ! wet deposition map to CGRID
      INTEGER, SAVE :: QAE( N_AE_SPCD ) ! CGRID pointer to ug/m3 species
      INTEGER, SAVE :: NAE( N_AE_SPCD ) ! CGRID pointer to #/m3 species
      INTEGER, SAVE :: SAE( N_AE_SPCD ) ! CGRID pointer to m2/m3 species
#endif

      INTEGER, SAVE :: LOGDEV           ! output log unit number


!     REAL          DENS    ( NCOLS,NROWS,NLAYS ) ! air density (kg/m3)
!     REAL          CONV_DEP( NCOLS,NROWS,N_SPC_WDEP+8 ) ! convective wdep only
!     REAL          TOT_DEP ( NCOLS,NROWS,N_SPC_WDEP+1 ) ! total wdep
      REAL, ALLOCATABLE, SAVE :: DENS    ( :,:,: ) ! air density (kg/m3)
      REAL, ALLOCATABLE, SAVE :: CONV_DEP( :,:,: ) ! convective wdep only
      REAL, ALLOCATABLE, SAVE :: TOT_DEP ( :,:,: ) ! total wdep

      REAL FAC                       ! temp conversion factor
      REAL CCMIN

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C...........EXTERNAL FUNCTIONS

      LOGICAL, EXTERNAL :: ENVYN

      INTERFACE
         SUBROUTINE RESCLD ( CGRID, JDATE, JTIME, TSTEP,
     &                       N_SPC_WDEP, WDEP_MAP, DEP )
            IMPLICIT NONE
            REAL, POINTER          :: CGRID( :,:,:,: )
!            REAL :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )  :: JDATE, JTIME, TSTEP( 2 )
            INTEGER, INTENT( IN )  :: N_SPC_WDEP
            INTEGER, INTENT( IN )  :: WDEP_MAP( : )
            REAL,    INTENT( OUT ) :: DEP( :,:,: )
         END SUBROUTINE RESCLD
         SUBROUTINE CONVCLD_ACM ( CGRID, JDATE, JTIME, TSTEP,
     &                            N_SPC_WDEP, WDEP_MAP, DEP )
            IMPLICIT NONE
            REAL, POINTER          :: CGRID( :,:,:,: )
!            REAL :: CGRID( :,:,:,: )
            INTEGER  :: JDATE, JTIME, TSTEP( 2 )
            INTEGER  :: N_SPC_WDEP
            INTEGER  :: WDEP_MAP( : )
            REAL :: DEP( :,:,: )
         END SUBROUTINE CONVCLD_ACM
      END INTERFACE
C-----------------------------------------------------------------------
C   begin body of subroutine  CLDPROC

C...Initialization

      IF ( FIRSTIME ) THEN
        LOGDEV = INIT3()
        FIRSTIME = .FALSE.

C...first check to make sure that some species in CGRID were specified
C...for output in the wet deposition array, otherwise notify the user
C...and return

#ifndef Incl
        N_SPC_WDEP = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + N_TR_WDEP
        ALLOCATE ( WDEP_MAP( N_SPC_WDEP ), STAT = ALLOCSTAT )
        IF ( ALLOCSTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating WDEP_MAP'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF
#endif
        IF ( N_SPC_WDEP .LE. 0 ) THEN

          XMSG = 'No species were specified for wet deposition ' //
     &           'tracking'
          CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )

          XMSG = 'ONLY CLOUD DIAGNOSTICS WILL BE WRITTEN TO THE ' //
     &           'WETDEP OUTPUT FILE!'
          CALL M3MESG ( XMSG )

        END IF

#ifdef Incl
        CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )
#endif

C...check to see if user wants to output extra diagnostic files

        CLD_DIAG = .FALSE.
        XMSG = 'Output cloud diagnostic files? (Y/N)'
        CLD_DIAG = ENVYN( 'CLD_DIAG', XMSG, CLD_DIAG, STATUS )

        IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) XMSG

        IF ( STATUS .EQ. 1 ) THEN
          XMSG = 'Environment variable improperly formatted'
          CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
        ELSE IF ( STATUS .EQ. -1 ) THEN
          XMSG = 'Environment variable set, but empty ... Using default:'
          WRITE( LOGDEV, '(5X, A, L9)' ) XMSG, CLD_DIAG
        ELSE IF ( STATUS .EQ. -2 ) THEN
          XMSG = 'Environment variable not set ... Using default:'
          WRITE( LOGDEV, '(5X, A, L9)' ) XMSG, CLD_DIAG
        END IF

        SPC = 0
        STRT = 1
        FINI = N_GC_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = GC_STRT - 1 + GC_WDEP_MAP( SPC )
        END DO

        SPC = 0
        STRT = N_GC_WDEP + 1
        FINI = N_GC_WDEP + N_AE_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = AE_STRT - 1 + AE_WDEP_MAP( SPC )
        END DO

        SPC = 0
        STRT = N_GC_WDEP + N_AE_WDEP + 1
        FINI = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = NR_STRT - 1 + NR_WDEP_MAP( SPC )
        END DO

        SPC = 0
        STRT = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + 1
        FINI = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + N_TR_WDEP
        DO VAR = STRT, FINI
          SPC = SPC + 1
          WDEP_MAP( VAR ) = TR_STRT - 1 + TR_WDEP_MAP( SPC )
        END DO

C Try to open existing wet deposition file for update

        IF ( MYPE .EQ. 0 ) THEN

          IF ( .NOT. OPEN3( CTM_WET_DEP_1, FSRDWR3, PNAME ) ) THEN

            XMSG = 'Could not open ' // CTM_WET_DEP_1 // ' file for update - '
     &           // 'try to open new'
            CALL M3MESG( XMSG )

C...Create wet deposition file:
C...   copy most of header description from
C...   MET_CRO_2D and initialize depositions to zero.

C...open MET_CRO_2D

            IF ( .NOT. OPEN3( MET_CRO_2D, FSREAD3, PNAME ) ) THEN
              XMSG = 'Could not open '// MET_CRO_2D // ' file'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

C...get description from the met file

            IF ( .NOT. DESC3( MET_CRO_2D ) ) THEN
              XMSG = 'Could not get ' // MET_CRO_2D //' file description'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            NCOLS3D = GL_NCOLS
            NROWS3D = GL_NROWS

            XORIG3D = XORIG_GD
            YORIG3D = YORIG_GD

C...advance to time variable for the output file to the next hour

            SDATE3D = JDATE
            STIME3D = 10000 * ( JTIME / 10000 )     !  clear MMSS fields
            TSTEP3D = 010000   ! 1-hour timesteps (hhmmss)
            CALL NEXTIME ( SDATE3D, STIME3D, TSTEP3D ) !  on the next hour

            NVARS3D = N_SPC_WDEP + 1

C...define wet deposition output descriptors for gaseous species

            SPC = 0
            STRT = 1
            FINI = N_GC_WDEP
            DO VAR = STRT, FINI
              SPC = SPC + 1
              VTYPE3D ( VAR ) = M3REAL
              VNAME3D ( VAR ) = GC_WDEP( SPC )
              UNITS3D ( VAR ) = 'kg/hectare'
              VDESC3D ( VAR ) = 'hourly wet deposition values'
            END DO

C...define wet deposition output descriptors for aerosol species

            SPC = 0
            STRT = N_GC_WDEP + 1
            FINI = N_GC_WDEP + N_AE_WDEP
            DO VAR = STRT, FINI
              SPC = SPC + 1
              VTYPE3D ( VAR ) = M3REAL
              VNAME3D ( VAR ) = AE_WDEP( SPC )
              IF ( INDEX( AE_WDEP( SPC ), 'NUM' ) .GT. 0 ) THEN
                UNITS3D ( VAR ) = 'number/hectare'
              ELSE IF ( INDEX( AE_WDEP( SPC ), 'SRF' ) .GT. 0 ) THEN
                UNITS3D ( VAR ) = 'm2/hectare'
              ELSE
                UNITS3D ( VAR ) = 'kg/hectare'
              END IF
              VDESC3D ( VAR ) = 'hourly wet deposition values'
            END DO

C...define wet deposition output descriptors for non-reactive species

            SPC = 0
            STRT = N_GC_WDEP + N_AE_WDEP + 1
            FINI = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP
            DO VAR = STRT, FINI
              SPC = SPC + 1
              VTYPE3D ( VAR ) = M3REAL
              VNAME3D ( VAR ) = NR_WDEP( SPC )
              UNITS3D ( VAR ) = 'kg/hectare'
              VDESC3D ( VAR ) = 'hourly wet deposition values'
            END DO

C...define wet deposition output descriptors for tracer species

            SPC = 0
            STRT = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + 1
            FINI = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + N_TR_WDEP
            DO VAR = STRT, FINI
              SPC = SPC + 1
              VTYPE3D ( VAR ) = M3REAL
              VNAME3D ( VAR ) = TR_WDEP( SPC )
              UNITS3D ( VAR ) = 'kg/hectare'
              VDESC3D ( VAR ) = 'hourly wet deposition values'
            END DO

C...define output descriptor for hydrogen concentration

            VAR = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + N_TR_WDEP + 1
            VTYPE3D( VAR ) = M3REAL
            VNAME3D( VAR ) = 'HPLUS'
            UNITS3D( VAR ) = 'kg/hectare'
            VDESC3D( VAR ) = 'hourly wet deposition values'

C...initialize the file descriptors

            FDESC3D( 1 ) = 'hourly 1-layer cross-point wet deposition '
     &                       // 'from all cloud types'
            DO L = 2 , MXDESC3
              FDESC3D( L ) = ' '
            END DO

            GDNAM3D = GRID_NAME  ! from HGRD_DEFN

C...try to open the output file new

            IF ( .NOT. OPEN3( CTM_WET_DEP_1, FSNEW3, PNAME ) ) THEN
              XMSG = 'Could not create '// CTM_WET_DEP_1 // ' file'
              CALL M3EXIT ( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
            END IF

          END IF   ! open CTM_WET_DEP_1 read/write

C...open the diagnostic file if user asked for it

          IF ( CLD_DIAG ) THEN

C...Try to open existing file for update

            IF ( .NOT. OPEN3( CTM_WET_DEP_2, FSRDWR3, PNAME ) ) THEN

              XMSG = 'Could not open ' // CTM_WET_DEP_2 // ' file for update - '
     &             // 'try to open new'
              CALL M3MESG( XMSG )

C...copy most of header description from
C...CTM_WET_DEP_1 and initialize depositions to zero.

C...get description from the other wet dep file

              IF ( .NOT. DESC3( CTM_WET_DEP_1 ) ) THEN
                XMSG = 'Could not get ' // CTM_WET_DEP_1 //' file description'
                CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
              END IF

C...advance to time variable for the output file to the next hour

              SDATE3D = JDATE
              STIME3D = 10000 * ( JTIME / 10000 )     !  clear MMSS fields
              TSTEP3D = 010000   ! 1-hour timesteps (hhmmss)
              CALL NEXTIME ( SDATE3D, STIME3D, TSTEP3D ) !  on the next hour

              NVARS3D = N_SPC_WDEP + 8

C...define output descriptor for pressure at the lifting condensation
C...level (LCL)

              VAR = N_GC_WDEP + N_AE_WDEP + N_NR_WDEP + N_TR_WDEP + 2
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'PLCL'
              UNITS3D( VAR ) = 'Pa'
              VDESC3D( VAR ) = 'pressure at lifting condensation level, '
     &                         // 'or MISSING = -9999.9'

C...define output descriptor for cloud bottom layer number

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'CLOUD_BOTTOM'
              UNITS3D( VAR ) = 'layer-number'
              VDESC3D( VAR ) =
     &              'layer containing bottom of cloud, or MISSING = -9999.9'

C...define output descriptor for precipitating cloud top layer number

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'PR_CLOUDTOP'
              UNITS3D( VAR ) = 'layer-number'
              VDESC3D( VAR ) =
     &              'layer containing top of PR cloud, or MISSING = -9999.9'

C...define output descriptor for non-precipitating cloud top layer number

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'NP_CLOUDTOP'
              UNITS3D( VAR ) = 'layer-number'
              VDESC3D( VAR ) =
     &              'layer containing top of NP cloud, or MISSING = -9999.9'

C...define output descriptor for rain event flag

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'RAIN_FLAG'
              UNITS3D( VAR ) = '1 or 0'
              VDESC3D( VAR ) = 'Rain-event flag'
 
C...define output descriptor for precitating cloud fraction

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'PR_CLDFRAC'
              UNITS3D( VAR ) = ' '
              VDESC3D( VAR ) = 'PR cloud fraction'

C...define output descriptor for non-precipitating cloud fraction

              VAR = VAR + 1
              VTYPE3D( VAR ) = M3REAL
              VNAME3D( VAR ) = 'NP_CLDFRAC'
              UNITS3D( VAR ) = ' '
              VDESC3D( VAR ) = 'NP cloud fraction'

C...initialize the file descriptors

              FDESC3D( 1 ) = 'hourly 1-layer cross-point wet deposition '
     &                           // 'from sub-grid clouds and cloud data'
              DO L = 2 , MXDESC3
                FDESC3D( L ) = ' '
              END DO

C...try to open the diagnostic output file

              IF ( .NOT. OPEN3( CTM_WET_DEP_2, FSNEW3, PNAME ) ) THEN
                XMSG = 'Could not create '// CTM_WET_DEP_2 // ' file'
                CALL M3EXIT ( PNAME, SDATE3D, STIME3D, XMSG, XSTAT1 )
              END IF

            END IF   ! open CTM_WET_DEP_2 read/write

          END IF   ! CLD_DIAG

        END IF   ! MYPE = 0

        IF ( N_AE_SPC .GT. 0 ) THEN
C...create aerosol species pointers to distinguish micro-grams / m**3
C...  # / m**3 (number density), and m**2 / m**3 (surface area) units

#ifndef Incl
             ALLOCATE ( QAE( N_AE_SPC ),
     &                  NAE( N_AE_SPC ),
     &                  SAE( N_AE_SPC ), STAT = ALLOCSTAT )
             IF ( ALLOCSTAT .NE. 0 ) THEN
                XMSG = 'Failure allocating QAE, NAE, or SAE'
                CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
             END IF
#endif

           NQAE = 0       ! number of ug/m3 species
           NNAE = 0       ! number of #/m3 species
           NSAE = 0       ! number of m2/m3 species

           DO VAR = 1, N_AE_SPC
             IF ( AE_SPC( VAR )( 1:3 ) .EQ. 'NUM' ) THEN
               NNAE = NNAE + 1
               NAE( NNAE ) = AE_STRT - 1 + VAR
             ELSE IF ( AE_SPC( VAR )( 1:3 ) .EQ. 'SRF' ) THEN
               NSAE = NSAE + 1
               SAE( NSAE ) = AE_STRT - 1 + VAR
             ELSE
               NQAE = NQAE + 1
               QAE( NQAE ) = AE_STRT - 1 + VAR
             END IF
           END DO

        END IF

C...initialize the deposition array before processing clouds

        ALLOCATE ( TOT_DEP ( MY_NCOLS,MY_NROWS,N_SPC_WDEP+1 ),
     &             STAT = ALLOCSTAT )
        IF ( ALLOCSTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating TOT_DEP'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        TOT_DEP = 0.0

        ALLOCATE ( CONV_DEP( MY_NCOLS,MY_NROWS,N_SPC_WDEP+8 ),
     &             STAT = ALLOCSTAT )
        IF ( ALLOCSTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating CONV_DEP'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        CONV_DEP = 0.0

        ALLOCATE ( DENS( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
        IF ( ALLOCSTAT .NE. 0 ) THEN
           XMSG = 'Failure allocating DENS'
           CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                  STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C Create new adjoint file and open for update
C        CALL CHK_FILE_OPEN ( ADJ_CLD_CHK, JDATE, JTIME, .TRUE., TSTEP )
      END IF   ! FIRSTIME

      MDATE = JDATE
      MTIME = JTIME

      VNAME = 'DENS'

!     IF ( .NOT. INTERP3( MET_CRO_3D, VNAME, PNAME, MDATE, MTIME,
!    &                    NCOLS * NROWS * NLAYS, DENS ) ) THEN
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, DENS ) ) THEN
        XMSG = 'Could not interpolate DENS'
        CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      END IF

C Checkpoint concentrations
C      CALL CHK_FILE_WRITE ( ADJ_CLD_CHK, JDATE, JTIME, CGRID )

C...convert units on cgrid to molar mixing ratio (by volume)
C...  determine no. of moles of gas per moles of air by volume

      STRT = GC_STRT
      FINI = GC_STRT + N_GC_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = MAX(CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E-6,  1.0E-30)
            END DO
          END DO
        END DO
      END DO
	
C...for aerosol mass concentration
C...  convert to moles of aerosol per mole of air by volume

      STRT = 1
      FINI = NQAE
      DO VAR = STRT, FINI
        SPC = QAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              FAC = CNV1 / DENS( COL, ROW, LAY ) / AE_MOLWT( SPC - AE_STRT + 1 )
              CCMIN = MAX(CGRID( COL, ROW, LAY, SPC ), 1.0E-30/FAC)
              CGRID( COL, ROW, LAY, SPC ) = FAC * CCMIN
            END DO
          END DO
        END DO
      END DO
C...for aerosol number concentration
C...    convert to # per mole of air by volume

      STRT = 1
      FINI = NNAE
      DO VAR = STRT, FINI
        SPC = NAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
!             CGRID( COL, ROW, LAY, SPC ) = CGRID( COL, ROW, LAY, SPC )
!    &                                    * MWAIR / ( 1000.0
!    &                                    * DENS( COL, ROW, LAY ) )
              FAC = CNV3 / DENS( COL, ROW, LAY )
              CCMIN = MAX(CGRID( COL, ROW, LAY, SPC ), 1.0E-30/FAC)
              CGRID( COL, ROW, LAY, SPC ) = FAC * CCMIN
            END DO
          END DO
        END DO
      END DO

C...for aerosol surface area
C...    convert to m2 per mole of air by volume

      STRT = 1
      FINI = NSAE
      DO VAR = STRT, FINI
        SPC = SAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              FAC = CNV2 / DENS( COL, ROW, LAY )
              CCMIN = MAX(CGRID( COL, ROW, LAY, SPC ), 1.0E-30/FAC)
              CGRID( COL, ROW, LAY, SPC ) = FAC * CCMIN
            END DO
          END DO
        END DO
      END DO

C...determine no. of moles of non-reactive gas per moles of air by volume

      STRT = NR_STRT
      FINI = NR_STRT + N_NR_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = MAX(CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E-6, 1.0E-30)
            END DO
          END DO
        END DO
      END DO

C...determine no. of moles of tracer gas per moles of air by volume

      STRT = TR_STRT
      FINI = TR_STRT + N_TR_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = MAX(CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E-6, 1.0E-30)
            END DO
          END DO
        END DO
      END DO

C...compute cloud effects for the resolved clouds and resolved rainwater

      CALL RESCLD ( CGRID, JDATE, JTIME, TSTEP,
     &              N_SPC_WDEP, WDEP_MAP, TOT_DEP )

C...compute cloud effects for convective (subgrid) clouds

      CALL CONVCLD_ACM ( CGRID, JDATE, JTIME, TSTEP,
     &                   N_SPC_WDEP, WDEP_MAP, CONV_DEP )

C...now check to see if it's time to write the deposition file

      CALL NEXTIME ( MDATE, MTIME, TSTEP( 2 ) ) ! set mdate:mtime to next tstep
      NDATE = JDATE
      NTIME = 10000 * ( JTIME / 10000 )     !  clear MMSS fields
      CALL NEXTIME ( NDATE, NTIME, 10000 )   ! set Ndate:Ntime to next hour

      IF ( ( MDATE .EQ. NDATE ) .AND. ( MTIME .GE. NTIME ) ) THEN

        DO VAR = 1, N_SPC_WDEP + 1
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              TOT_DEP( COL, ROW, VAR ) = TOT_DEP ( COL, ROW, VAR )
     &                                 + CONV_DEP( COL, ROW, VAR )
            END DO
          END DO
        END DO

C...write data to the normal wet deposition file

        IF ( .NOT. WRITE3( CTM_WET_DEP_1, ALLVAR3, NDATE,
     &                     NTIME, TOT_DEP ) ) THEN
          XMSG = 'Could not write ' // CTM_WET_DEP_1 // ' file'
          CALL M3EXIT ( PNAME, NDATE, NTIME, XMSG, XSTAT1 )
        END IF

        WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CTM_WET_DEP_1,
     &         'for date and time', NDATE, NTIME

C...write data to the diagnostic file if requested by the user

        IF ( CLD_DIAG ) THEN

          IF ( .NOT. WRITE3( CTM_WET_DEP_2, ALLVAR3, NDATE,
     &                       NTIME, CONV_DEP ) ) THEN
          XMSG = 'Could not write ' // CTM_WET_DEP_2 // ' file'
            CALL M3EXIT ( PNAME, NDATE, NTIME, XMSG, XSTAT1 )
          END IF

          WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &           'Timestep written to', CTM_WET_DEP_2,
     &           'for date and time', NDATE, NTIME

        END IF   ! CLD_DIAG

C...reinitialize deposition array

        TOT_DEP  = 0.0
        CONV_DEP = 0.0

      END IF   ! time to write

C...convert units on cgrid back to original units
C...  convert to ppmV gas

      STRT = GC_STRT
      FINI = GC_STRT + N_GC_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E+6
            END DO
          END DO
        END DO
      END DO

C...convert to ug/m3 of aerosol mass

      STRT = 1
      FINI = NQAE
      DO VAR = STRT, FINI
        SPC = QAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              FAC = CNV1I * DENS( COL, ROW, LAY )
              CGRID( COL, ROW, LAY, SPC ) = FAC * CGRID( COL, ROW, LAY, SPC )
     &                                    * AE_MOLWT( SPC - AE_STRT + 1 )
            END DO
          END DO
        END DO
      END DO

C...convert to #/m3 of aerosol number

      STRT = 1
      FINI = NNAE
      DO VAR = STRT, FINI
        SPC = NAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
!             CGRID( COL, ROW, LAY, SPC ) = CGRID( COL, ROW, LAY, SPC )
!    &                                    * DENS( COL, ROW, LAY )
!    &                                    * 1000.0 / MWAIR
              FAC = CNV3I * DENS( COL, ROW, LAY )
              CGRID( COL, ROW, LAY, SPC ) = FAC * CGRID( COL, ROW, LAY, SPC )
            END DO
          END DO
        END DO
      END DO

C...convert to m2/m3 of aerosol surface area

      STRT = 1
      FINI = NSAE
      DO VAR = STRT, FINI
        SPC = SAE( VAR )
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              FAC = CNV2I * DENS( COL, ROW, LAY )
              CGRID( COL, ROW, LAY, SPC ) = FAC * CGRID( COL, ROW, LAY, SPC )
            END DO
          END DO
        END DO
      END DO

C...convert to ppmV non-reactive gas

      STRT = NR_STRT
      FINI = NR_STRT + N_NR_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E+6
            END DO
          END DO
        END DO
      END DO

C...convert to ppmV tracer gas

      STRT = TR_STRT
      FINI = TR_STRT + N_TR_SPC - 1
      DO SPC = STRT, FINI
        DO LAY = 1, NLAYS
          DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
              CGRID( COL, ROW, LAY, SPC ) = CGRID( COL, ROW, LAY, SPC )
     &                                    * 1.0E+6
            END DO
          END DO
        END DO
      END DO

      RETURN

      END
