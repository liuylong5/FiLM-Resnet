
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/wr_conc.F,v 1.1 2010/11/17 19:22:55 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     SUBROUTINE WR_CONC ( SGRID, JDATE, JTIME, TSTEP )
      SUBROUTINE WR_CONC ( JDATE, JTIME, TSTEP )

C Revision History:
C   10/13/99 David Wong at LM
C      -- Called from driver, where CGRID is a pointer (subset) of PCGRID.
C         Necessary, to keep from referencing parts of PCGRID that don't
C         belong to CGRID.
C    1/31/2000 Jeff Young
C      -- f90 memory mgmt
C   Jeff - Dec 00 - move CGRID_MAP into f90 module
C   Jeff - Feb 01 - assumed shape arrays
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; DBUFF for WRITE3
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C   30 May 05 J.Young: optional save derived vert. vel. to conc file
C   29 Aug 06 J.Young: in conjunction with state CGRID file
C   30 Jan 10 David Wong: modified the code to write one species at a time
C                         in case ALLVAR3 cannot deal with large CONC size
C   Aug 2010 by Peter Percell at UH-IMAQS:
C   - tried to improve robustness by opening read-only except on processor 0
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE STD_CONC              ! standard CONC
      USE WVEL_DEFN             ! derived vertical velocity component

      IMPLICIT NONE

C Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
#include      SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! I/O definitions and declarations

!     REAL      :: SGRID( :,:,:,: )
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! output timestep (HHMMSS)

      INTEGER, EXTERNAL :: TRIMLEN  ! string length, excl. trailing blanks

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'WR_CONC'
      CHARACTER( 96 ) :: XMSG = ' '

      REAL, ALLOCATABLE, SAVE :: JACF( :,:,: ) ! full-layer Jacobian
      REAL, ALLOCATABLE, SAVE :: MSFX2( :,: )  ! map scale factor ** 2
      REAL, ALLOCATABLE :: DBUFF ( :,:,: )     ! input for W_VEL
      REAL, ALLOCATABLE :: WBUFF ( :,:,: )     ! output for W_VEL
      INTEGER      ALLOCSTAT

      INTEGER, SAVE :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

      INTEGER      GXOFF, GYOFF              ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3
      INTEGER       :: STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2

      INTEGER      C, R, K, L, V   ! loop induction variables

      CHARACTER( 16 ), ALLOCATABLE, SAVE :: VNAME( : )
      INTEGER, SAVE :: NVARS

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         IF ( W_VEL ) THEN

            CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                      STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

            ALLOCATE ( JACF( NCOLS,NROWS,NLAYS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating JACF'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            CALL SUBHFILE ( GRID_CRO_2D, GXOFF, GYOFF,
     &                      STRTCOLGC2, ENDCOLGC2, STRTROWGC2, ENDROWGC2 )

            ALLOCATE ( MSFX2( NCOLS,NROWS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating MSFX2'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            IF ( .NOT. INTERPX( GRID_CRO_2D, 'MSFX2', PNAME,
     &                          STRTCOLGC2,ENDCOLGC2, STRTROWGC2,ENDROWGC2, 1,1,
     &                          JDATE, JTIME, MSFX2 ) ) THEN
               XMSG = 'Could not interpolate MSFX2 from ' // GRID_CRO_2D
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         END IF

C open conc file for update

         IF ( MYPE .EQ. 0 ) THEN
            IF ( .NOT. OPEN3( CTM_CONC_1, FSRDWR3, PNAME ) ) THEN
               XMSG = 'Could not open ' // CTM_CONC_1 // ' file for update'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         ELSE
            IF ( .NOT. OPEN3( CTM_CONC_1, FSREAD3, PNAME ) ) THEN
               XMSG = 'Could not open ' // CTM_CONC_1 // ' file as OLD:READ-ONLY'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         IF ( .NOT. DESC3( CTM_CONC_1 ) ) THEN
            XMSG = 'Could not get file description from ' // CTM_CONC_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( W_VEL ) THEN
            NVARS = NVARS3D - 1
         ELSE
            NVARS = NVARS3D
         END IF

         ALLOCATE ( VNAME( NVARS ), STAT=ALLOCSTAT )

         VNAME = VNAME3D( 1:NVARS )

      END IF   ! firstime

!     IF ( .NOT. WRITE3( CTM_CONC_1, ALLVAR3, JDATE, JTIME, SGRID ) ) THEN
!        XMSG = 'Could not write ' // CTM_CONC_1
!        CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!     END IF

      DO V = 1, NVARS
         IF ( .NOT. WRITE3( CTM_CONC_1, VNAME( V ), JDATE, JTIME,
     &                      SGRID( :,:,:,V ) ) ) THEN
            XMSG = 'Could not write ' // CTM_CONC_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
      END DO

      IF ( W_VEL ) THEN

         IF ( .NOT. INTERPX( MET_CRO_3D, 'JACOBF', PNAME,
     &                       STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                       1,NLAYS, JDATE, JTIME, JACF ) ) THEN
            XMSG = 'Could not interpolate JACOBF from MET_CRO_3D - '
     &           // 'Using JACOBM <- KLUDGE!'
            CALL M3WARN( PNAME, JDATE, JTIME, XMSG )

            IF ( .NOT. ALLOCATED( DBUFF ) ) THEN
               ALLOCATE ( DBUFF( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
               IF ( ALLOCSTAT .NE. 0 ) THEN
                  XMSG = 'Failure allocating DBUFF'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF

            IF ( .NOT. INTERPX( MET_CRO_3D, 'JACOBM', PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                          1,NLAYS, JDATE, JTIME, DBUFF ) ) THEN
               XMSG = 'Could not interpolate JACOBM from MET_CRO_3D'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            DO L = 1, NLAYS - 1
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     JACF( C,R,L ) = 0.5 * ( DBUFF( C,R,L ) + DBUFF( C,R,L+1 ) )
                  END DO
               END DO
            END DO
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  JACF( C,R,NLAYS ) = 0.4 * JACF( C,R,NLAYS-1 )
     &                              +      DBUFF( C,R,NLAYS )
               END DO
            END DO

         END IF   ! INTERPX for JACOBF

C convert from contravariant vertical velocity component to true wind

         IF ( .NOT. ALLOCATED( WBUFF ) ) THEN
            ALLOCATE ( WBUFF( MY_NCOLS,MY_NROWS,C_NLAYS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating WBUFF'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         K = 0
         DO L = CONC_BLEV, CONC_ELEV
            K = K + 1
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  WBUFF( C,R,K ) = JACF( C,R,L ) * MSFX2( C,R ) * WVEL( C,R,L )
               END DO
            END DO
         END DO

         IF ( .NOT. WRITE3( CTM_CONC_1, 'W_VEL', JDATE, JTIME, WBUFF ) ) THEN
            XMSG = 'Could not write W_VEL to ' // CTM_CONC_1
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

      END IF   ! W_VEL

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &      'Timestep written to', CTM_CONC_1,
     &      'for date and time', JDATE, JTIME

      RETURN 
      END
