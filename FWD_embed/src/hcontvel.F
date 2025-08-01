
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/hcontvel.F,v 1.1 2010/11/17 19:22:55 sergey Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)hcontvel.F        1.1 /project/mod3/CMAQ/src/util/util/SCCS/s.hcontvel.F 03 Jun 1997 12:08:52

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE HCONTVEL ( JDATE, JTIME, TSTEP, LVL, UORV, WIND )
      
C-----------------------------------------------------------------------
C Function:
C    This subroutine reads physical velocities in the x1 or x2 directions
C    and returns the contravariant velocities.
 
C Preconditions:
C    This routine can be used only for conformal map coordinates 
C    in the horizontal.
C    Dates and times should be represented YYYYDDD:HHMMSS.
 
C Subroutines and functions called:
C    INTERPX, INTERPB, M3EXIT, TIME2SEC, SEC2TIME, NEXTIME
      
C Revision history:
C   January 30, 1996 by Clint L. Ingram at NCSC: created for
C   RADM-coordinates

C   22 Apr 97 Jeff:
C    7 Aug 97 Jeff: for NTHIK = 1
C    4 Feb 98 Jeff: deal with end-of-scenario
C   20 Sep 98 David Wong: parallelized the code
C                         -- adjust the data declaration for DENSJ
C                         -- remove indirect index reference, and re-adapt to
C                            a general case
C                         -- invoke stencil exchange library
C   21 Nov 00 J.Young: PE_COMM3 -> Dave Wong's f90 stenex COMM
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN; replace INTERP3 with INTERPX
C    6 Apr 01 J.Young: Eliminate NTHIN confusion (assumes NTHIK = 1)
C   12 Apr 01 J.Young: Use PINTERPB for boundary data
C   23 Jun 03 J.Young: for layer dependent advection tstep
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE SUBST_MODULES         ! stenex
!     USE SUBST_COMM_MODULE     ! stenex

      IMPLICIT NONE

C Includes:

!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! horizontal dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structuer
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
      INCLUDE SUBST_PE_COMM     ! PE communication displacement and direction
 
C Parameters

      CHARACTER( 16 ), PARAMETER :: UFLAG = 'X1VEL'
      CHARACTER( 16 ), PARAMETER :: VFLAG = 'X2VEL'
      
C Arguments:
      
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! time step (HHMMSS)
      INTEGER      LVL          ! layer
      CHARACTER( 16 ) :: UORV   ! flag for velocity component
!     REAL         WIND( NCOLS+1,NROWS+1 ) ! CX xi-velocity 
      REAL      :: WIND( :,: )  ! CX xi-velocity 
      
C Parameters:

C file variables:
      
!     REAL      DENSJ_BUF( NCOLS,NROWS,NLAYS )     ! Jacobian * air density
      REAL, ALLOCATABLE, SAVE :: DENSJ_BUF( :,:,: ) ! Jacobian * air density
!     REAL, ALLOCATABLE, SAVE :: DENSJ_BUF( :,: )   ! Jacobian * air density
!     REAL      DENSJ_BND( NBNDY,NLAYS )           ! bndy Jacobian * air density
      REAL, ALLOCATABLE, SAVE :: DENSJ_BND( :,: )  ! bndy Jacobian * air density
!     REAL    DENSJ( 0:NCOLS+1,0:NROWS+1,NLAYS )
      REAL, ALLOCATABLE, SAVE :: DENSJ( :,: )      ! Jacobian * air density
 
C External Functions (not already declared by IODECL3.EXT):
 
      INTEGER, EXTERNAL :: SEC2TIME, TIME2SEC, TRIMLEN
#ifdef parallel
      LOGICAL, EXTERNAL :: INTERPB
#endif
      
C local variables:
      
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV
       
      INTEGER   ROW               ! Row index
      INTEGER   COL               ! Column index
      INTEGER   MDATE             ! mid-advection date
      INTEGER   MTIME             ! mid-advection time
      INTEGER   STEP              ! advection time step in seconds
      INTEGER, SAVE :: LDATE( 3 ) ! last date for data on file
      INTEGER, SAVE :: LTIME( 3 ) ! last time for data on file
      LOGICAL   REVERT            ! recover last time step if true
      REAL      DJ                ! temporary Jacobian * air density
      INTEGER   ALLOCSTAT
 
      CHARACTER( 16 ) :: VNAME
      CHARACTER( 16 ) :: PNAME = 'HCONTVEL'
      CHARACTER( 16 ) :: AMSG
      CHARACTER( 96 ) :: XMSG = ' '
 
      CHARACTER( 8 ), SAVE :: COMMSTR

      INTEGER COUNT     ! Counter for constructing density array.

      INTEGER    GXOFF, GYOFF ! global origin offset from file
      LOGICAL, SAVE :: WINDOW = .FALSE. ! posit same file and global
                                        ! processing domain
      INTEGER, SAVE :: NCOLSDENS, NROWSDENS       ! local for DENSJ_BUF
C for INTERPX
      INTEGER, SAVE :: STRTCOL,   ENDCOL,   STRTROW,   ENDROW
      INTEGER       :: STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC
      INTEGER, SAVE :: STRTCOLMD, ENDCOLMD, STRTROWMD, ENDROWMD

C-----------------------------------------------------------------------
 

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         CALL LSTEPF( MET_CRO_3D, LDATE( 1 ), LTIME( 1 ) )
         CALL LSTEPF( MET_BDY_3D, LDATE( 2 ), LTIME( 2 ) )
         CALL LSTEPF( MET_DOT_3D, LDATE( 3 ), LTIME( 3 ) )
 
         LDATE( 1 ) = MIN( LDATE( 1 ), LDATE( 2 ), LDATE( 3 ) )
         LTIME( 1 ) = SEC2TIME( MIN(
     &                         TIME2SEC( LTIME( 1 ) ),
     &                         TIME2SEC( LTIME( 2 ) ),
     &                         TIME2SEC( LTIME( 3 ) )
     &                         ) )

         WRITE( COMMSTR,'(4I2)' )  1, 1-NTHIK, 2, 1-NTHIK  ! ' 1 0 2 0'

         CALL SUBHFILE ( MET_DOT_3D, GXOFF, GYOFF,
     &                   STRTCOLMD, ENDCOLMD, STRTROWMD, ENDROWMD )
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC, ENDCOLMC, STRTROWMC, ENDROWMC )

         NCOLSDENS = ENDCOLMC - STRTCOLMC + 1
         NROWSDENS = ENDROWMC - STRTROWMC + 1

         IF ( NCOLSDENS .NE. MY_NCOLS .OR.
     &        NROWSDENS .NE. MY_NROWS ) THEN
            WRITE( XMSG,'( A, 4I8 )' ) 'Local Columns or Rows incorrect',
     &      NCOLSDENS, MY_NCOLS, NROWSDENS, MY_NROWS
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

C currently not implemented: case where only one origin component matches file's
         IF ( GXOFF .NE. 0 .AND. GYOFF .NE. 0 ) THEN
            WINDOW = .TRUE. ! windowing from file
            STRTCOL = STRTCOLMC - 1
            ENDCOL  = ENDCOLMC  + 1
            STRTROW = STRTROWMC - 1
            ENDROW  = ENDROWMC  + 1
            ELSE
            STRTCOL = STRTCOLMC
            ENDCOL  = ENDCOLMC
            STRTROW = STRTROWMC
            ENDROW  = ENDROWMC
            END IF

         ALLOCATE ( DENSJ( 0:MY_NCOLS+1,0:MY_NROWS+1 ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DENSJ'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         IF ( .NOT. WINDOW ) THEN

            ALLOCATE ( DENSJ_BUF( NCOLSDENS,NROWSDENS,NLAYS ),STAT = ALLOCSTAT )
!           ALLOCATE ( DENSJ_BUF( NCOLSDENS,NROWSDENS ),      STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DENSJ_BUF'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            ALLOCATE ( DENSJ_BND( NBNDY,NLAYS ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating DENSJ_BND'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            END IF

         END IF                    ! if firstime
 
      MDATE  = JDATE
      MTIME  = JTIME
      STEP   = TIME2SEC( TSTEP )
      CALL NEXTIME( MDATE, MTIME, SEC2TIME( STEP / 2 ) )

      IF ( MDATE .LT. LDATE( 1 ) ) THEN
         REVERT = .FALSE.
         ELSE IF ( MDATE .EQ. LDATE( 1 ) ) THEN
            IF ( MTIME .LE. LTIME( 1 ) ) THEN
               REVERT = .FALSE.
               ELSE
               REVERT = .TRUE.
               END IF
         ELSE   ! MDATE .GT. LDATE
         REVERT = .TRUE.
         END IF
 
      IF ( REVERT ) THEN
         XMSG = 'Current scenario interpolation step not available in all of '
     &        // MET_CRO_3D(1:TRIMLEN( MET_CRO_3D ) ) // ', '
     &        // MET_BDY_3D(1:TRIMLEN( MET_BDY_3D ) ) // ' and '
     &        // MET_DOT_3D(1:TRIMLEN( MET_DOT_3D ) )
         CALL M3MESG( XMSG )
!        CALL NEXTIME( MDATE, MTIME, -SEC2TIME( STEP / 2 ) )
         WRITE( AMSG,'( 2I8 )' ) LDATE( 1 ), LTIME( 1 )
         XMSG = 'Using data for last file step: ' // AMSG
         CALL M3MESG( XMSG )
         MDATE = LDATE( 1 )
         MTIME = LTIME( 1 )
         END IF
 
      IF ( LVL .GT. NLAYS ) THEN
         LVL = NLAYS
         XMSG = 'Adjusting layer down to NLAYS'
         CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
         END IF

C Interpolate Jacobian X Air Density
 
      IF ( WINDOW ) THEN

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPX ( MET_CRO_3D, VNAME, PNAME,
     &                        STRTCOL,ENDCOL, STRTROW,ENDROW, LVL,LVL,
     &                        MDATE, MTIME, DENSJ ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         ELSE ! need to extend data from bndy file

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPX ( MET_CRO_3D, VNAME, PNAME,
!    &                        STRTCOL,ENDCOL, STRTROW,ENDROW, LVL,LVL,
!    &                        STRTCOL,ENDCOL, STRTROW,ENDROW, LVL,LVL+1,
     &                        STRTCOL,ENDCOL, STRTROW,ENDROW, 1,NLAYS,
     &                        MDATE, MTIME, DENSJ_BUF ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_CRO_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         VNAME = 'DENSA_J'
         IF ( .NOT. INTERPB ( MET_BDY_3D, VNAME, PNAME,
     &                        MDATE, MTIME, NBNDY*NLAYS,
     &                        DENSJ_BND ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_BDY_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

C Load DENSJ array

         DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS
               DENSJ( COL,ROW ) = DENSJ_BUF( COL,ROW,LVL )
!              DENSJ( COL,ROW ) = DENSJ_BUF( COL,ROW )
               END DO
            END DO

C Fill in DENSJ array for boundaries

         COUNT = 0
         DO ROW = 0, 0
            DO COL = 1, MY_NCOLS+1
               COUNT = COUNT + 1
               DENSJ( COL,ROW ) = DENSJ_BND( COUNT,LVL )  ! South
               END DO
            END DO
         DO ROW = 1, MY_NROWS+1
            DO COL = MY_NCOLS+1, MY_NCOLS+1
               COUNT = COUNT + 1
               DENSJ( COL,ROW ) = DENSJ_BND( COUNT,LVL )  ! East
               END DO
            END DO
         DO ROW = MY_NROWS+1, MY_NROWS+1
            DO COL = 0, MY_NCOLS
               COUNT = COUNT + 1
               DENSJ( COL,ROW ) = DENSJ_BND( COUNT,LVL )  ! North
               END DO
            END DO
         DO ROW = 0, MY_NROWS
            DO COL = 0, 0
               COUNT = COUNT + 1
               DENSJ( COL,ROW ) = DENSJ_BND( COUNT,LVL )  ! West
               END DO
            END DO

!        DEALLOCATE ( DENSJ_BUF )
!        DEALLOCATE ( DENSJ_BND )

         END IF   ! WINDOW

C Interpolate Contravariant Velocity components (already at flux points)
C X Jacobian X Air Density

      IF ( UORV .EQ. UFLAG ) THEN

         VNAME = 'UHAT_JD'   ! x1 component of CX-vel * Jacobian * air density
         IF ( .NOT. INTERPX ( MET_DOT_3D, VNAME, PNAME,
     &                        STRTCOLMD,ENDCOLMD, STRTROWMD,ENDROWMD, LVL,LVL,
     &                        MDATE, MTIME, WIND ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_DOT_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         ELSE IF ( UORV .EQ. VFLAG ) THEN

         VNAME = 'VHAT_JD'   ! x2 component of CX-vel * Jacobian * air density
         IF ( .NOT. INTERPX ( MET_DOT_3D, VNAME, PNAME,
     &                        STRTCOLMD,ENDCOLMD, STRTROWMD,ENDROWMD, LVL,LVL,
     &                        MDATE, MTIME, WIND ) ) THEN
            XMSG = 'Could not read ' // VNAME // ' from ' // MET_DOT_3D
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         ELSE
         XMSG = 'UORV Flag, ' // VNAME // ' not set correctly'
         CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
         END IF


C Obtain flux point values of Jacobian * air density and retrieve
C contravariant velocities 

      IF ( UORV .EQ. UFLAG ) THEN

C store actual north, east, south, and west displacement
C define communication pattern

         CALL SUBST_COMM ( DENSJ, DSPL_N0_E0_S0_W1, DRCN_W, COMMSTR )
         DO ROW = 1, MY_NROWS
            DO COL = 1, MY_NCOLS+1
               DJ = 0.5 * ( DENSJ( COL,ROW) + DENSJ( COL-1,ROW ) )
               WIND( COL,ROW ) = WIND( COL,ROW ) / DJ
               END DO
            END DO

         ELSE

C store actual north, east, south, and west displacement
C define communication pattern

         CALL SUBST_COMM ( DENSJ, DSPL_N0_E0_S1_W0, DRCN_S, COMMSTR )
         DO ROW = 1, MY_NROWS+1
            DO COL = 1, MY_NCOLS
               DJ = 0.5 * ( DENSJ( COL,ROW ) + DENSJ( COL,ROW-1 ) )
               WIND( COL,ROW ) = WIND( COL,ROW ) / DJ
               END DO
            END DO
         END IF

!     DEALLOCATE ( DENSJ )

      RETURN
      END
