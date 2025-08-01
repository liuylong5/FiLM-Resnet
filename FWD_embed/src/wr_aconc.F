
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/wr_aconc.F,v 1.1 2010/11/17 19:22:55 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WR_ACONC ( AGRID, JDATE, JTIME, TSTEP )

C Revision History:
C   Jeff - July 01
C   Note: If previous A_CONC exists, check that user hasn't changed what
C         species/layers to save (or domain).
C   30 Mar 01 J.Young: dyn alloc - Use HGRD_DEFN
C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE AVG_CONC              ! integral average CONC

      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex

      IMPLICIT NONE

C Include Files:

      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
#include      SUBST_IODECL      ! I/O definitions and declarations
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_FILES_ID    ! file name parameters
!     INCLUDE SUBST_COORD_ID    ! coord. and domain definitions (req IOPARMS)

!     CHARACTER( 16 ) :: A_CONC_1 = 'A_CONC_1'

      REAL      :: AGRID( :,:,:,: )
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! output timestep (HHMMSS)

C Local variables:

      INTEGER      MDATE        ! modified model date, coded YYYYDDD
      INTEGER      MTIME        ! modified model time, coded HHMMSS

      CHARACTER( 16 ) :: PNAME = 'WR_ACONC'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER, EXTERNAL :: TRIMLEN  ! string length, excl. trailing blanks

      INTEGER, SAVE :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL OK

      INTEGER      L, K, KD, VAR, SPC ! loop counters

      INTEGER, SAVE :: A_NLYS

      INTEGER TSTEP_RF, NTHIK_RF, NCOLS_RF, NROWS_RF, GDTYP_RF
      REAL( 8 ) :: P_ALP_RF, P_BET_RF, P_GAM_RF
      REAL( 8 ) :: XCENT_RF, YCENT_RF
      REAL( 8 ) :: XORIG_RF, YORIG_RF
      REAL( 8 ) :: XCELL_RF, YCELL_RF
      INTEGER VGTYP_RF
      REAL VGTOP_RF
C-----------------------------------------------------------------------

C Change output date/time to starting date/time - e.g. timestamp 1995196:090000
C represents data computed from time 1995196:090000 to 1995196:100000

      MDATE = JDATE
      MTIME = JTIME
      CALL NEXTIME ( MDATE, MTIME, -TSTEP )

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

C Try to open existing file for update

         CALL SUBST_BARRIER

         A_NLYS = ACONC_ELEV - ACONC_BLEV + 1

         OK = OPEN3( A_CONC_1, FSRDWR3, PNAME )
         CALL SUBST_GLOBAL_LOGICAL( OK, 'AND' )
         IF ( .NOT. OK ) THEN

            XMSG = 'Could not open ' // A_CONC_1( 1:TRIMLEN( A_CONC_1 ) )
     &           // ' file for update - try to open new'
            CALL M3MESG( XMSG )

            IF ( MYPE .EQ. 0 ) THEN

C Get default file header attibutes from CONC file (assumes file already open)

               IF ( .NOT. DESC3( CTM_CONC_1 ) ) THEN
                  XMSG = 'Could not get '
     &                 // CTM_CONC_1( 1:TRIMLEN( CTM_CONC_1 ) )
     &                 // ' file description'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF

C Set file header attributes that differ from CONC and open the file

               SDATE3D = MDATE
               STIME3D = MTIME
               TSTEP3D = TSTEP
               NVARS3D = N_ASPCS
               NLAYS3D = A_NLYS

               L = 0
               DO K = ACONC_BLEV, ACONC_ELEV
                  L = L + 1
                  VGLVS3D( L ) = VGLVS_GD( K )
                  END DO
               VGLVS3D( A_NLYS + 1 ) = VGLVS_GD( NLAYS + 1 )
!              GDNAM3D = GDNAME_GD
               GDNAM3D = GRID_NAME  ! from HGRD_DEFN

               FDESC3D( 1 ) = 'Concentration file output '
               FDESC3D( 2 ) = 'Averaged over the synchronization time steps '
               FDESC3D( 3 ) = 'Timestamp represents beginning computed date/time '
               FDESC3D( 4 ) = 'Layer mapping (CGRID to AGRID):'
               KD = 4
               VAR = ACONC_BLEV
               L = 0
               DO K = KD + 1, MIN ( A_NLYS + KD, MXDESC3 )
                  L = L + 1
                  WRITE( FDESC3D( K ),'( "Layer", I3, " to", I3, " " )' )
     &            VAR + L - 1, L
                  END DO
               IF ( ( KD + 1 + L ) .LT. MXDESC3 ) THEN
                  DO K = KD + 1 + L, MXDESC3
                     FDESC3D( K ) = ' '
                     END DO
                  END IF

               WRITE( LOGDEV,* ) ' '
               WRITE( LOGDEV,* ) '      Avg Conc File Header Description:'
               DO K = 1, KD + L
                  WRITE( LOGDEV,* ) '    => ',
     &            FDESC3D( K )( 1:TRIMLEN( FDESC3D( K ) ) )
                  END DO

               VAR = 0

               DO SPC = 1, N_A_GC_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = A_GC_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
                  END DO

               DO SPC = 1, N_A_AE_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = A_AE_SPC( SPC )
                  IF ( VNAME3D( VAR )(1:3) .EQ. 'NUM' ) THEN
                     UNITS3D( VAR ) = 'number/m**3'
                     ELSE IF ( VNAME3D( VAR )(1:3) .EQ. 'SRF' ) THEN
                     UNITS3D( VAR ) = 'm**2/m**3'
                     ELSE
                     UNITS3D( VAR ) = 'micrograms/m**3'
                     END IF
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
                  END DO

               DO SPC = 1, N_A_NR_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = A_NR_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
                  END DO

               DO SPC = 1, N_A_TR_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = A_TR_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
                  END DO

            DO SPC = 1, VAR
               WRITE( LOGDEV,'( 5X, "=> VNAME3D(", I3, " ): ", A )' )
     &                SPC, VNAME3D( SPC )
               END DO

               IF ( .NOT. OPEN3( A_CONC_1, FSNEW3, PNAME ) ) THEN
                  XMSG = 'Could not open '
     &                 // A_CONC_1( 1:TRIMLEN( A_CONC_1 ) ) // ' file'
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                  END IF

               END IF   ! MYPE = 0

            ELSE

C File exists. Make sure it matches requested output.

            IF ( .NOT. DESC3( CTM_CONC_1 ) ) THEN
               XMSG = 'Could not get '
     &              // CTM_CONC_1( 1:TRIMLEN( CTM_CONC_1 ) )
     &              // ' file description'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            TSTEP_RF = TSTEP3D
            NTHIK_RF = NTHIK3D
            NCOLS_RF = NCOLS3D
            NROWS_RF = NROWS3D
            GDTYP_RF = GDTYP3D
            P_ALP_RF = P_ALP3D
            P_BET_RF = P_BET3D
            P_GAM_RF = P_GAM3D
            XCENT_RF = XCENT3D
            YCENT_RF = YCENT3D
            XORIG_RF = XORIG3D
            YORIG_RF = YORIG3D
            XCELL_RF = XCELL3D
            YCELL_RF = YCELL3D
            VGTYP_RF = VGTYP3D
            VGTOP_RF = VGTOP3D

            IF ( .NOT. DESC3( A_CONC_1 ) ) THEN
               XMSG = 'Could not get '
     &              // A_CONC_1( 1:TRIMLEN( A_CONC_1 ) )
     &              // ' file description'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            IF ( N_ASPCS .NE. NVARS3D ) THEN
               WRITE( XMSG, '( A, 2I6 )' )
     &         'Number of variables don''t match file: ', N_ASPCS, NVARS3D
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF

            DO SPC = 1, N_ASPCS
               DO VAR = 1, NVARS3D
                  IF ( AVG_SPCS( SPC ) .EQ. VNAME3D( VAR ) ) GO TO 101
                  END DO
               XMSG = 'Could not find ' // AVG_SPCS( SPC )
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
101            CONTINUE
               END DO

            IF ( A_NLYS .NE. NLAYS3D ) THEN
               WRITE( XMSG, '( A, 2I6 )' )
     &         'Number of layers don''t match file: ', A_NLYS, NLAYS3D
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF

C Check other header data with CONC file as reference

            IF ( TSTEP_RF .NE. TSTEP3D .OR.
     &           NTHIK_RF .NE. NTHIK3D .OR.
     &           NCOLS_RF .NE. NCOLS3D .OR.
     &           NROWS_RF .NE. NROWS3D .OR.
     &           GDTYP_RF .NE. GDTYP3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( P_ALP_RF .NE. P_ALP3D .OR.
     &           P_BET_RF .NE. P_BET3D .OR.
     &           P_GAM_RF .NE. P_GAM3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( XCENT_RF .NE. XCENT3D .OR.
     &           YCENT_RF .NE. YCENT3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( XORIG_RF .NE. XORIG3D .OR.
     &           YORIG_RF .NE. YORIG3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( XCELL_RF .NE. XCELL3D .OR.
     &           YCELL_RF .NE. YCELL3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( VGTYP_RF .NE. VGTYP3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF
            IF ( VGTOP_RF .NE. VGTOP3D ) THEN
                 XMSG = 'Header inconsistent on existing A_CONC_1'
                 CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                 END IF

            END IF   ! .NOT. OPEN A_CONC_1

         CALL SUBST_BARRIER

         END IF   ! FIRSTIME
    
      VAR = 0
    
      DO SPC = 1, N_A_GC_SPC
         VAR = VAR + 1
         IF ( .NOT. WRITE3( A_CONC_1, A_GC_SPC( SPC ),
     &      MDATE, MTIME, AGRID( :,:,:,VAR ) ) ) THEN
            XMSG = 'Could not write '
     &           // A_GC_SPC( SPC )( 1:TRIMLEN( A_GC_SPC( SPC ) ) )
     &           // ' to ' // A_CONC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
!     write( logdev,* ) ' <>var, A_GC_SPC: ', VAR, A_GC_SPC( SPC )
         END DO
 
      DO SPC = 1, N_A_AE_SPC
         VAR = VAR + 1
         IF ( .NOT. WRITE3( A_CONC_1, A_AE_SPC( SPC ),
     &      MDATE, MTIME, AGRID( :,:,:,VAR ) ) ) THEN
            XMSG = 'Could not write '
     &           // A_AE_SPC( SPC )( 1:TRIMLEN( A_AE_SPC( SPC ) ) )
     &           // ' to ' // A_CONC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
!     write( logdev,* ) ' <>var, A_AE_SPC: ', VAR, A_AE_SPC( SPC )
         END DO
 
      DO SPC = 1, N_A_NR_SPC
         VAR = VAR + 1
         IF ( .NOT. WRITE3( A_CONC_1, A_NR_SPC( SPC ),
     &      MDATE, MTIME, AGRID( :,:,:,VAR ) ) ) THEN
            XMSG = 'Could not write '
     &           // A_NR_SPC( SPC )( 1:TRIMLEN( A_NR_SPC( SPC ) ) )
     &           // ' to ' // A_CONC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
!     write( logdev,* ) ' <>var, A_NR_SPC: ', VAR, A_NR_SPC( SPC )
         END DO
 
      DO SPC = 1, N_A_TR_SPC
         VAR = VAR + 1
         IF ( .NOT. WRITE3( A_CONC_1, A_TR_SPC( SPC ),
     &      MDATE, MTIME, AGRID( :,:,:,VAR ) ) ) THEN
            XMSG = 'Could not write '
     &           // A_TR_SPC( SPC )( 1:TRIMLEN( A_TR_SPC( SPC ) ) )
     &           // ' to ' // A_CONC_1
            CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
         END DO

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &      'Timestep written to', A_CONC_1,
     &      'for date and time', MDATE, MTIME

      RETURN 
      END
