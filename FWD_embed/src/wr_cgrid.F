
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_fwd/wr_cgrid.F,v 1.1 2010/11/17 19:22:55 sergey Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE WR_CGRID ( CGRID, JDATE, JTIME, TSTEP )

C Save the instantaneous state of CGRID in an I/O-API "circular buffer"
C file to use for restart/continuation for subsequent simulation. This
C file will replace the CONC file for such use, allowing users to reduce
C the number of variables and/or layers saved to the CONC file.

C Revision History:
C      May 06 J.Young: initial
C      Feb 08 J.Young: fix VDESC3D bug
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets

      USE SUBST_MODULES         ! stenex
!     USE SUBST_UTIL_MODULE     ! stenex

      IMPLICIT NONE

C Include Files:

      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structure
#include      SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP        ! output timestep (HHMMSS)

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'WR_CGRID'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER, EXTERNAL :: TRIMLEN  ! string length, excl. trailing blanks

      INTEGER K, MXK, SPC, VAR
      INTEGER ALLOCSTAT

      INTEGER TSTEP_RF, NTHIK_RF, NCOLS_RF, NROWS_RF, GDTYP_RF
      REAL( 8 ) :: P_ALP_RF, P_BET_RF, P_GAM_RF
      REAL( 8 ) :: XCENT_RF, YCENT_RF
      REAL( 8 ) :: XORIG_RF, YORIG_RF
      REAL( 8 ) :: XCELL_RF, YCELL_RF
      INTEGER VGTYP_RF
      REAL VGTOP_RF

      INTEGER, SAVE :: LOGDEV       ! FORTRAN unit number for log file
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      LOGICAL OK

#ifdef parallel
      INTERFACE
         LOGICAL FUNCTION PTRWRITE3( S_CGRID, ALLVAR3, JDATE, JTIME, CGRID )
            IMPLICIT NONE
            CHARACTER( 16 ), INTENT( IN ) :: S_CGRID, ALLVAR3
            INTEGER, INTENT( IN )         :: JDATE, JTIME
            REAL, POINTER                 :: CGRID( :,:,:,: )
         END FUNCTION PTRWRITE3
      END INTERFACE
#endif

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

C Try to open existing file for update

         CALL SUBST_BARRIER

         OK = OPEN3( S_CGRID, FSRDWR3, PNAME )
         CALL SUBST_GLOBAL_LOGICAL( OK, 'AND' )
         IF ( .NOT. OK ) THEN

            XMSG = 'Could not open ' // S_CGRID( 1:TRIMLEN( S_CGRID ) )
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

C Set tstep negative to create "circular buffer" type file

!              TSTEP3D = -TSTEP
               SDATE3D = JDATE
               STIME3D = JTIME

C Get CGRID nvars

               CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

               NVARS3D = NSPCSD
               NLAYS3D = NLAYS

C Set file header attributes that differ from CONC and open the file

               FDESC3D = ' '
               FDESC3D( 1 ) = 'Computational grid instantaneous concentrations'
               FDESC3D( 2 ) = '- for scenario continuation.'

               WRITE( LOGDEV,* ) ' '
               WRITE( LOGDEV,* ) '       State CGRID File Header Description:'
               DO K = 1, 2
                  WRITE( LOGDEV,* ) '    => ',
     &            FDESC3D( K )( 1:TRIMLEN( FDESC3D( K ) ) )
               END DO

               VAR = 0

               DO SPC = 1, N_GC_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = GC_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
               END DO

C pick up transported RHOJ

               VAR = VAR + 1
               VTYPE3D( VAR ) = M3REAL
               VNAME3D( VAR ) = 'RHOJ'
               UNITS3D( VAR ) = 'm*Kg/m**3'
               VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )

               DO SPC = 1, N_AE_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = AE_SPC( SPC )
                  IF ( VNAME3D( VAR )(1:3) .EQ. 'NUM' ) THEN
                     UNITS3D( VAR ) = 'number/m**3'
                  ELSE IF ( VNAME3D( VAR )(1:3) .EQ. 'SRF' ) THEN
                     UNITS3D( VAR ) = 'm**2/m**3'
                  ELSE
                     UNITS3D( VAR ) = 'micrograms/m**3'
                  END IF
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
               END DO

               DO SPC = 1, N_NR_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = NR_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
               END DO

               DO SPC = 1, N_TR_SPC
                  VAR = VAR + 1
                  VTYPE3D( VAR ) = M3REAL
                  VNAME3D( VAR ) = TR_SPC( SPC )
                  UNITS3D( VAR ) = 'ppmV'
                  VDESC3D( VAR ) = 'Variable ' // VNAME3D( VAR )
               END DO

               WRITE( LOGDEV,* ) ' '
               WRITE( LOGDEV,* ) '       State CGRID File Variable List:'
               DO SPC = 1, VAR
                  WRITE( LOGDEV,'( 5X, "=> VNAME3D(", I3, " ): ", A )' )
     &                   SPC, VNAME3D( SPC )
               END DO

               IF ( .NOT. OPEN3( S_CGRID, FSNEW3, PNAME ) ) THEN
                  XMSG = 'Could not open '
     &                 // S_CGRID( 1:TRIMLEN( S_CGRID ) ) // ' file'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            END IF   ! MYPE = 0

         ELSE

C File exists. Check header data with CONC file as reference. Currently only
C proc 0 has CTM_CONC_1 open

            IF ( MYPE .EQ. 0 ) THEN

               IF ( .NOT. DESC3( CTM_CONC_1 ) ) THEN
                  XMSG = 'Could not get '
     &                 // CTM_CONC_1( 1:TRIMLEN( CTM_CONC_1 ) )
     &                 // ' file description'
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

               IF ( .NOT. DESC3( S_CGRID ) ) THEN
                  XMSG = 'Could not get '
     &                 // S_CGRID( 1:TRIMLEN( S_CGRID ) )
     &                 // ' file description'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

               IF ( TSTEP_RF .NE. ABS( TSTEP3D ) .OR.
     &              NTHIK_RF .NE. NTHIK3D .OR.
     &              NCOLS_RF .NE. NCOLS3D .OR.
     &              NROWS_RF .NE. NROWS3D .OR.
     &              GDTYP_RF .NE. GDTYP3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( P_ALP_RF .NE. P_ALP3D .OR.
     &              P_BET_RF .NE. P_BET3D .OR.
     &              P_GAM_RF .NE. P_GAM3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( XCENT_RF .NE. XCENT3D .OR.
     &              YCENT_RF .NE. YCENT3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( XORIG_RF .NE. XORIG3D .OR.
     &              YORIG_RF .NE. YORIG3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( XCELL_RF .NE. XCELL3D .OR.
     &              YCELL_RF .NE. YCELL3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( VGTYP_RF .NE. VGTYP3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
               IF ( VGTOP_RF .NE. VGTOP3D ) THEN
                    XMSG = 'Header inconsistent on existing S_CGRID'
                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF

            END IF   ! MYPE = 0

         END IF   ! .NOT. OPEN S_CGRID

         CALL SUBST_BARRIER

      END IF   ! FIRSTIME

#ifdef parallel
      IF ( .NOT. PTRWRITE3( S_CGRID, ALLVAR3, JDATE, JTIME, CGRID ) ) THEN
         XMSG = 'Could not write S_CGRID'
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
#else
      IF ( .NOT. WRITE3( S_CGRID, ALLVAR3, JDATE, JTIME, CGRID ) ) THEN
         XMSG = 'Could not write S_CGRID'
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF
#endif

      WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &      'Timestep written to', S_CGRID,
     &      'for date and time', JDATE, JTIME

      RETURN 
      END
