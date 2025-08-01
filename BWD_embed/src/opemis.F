
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/opemis.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE OPEMIS ( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM, EMLAYS )

C  7 Mar 02 - J.Young: add units string variations
C 29 Oct 05 - J.Young: dyn. layers
C 19 Feb 08 - David Wong: add EM_TRAC = .TRUE. when EMIS_TRAC_1 exists

      USE VGRD_DEFN           ! vertical layer specifications

      IMPLICIT NONE

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
      INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
!     INCLUDE SUBST_EMLYRS_ID ! emissions layers parameter
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

      INTEGER      JDATE      ! current model date, coded YYYYDDD
      INTEGER      JTIME      ! current model time, coded HHMMSS
      INTEGER      NEMIS      ! no. of gas chem emissions species for vdiff
      LOGICAL      EM_TRAC    ! are there tracer emissions?
      REAL         CONVEM     ! conversion for emissions rates
      INTEGER      EMLAYS     ! emissions layers from file

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: INDEX1, TRIMLEN, SETUP_LOGDEV

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'OPEMIS'
      CHARACTER( 96 ) :: XMSG
      CHARACTER( 16 ) :: UNITSCK

      LOGICAL ::   WRFLG = .FALSE.
      INTEGER      LOGDEV
      INTEGER      V, N, S        ! induction variables

C-----------------------------------------------------------------------
 
!     LOGDEV = INIT3()
      LOGDEV = SETUP_LOGDEV()

C Open the tracer emissions file

      IF ( N_TR_EMIS .GT. 0 ) THEN

         IF ( .NOT. OPEN3( EMIS_TRAC_1, FSREAD3, PNAME ) ) THEN

            XMSG = 'Could not open '// EMIS_TRAC_1 // ' file'
            CALL M3MESG( XMSG )
            EM_TRAC = .FALSE.

         ELSE

            EM_TRAC = .TRUE.
            IF ( .NOT. DESC3( EMIS_TRAC_1 ) ) THEN
               XMSG = 'Could not get '// EMIS_TRAC_1 // ' file description'
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

            EMLAYS = NLAYS3D
 
            V = INDEX1( TR_EMIS( 1 ), NVARS3D, VNAME3D )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D( V )
            ELSE
               XMSG = 'Emissions species '
     &              // TR_EMIS( 1 )( 1:TRIMLEN( TR_EMIS( 1 ) ) )
     &              // ' not found on ' // EMIS_1
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

            DO N = 2, N_TR_EMIS
               V = INDEX1( TR_EMIS( N ), NVARS3D, VNAME3D )
               IF ( V .NE. 0 ) THEN
                  IF ( UNITS3D( V ) .NE. UNITSCK ) THEN
                     XMSG = 'Units not uniform on ' // EMIS_TRAC_1
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
                  END IF
               ELSE
                  XMSG = 'Emissions species '
     &                 // TR_EMIS( N )( 1:TRIMLEN( TR_EMIS( 1 ) ) )
     &                 // ' not found on ' // EMIS_TRAC_1
                  CALL M3MESG( XMSG )
               END IF
            END DO

         END IF

      END IF   ! tracer emissions

!     IF ( NEMIS + N_AE_EMIS + N_NR_EMIS + N_TR_EMIS .GT. 0 ) THEN
      IF ( NEMIS + N_NR_EMIS .GT. 0 ) THEN

C Open the emissions file (for gas chem, aerosols and non-reactive species)

         IF ( .NOT. OPEN3( EMIS_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open '// EMIS_1 // ' file'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( .NOT. DESC3( EMIS_1 ) ) THEN
            XMSG = 'Could not get '// EMIS_1 // ' file description'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         END IF

         IF ( EMLAYS .EQ. 0 ) THEN
            EMLAYS = NLAYS3D
         ELSE
!slzx            IF ( EMLAYS .NE. NLAYS3D ) THEN
!slzx               WRITE( XMSG,'(A, I3, A, I3)' ) 'Tracer No. of emis layers ', EMLAYS,
!slzx     &         ' Different than Gas No. of emis layers: ', NLAYS3D 
!slzx               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!slzx            END IF
         END IF

         DO N = 1, N_GC_EMIS
            V = INDEX1( GC_EMIS( N ), NVARS3D, VNAME3D )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D( V )
               S = N + 1
               GO TO 101
            ELSE
               XMSG = 'Emissions species '
     &              // GC_EMIS( 1 )( 1:TRIMLEN( GC_EMIS( 1 ) ) )
     &              // ' not found on ' // EMIS_1
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF
         END DO
         XMSG = ' No emissions species ' // ' found on ' // EMIS_1
         CALL M3WARN( PNAME, JDATE, JTIME, XMSG )

101      CONTINUE

         DO N = S, N_GC_EMIS
            V = INDEX1( GC_EMIS( N ), NVARS3D, VNAME3D )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D( V )
               CALL UPCASE ( UNITSCK )
               IF ( UNITSCK .NE. 'MOLES/S'   .AND.
     &              UNITSCK .NE. 'MOLE/S'    .AND.
     &              UNITSCK .NE. 'MOL/S'     .AND.
     &              UNITSCK .NE. 'MOLES/SEC' .AND.
     &              UNITSCK .NE. 'MOLE/SEC'  .AND.
     &              UNITSCK .NE. 'MOL/SEC' ) THEN
                  XMSG = 'GC units incorrect on ' // EMIS_1
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            ELSE
               XMSG = 'Emissions species '
     &              // GC_EMIS( N )( 1:TRIMLEN( GC_EMIS( N ) ) )
     &              // ' not found on ' // EMIS_1
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF
         END DO

         DO N = 1, N_NR_EMIS
            V = INDEX1( NR_EMIS( N ), NVARS3D, VNAME3D )
            IF ( V .NE. 0 ) THEN
               UNITSCK = UNITS3D( V )
               CALL UPCASE ( UNITSCK )
               IF ( UNITSCK .NE. 'MOLES/S'   .AND.
     &              UNITSCK .NE. 'MOLE/S'    .AND.
     &              UNITSCK .NE. 'MOL/S'     .AND.
     &              UNITSCK .NE. 'MOLES/SEC' .AND.
     &              UNITSCK .NE. 'MOLE/SEC'  .AND.
     &              UNITSCK .NE. 'MOL/SEC' ) THEN
                  XMSG = 'NR units incorrect on ' // EMIS_1
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF
            ELSE
               XMSG = 'Emissions species '
     &              // NR_EMIS( N )( 1:TRIMLEN( NR_EMIS( N ) ) )
     &              // ' not found on ' // EMIS_1
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
            END IF
         END DO

C Assume units = mol/s across gas and non-reactive species classes

         CONVEM = 1.0E-03  ! assuming gram-moles

      END IF  ! (gas chem or non-react emissions)

      IF ( NEMIS .GT. 0 ) THEN
         WRFLG = .TRUE.
         WRITE( LOGDEV,1001 )
1001     FORMAT( / 10X, 'Gas Chemistry Emissions Processing in',
     &              1X, 'Vertical diffusion ...' )
      END IF

!     IF ( N_AE_EMIS .GT. 0 ) THEN
!        WRFLG = .TRUE.
!        WRITE( LOGDEV,1003 )
1003     FORMAT( / 10X, 'Aerosol Emissions Processing in',
     &              1X, 'Vertical diffusion ...' )
!     END IF

      IF ( N_NR_EMIS .GT. 0 ) THEN
         WRFLG = .TRUE.
         WRITE( LOGDEV,1005 )
1005     FORMAT( / 10X, 'Non-reactives Emissions Processing in',
     &              1X, 'Vertical diffusion ...' )
      END IF

      IF ( N_TR_EMIS .GT. 0 .AND. EM_TRAC ) THEN
         WRFLG = .TRUE.
         WRITE( LOGDEV,1007 )
1007     FORMAT( / 10X, 'Tracer Emissions Processing in',
     &              1X, 'Vertical diffusion ...' )
      END IF

!slzdbg      EMLAYS = 1

      RETURN
      END
