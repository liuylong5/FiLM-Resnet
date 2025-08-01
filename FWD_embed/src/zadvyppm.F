
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vadv/vyamo_cadj_fwd/zadvyppm.F,v 1.4 2013/05/14 21:06:52 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ZADV ( CGRID, JDATE, JTIME, TSTEP )

C-----------------------------------------------------------------------
C Function:
C   Advection in the vertical, x3-direction:
C   The process time step is set equal to TSTEP

C Preconditions:
C   Dates and times represented YYYYDDD:HHMMSS.
C   No "skipped" dates and times. Process time step divides TSTEP exactly
C   CGRID in transport units: SQRT{DET[metric tensor]}*concentration (Mass/Vol)

C Subroutines and functions called:
C   TIME2SEC

C Revision history:
C   02/19/93 by M. Talat Odman  at NCSC
C   05/17/93 by Carlie J. Coats at NCSC:  now uses INTERP3()
C   06/14/94 by Dongming Hwang at NCSC:
C              include statement and subroutine name template
C   10/15/95 by M. Talat Odman at NCSC: generalized coordinates

C   Sep 97 Jeff
C   Aug 98 Jeff better Courant condition tstep limit

C    David Wong, Sep. 1998
C      -- parallelized the code

C    15 Dec 00 J.Young: move CGRID_MAP into f90 module
C                       GLOBAL_RSUM -> Dave Wong's f90 stenex GLOBAL_SUM
C                       GLOBAL_ISUM -> Dave Wong's f90 stenex GLOBAL_SUM

C    28 Jul 01 J.Young: allocatable arrays ...
C                       Since F90 does not preserve dummy argument array
C                       indices, the 3rd dimension of WHAT has been changed
C                       from 0:NLAYS to 1:NLAYS+1 for the sake of vcontvel

C    03 Sep 01 David Wong
C      -- inserted F90 DEALLOCATE statement for NX3
C
C   1/03 - JP modified for Yamo mass conservation
C          Vertical velocity is diagnosed from mass continuity
C          vertical advection is upstream (no call to adv scheme)

C    31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module
C    27 Apr 07 J.Young: Talat's First-order upstream (donor cell) algorithm
C    30 Apr 09 J.Pleim, J.Young: Replace donor cell with ppm, adjust velocity
C                                accordingly
C    21 Aug 09 J.Young: Don't bypass VPPMY if ITER = 0
C    18 Nov 09 J.Young: Combine VPPMY and VPPM functionality
C    Aug 2010 by Peter Percell at UH-IMAQS: modified for adjoint/4dvar
C-----------------------------------------------------------------------

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE CGRID_SPCS            ! CGRID species number and offsets
      USE WVEL_DEFN             ! derived vertical velocity component
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines
      IMPLICIT NONE

C Includes:

      INCLUDE SUBST_GC_SPC      ! gas chemistry species table
      INCLUDE SUBST_AE_SPC      ! aerosol species table
      INCLUDE SUBST_NR_SPC      ! non-reactive species table
      INCLUDE SUBST_TR_SPC      ! tracer species table
      INCLUDE SUBST_GC_ADV      ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV      ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV      ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV      ! tracer advection species and map table
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters

C Arguments:

      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER     JDATE         ! current model date, coded YYYYDDD
      INTEGER     JTIME         ! current model time, coded HHMMSS
      INTEGER     TSTEP( 2 )    ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C External Functions not declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: TIME2SEC

C Parameters:

      INTEGER, PARAMETER :: MAXSS = 100 ! error exit limit

C Advected species dimension

      INTEGER, PARAMETER :: N_SPC_ADV = N_GC_ADV
     &     + N_AE_ADV
     &     + N_NR_ADV
     &     + N_TR_ADV
     &     + 1                  ! for advecting air

C File Variables:

      REAL :: RHOJM( NCOLS,NROWS,NLAYS ) ! RhoJ (Kg/m**2) from Met file
      REAL :: WY   ( NLAYS,NCOLS,NROWS ) ! Diagnosed vert vel ala yamo

C Local variables:

      CHARACTER( 16 ) :: PNAME = 'ZADVYPPM'
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

C for INTERPX
      INTEGER      GXOFF, GYOFF ! global origin offset from file
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      INTEGER       ETIME, EDATE
      REAL          VEL ( NLAYS+1, MAXSS+1 ) ! Velocities in a N-S column
      REAL          FLX ( NLAYS+1 ) ! upstream donor cell computed conc. flux
      REAL, ALLOCATABLE, SAVE :: DS ( : ) ! dx3 (dimensionless in sigma coord.)
      REAL          DTSEC       ! sync time step in seconds
      REAL          DELT        ! adjusted time step
      REAL          FLUX        ! intermediate flux

      INTEGER, SAVE :: ADV_MAP( N_SPC_ADV ) ! global adv map to CGRID

      INTEGER       COL, ROW, LVL, SPC, VAR ! loop counters
      INTEGER       A2C
      INTEGER       SS, NSS

      CHARACTER( 96 ) :: XMSG = ' '
      INTEGER, SAVE :: ASPC     ! pointer in CGRID to transported RHOJ
      REAL          RJT( NLAYS ) ! local adjusted RHOJ
      REAL          RJM( NLAYS ) ! local RHOJM at tstep + 1
      REAL          DSTM        ! subexpression
      REAL          CC          ! local Courant No.
      REAL          DTNEW       ! sub timestep
      REAL          DSDT        ! DS/DT
      REAL          NSTEP       ! integral DT divisor
      INTEGER       ALLOCSTAT

      INTEGER, SAVE :: LOGDEV

      REAL :: DT(MAXSS+1)

      INTERFACE
         SUBROUTINE VPPM(NI, DT, DS, VEL, CON, FLX)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NI
            REAL, INTENT(IN)    :: DT
            REAL, INTENT(IN)    :: DS(NI)
            REAL, INTENT(INOUT) :: VEL(NI+1)
            REAL, INTENT(INOUT) :: CON(NI)
            REAL, INTENT(IN), OPTIONAL :: FLX(NI+1)
         END SUBROUTINE VPPM
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3 ()

         ALLOCATE ( DS( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DS'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C Get default file header attibutes from MET_CRO_3D (assumes file already open)

C Get dx3 from COORD include file
         WRITE( LOGDEV,* ) ' '
         WRITE( LOGDEV,* ) '    layer    S (X3FACE_GD) Delta S'
         DO LVL = 1, NLAYS
            DS ( LVL ) = ABS ( X3FACE_GD( LVL ) - X3FACE_GD( LVL-1 ) )
            WRITE( LOGDEV,'(5X, I3, 3F14.7)' ) LVL, X3FACE_GD( LVL ),
     &           DS( LVL )
         END DO
         WRITE( LOGDEV,* ) ' '

C Get CGRID offsets
         CALL CGRID_MAP( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C Pointer to transported RHOJ
         ASPC = GC_STRT - 1 + N_GC_SPCD

C Create global map to CGRID
         SPC = 0
         DO VAR = 1, N_GC_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = GC_STRT - 1 + GC_ADV_MAP( VAR )
         END DO
         DO VAR = 1, N_AE_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = AE_STRT - 1 + AE_ADV_MAP( VAR )
         END DO
         DO VAR = 1, N_NR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = NR_STRT - 1 + NR_ADV_MAP( VAR )
         END DO
         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            ADV_MAP( SPC ) = TR_STRT - 1 + TR_ADV_MAP( VAR )
         END DO

         ADV_MAP( N_SPC_ADV ) = ASPC

C get file local domain offsets
         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &        STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C Create new adjoint file and open for update
C         CALL CHK_FILE_OPEN ( ADJ_VA_RHOJ_CHK, JDATE, JTIME, .TRUE., TSTEP )

      END IF                    ! FIRSTIME

C Checkpoint advected rhoj
C      CALL CHK_FILE_WRITE ( ADJ_VA_RHOJ_CHK, JDATE, JTIME, CGRID(:, :, :, ASPC:ASPC) )

C Time-stepped gridded computation for Z-direction advection.
      DTSEC  = FLOAT( TIME2SEC( TSTEP( 2 ) ) ) ! process time step (seconds)

C vertical velocities are at face centers, positive upward.
C No boundary conditions are needed because VEL(1) = VEL(NLAYS+1) = 0

C Get rho*J at end of sync step
      EDATE = JDATE
      ETIME = JTIME
      CALL NEXTIME( EDATE, ETIME, TSTEP( 2 ) )
      IF ( .NOT. INTERPX( MET_CRO_3D, 'DENSA_J', PNAME,
     &     STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &     EDATE, ETIME, RHOJM ) ) THEN
         XMSG = 'Could not read DENSA_J from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, EDATE, ETIME, XMSG, XSTAT1 )
      END IF

      DO ROW = 1, MY_NROWS
         DO COL = 1, MY_NCOLS

            DO LVL = 1, NLAYS
               RJM( LVL ) = RHOJM( COL,ROW,LVL )
               RJT( LVL ) = CGRID( COL,ROW,LVL,ASPC )
            END DO

            DELT = DTSEC
            VEL( 1,: ) = 0.0    ! impermeable boundary condition at the surface
            FLX( 1 ) = 0.0

            DO SS = 1, MAXSS + 1
               FLUX = 0.0
               DO LVL = 1, NLAYS
                  FLUX = FLUX - ( RJM( LVL ) - RJT( LVL ) ) * DS( LVL ) / DELT
                  FLX( LVL+1 ) = FLUX
               END DO

               DO LVL = 2, NLAYS
                  IF ( FLX( LVL ) .GE. 0.0 ) THEN
                     VEL( LVL,SS ) = FLX( LVL ) / RJT( LVL-1 )
                  ELSE
                     VEL( LVL,SS ) = FLX( LVL ) / RJT( LVL )
                  END If
               END DO

               VEL( NLAYS+1,SS ) = FLX( NLAYS+1 ) / RJT( NLAYS )

C Find Maximum Courant Number

               CC = 0.0

               DO LVL = 2, NLAYS
                  IF ( VEL( LVL,SS ) .GT. 0.0 ) THEN
                     CC = MAX ( CC, ( VEL( LVL,SS ) * DELT / DS( LVL-1 ) ) )
                  ELSE
                     CC = MAX ( CC, ( -VEL( LVL,SS ) * DELT / DS( LVL ) ) )
                  END IF
               END DO

               LVL = NLAYS+1
               IF ( VEL( LVL,SS ) .GT. 0.0 ) THEN
                  CC = MAX ( CC, ( VEL( LVL,SS ) * DELT / DS( LVL-1 ) ) )
               ELSE
                  CC = MAX ( CC, ( -VEL( LVL,SS ) * DELT / DS( LVL-1 ) ) )
               END IF

               IF ( CC .GT. 1.0 ) THEN ! courant number is larger than unity

C Calculate a sub-time step that satisfies the Courant stability limit.
C Perform vertical advection with the computed velocity and sub-time step.
C Then calculate the difference between the original and sub-time steps.
C The difference is the new sub-time step. Recompute vertical velocities
C that would bring the air density field back to being uniform. Note that
C if Courant number with the new velocity and sub-time step is larger than
C unity again, then the last sub-time step would be split into further
C sub-steps.

                  NSTEP = FLOAT( INT( CC ) + 1 )
                  DTNEW = DELT / NSTEP
                  DELT = DELT - DTNEW

                  DT(SS) = DTNEW
                  CALL VPPM ( NLAYS, DT(SS), DS, VEL(:,SS), RJT, FLX )

                  IF ( SS .GT. MAXSS ) THEN
                     WRITE( LOGDEV,2005 ) COL, ROW, CC, DELT, SS, JTIME
 2005                FORMAT( 'zadv col  row     CC', 8X, 'dt    iter   jtime'
     &                    / 'zzzz', 2I4, 1PE12.3, 0PF10.5, 1X, I4, I10.6
     &                    / 10X, 'MetRhoj', 3X, 'TrRhoj', 5X, 'Diff',
     &                    4X, 'adv_rhoj', 3X, 'vel(l)', 6X, 'vel(l+1)' )

                     DO LVL = 1, NLAYS
                        WRITE( LOGDEV,2009 ) LVL, RJM( LVL ), RJT( LVL ),
     &                       RJM( LVL ) - RJT( LVL ),
     &                       VEL( LVL,SS ), VEL( LVL+1,SS )
                     END DO
 2009                FORMAT( 'zzz2', I3, 3F10.2, 2(1PE12.3) )

                     WRITE( XMSG,2013 ) JTIME, TSTEP( 2 ), MAXSS
 2013                FORMAT( 'vert adv soln failed at', I7.6,  ' with adv step:',
     &                    I7.6, ' HHMMSS', 2X, 'Max Iterations =', I3 )
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                  END IF

               ELSE

                  DT(SS) = DELT
                  CALL VPPM ( NLAYS, DT(SS), DS, VEL(:,SS), RJT, FLX )
                  EXIT

               END IF

            END DO
            NSS = SS

            DO LVL = 1, NLAYS
               CGRID( COL,ROW,LVL,ASPC ) = RJT( LVL )
            END DO

C Update concentrations of trace species
            DO SPC = 1, N_SPC_ADV - 1
               A2C = ADV_MAP( SPC )
               DO SS = 1, NSS
                  CALL VPPM(NLAYS, DT(SS), DS, VEL(:,SS), CGRID(COL,ROW,:,A2C))
               END DO
            END DO

            DO LVL = 1, NLAYS
               WY( LVL,COL,ROW ) = VEL( LVL+1,NSS )
            END DO

         END DO                 ! COL
      END DO                    ! ROW

      IF ( W_VEL ) THEN
         DO LVL = 1, NLAYS
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS
                  WVEL( COL,ROW,LVL ) = WY( LVL,COL,ROW )
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE ZADV
