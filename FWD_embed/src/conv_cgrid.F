
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/conv_cgrid.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE CONV_CGRID ( CGRID, JDATE, JTIME, CNGRD )

C-----------------------------------------------------------------------
C Function:
C   Convert decoupled aerosol species to molar units (ppm and m**2/mol)
C   and reorder dimensions

C Revision History:
C   Written by: J.Young 21 Aug 03
C   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module
C-----------------------------------------------------------------------

      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE CGRID_SPCS          ! CGRID species number and offsets

      IMPLICIT NONE

!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_GC_SPC    ! gas chemistry species table
      INCLUDE SUBST_AE_SPC    ! aerosol species table
      INCLUDE SUBST_NR_SPC    ! non-reactive species table
      INCLUDE SUBST_TR_SPC    ! tracer species table
      INCLUDE SUBST_CONST     ! constants
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include      SUBST_IODECL    ! I/O definitions and declarations

C Arguments:

      REAL, POINTER :: CGRID( :,:,:,: )              !  concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      REAL :: CNGRD( :,:,:,: )  ! cgrid replacement

C Parameters:

      REAL, PARAMETER :: GPKG = 1.0E+03        ! g/Kg
      REAL, PARAMETER :: MAOGPKG = MWAIR / GPKG
      REAL, PARAMETER :: GPKGOMA = 1.0 / MAOGPKG
      REAL, PARAMETER :: MAOAVO1000 = 1.0E+03 * MWAIR / AVO
      REAL, PARAMETER :: AVOOMA_001 = 1.0 / MAOAVO1000

C Local Variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV
      CHARACTER( 16 ), SAVE :: PNAME = 'CONV_CGRID'
      CHARACTER( 96 ) :: XMSG = ' '

      REAL      DENS( NCOLS,NROWS,NLAYS )  ! air density

      INTEGER   NSPCS, OFF
      INTEGER   C, R, L, S, V            ! loop induction variables
      INTEGER, SAVE :: NQAE              ! number of micro-grams / m**3 species
      INTEGER, SAVE :: QAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NNAE              ! number of  # / m**3 species
      INTEGER, SAVE :: NAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      INTEGER, SAVE :: NSAE              ! number of  m**2 / m**3 species
      INTEGER, SAVE :: SAE( N_AE_SPCD )  ! CGRID pointer to "" "" ""
      REAL, SAVE    :: MOLWT( N_AE_SPCD ) ! only for "QAE" species

      REAL    CONV, FAC            ! temp var

      INTEGER     GXOFF, GYOFF          ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C Get number of species, and starting indices for CGRID array.

         CALL CGRID_MAP ( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

C create aerosol species pointers to distinguish micro-grams / m**3,
C # / m**3 (number density), and m**2 / m**3 (surface area) species

         NQAE = 0       ! no. of micro-grams / m**3 species
         NNAE = 0       ! no. of  # / m**3 species
         NSAE = 0       ! no. of  m**2 / m**3 species
         OFF = AE_STRT - 1
         DO S = 1, N_AE_SPC
            IF ( AE_SPC( S )( 1:3 ) .EQ. 'NUM' ) THEN
               NNAE = NNAE + 1
               NAE( NNAE ) = OFF + S
               ELSE IF ( AE_SPC( S )( 1:3 ) .EQ. 'SRF' ) THEN
               NSAE = NSAE + 1
               SAE( NSAE ) = OFF + S
               ELSE
               NQAE = NQAE + 1
               QAE( NQAE ) = OFF + S
               MOLWT( NQAE ) = AE_MOLWT( S )
               END IF
            END DO

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         END IF          !  if Firstime

      IF ( .NOT. INTERPX( MET_CRO_3D, 'DENS', PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME, DENS ) ) THEN
         XMSG = 'Could not interpolate DENS from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

C Convert non-molar mixing ratio species and re-order CGRID

C Gas - no conversion

      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  DO V = 1, NSPCS
                     CNGRD( OFF+V,L,C,R ) = CGRID( C,R,L,OFF+V )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C micro-grams/m**3 aerosol -> mol/mol air <- no
C micro-grams/m**3 aerosol -> ppmv
C (Don't divide by MGPG, then multiply by 1.0E+6: 1/MGPG = 10**-6 cancels out
C ppm = 10**6)

      NSPCS = NQAE
      IF ( NSPCS .GT. 0 ) THEN
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  FAC = MAOGPKG / DENS( C,R,L )
                  DO V = 1, NSPCS
                     CONV = FAC / MOLWT( V )
                     CNGRD( QAE( V ),L,C,R ) = CONV * CGRID( C,R,L,QAE( V ) )
                     END DO
                  END DO
               END DO
            END DO
         END IF

!     NSPCS = NQAE
!     IF ( NSPCS .GT. 0 ) THEN
!        DO V = 1, NSPCS
!           CONV = MAOGPKG / MOLWT( V )
!           DO L = 1, NLAYS
!              DO R = 1, MY_NROWS
!                 DO C = 1, MY_NCOLS
!                    CGRID( C,R,L,QAE( V ) ) = CGRID( C,R,L,QAE( V ) ) * CONV
!    &                                       / DENS( C,R,L )
!                    END DO
!                 END DO
!              END DO
!           END DO
!        END IF


C number/m**3 aerosol -> ppmv
C (Don't divide by MGPG, etc. See note above)

      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  CONV = MAOAVO1000 / DENS( C,R,L )
                  DO V = 1, NSPCS
                     CNGRD( NAE( V ),L,C,R ) = CONV * CGRID( C,R,L,NAE( V ) )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C m**2/m**3 aerosol -> m**2/mol air

      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  CONV = MAOGPKG / DENS( C,R,L )
                  DO V = 1, NSPCS
                     CNGRD( SAE( V ),L,C,R ) = CONV * CGRID( C,R,L,SAE( V ) )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C Non-reactives - no conversion

      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  DO V = 1, NSPCS
                     CNGRD( OFF+V,L,C,R ) = CGRID( C,R,L,OFF+V )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C Tracers - no conversion

      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
               DO L = 1, NLAYS
                  DO V = 1, NSPCS
                     CNGRD( OFF+V,L,C,R ) = CGRID( C,R,L,OFF+V )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ENTRY REV_CGRID ( CNGRD, JDATE, JTIME, CGRID )

C Revert non-molar mixing ratio species and re-order CGRID

      IF ( .NOT. INTERPX( MET_CRO_3D, 'DENS', PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    JDATE, JTIME, DENS ) ) THEN
         XMSG = 'Could not interpolate DENS from ' // MET_CRO_3D
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF

C Gas - no conversion

      NSPCS = N_GC_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = GC_STRT - 1
         DO V = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+V ) = CNGRD( OFF+V,L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C aerosol ppmv -> micro-grams/m**3
C (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
C ppm = 10**6)

      NSPCS = NQAE
      IF ( NSPCS .GT. 0 ) THEN
         DO V = 1, NSPCS
            FAC = GPKGOMA * MOLWT( V )
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CONV = FAC * DENS( C,R,L )
                     CGRID( C,R,L,QAE( V ) ) = CONV * CNGRD( QAE( V ),L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C aerosol ppmv -> number/m**3
C (Don't multiply by MGPG, etc. See note above)

      NSPCS = NNAE
      IF ( NSPCS .GT. 0 ) THEN
         DO V = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CONV = AVOOMA_001 * DENS( C,R,L )
                     CGRID( C,R,L,NAE( V ) ) = CONV * CNGRD( NAE( V ),L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C m**2/m**3 aerosol -> m**2/mol air

      NSPCS = NSAE
      IF ( NSPCS .GT. 0 ) THEN
         DO V = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CONV = GPKGOMA * DENS( C,R,L )
                     CGRID( C,R,L,SAE( V ) ) = CONV * CNGRD( SAE( V ),L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C Non-reactives - no conversion

      NSPCS = N_NR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = NR_STRT - 1
         DO V = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+V ) = CNGRD( OFF+V,L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

C Tracers - no conversion

      NSPCS = N_TR_SPC
      IF ( NSPCS .GT. 0 ) THEN
         OFF = TR_STRT - 1
         DO V = 1, NSPCS
            DO L = 1, NLAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     CGRID( C,R,L,OFF+V ) = CNGRD( OFF+V,L,C,R )
                     END DO
                  END DO
               END DO
            END DO
         END IF

      RETURN
      END

