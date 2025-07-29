
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

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ACMCLD_ADJ ( F, C, C_ADJ, SIGMAF, CBELOW, CBELOW_ADJ,
     &                        CLBASE, CLTOP, FRAC, NSP, NLAYS, TCLIFE, 
     &                        DTCLD )
C-----------------------------------------------------------------------
C
C  FUNCTION:  Adjoint of Subroutine to compute convective mixing in the 
C             CBL according to the Asymmetrical Convective Model (ACM).CBL
C             Ref: Pleim snd Chang (1992)
C
C  SUMMARY:
C   ACM is based on the Blackadar non-local convective model which is
C   used in HIRPBL where upward mixing similar to Blackadar but
C   downward mixing is to the next lower level representing more
C   realistic gradual subsidence.
C
C  REVISION  HISTORY:
C      Date   Who             What
C    -------- ---             -----------------------------------------
C     06/2005 J.Pleim         Initial version
C     07/2005 J.Young         Clean up for CMAQ-F
C     01/2013 M.Turner        Initial development of adjoint
C-----------------------------------------------------------------------

!slzdbg      USE FOURDVAR_MOD, ONLY: COL_SAVE, ROW_SAVE
!slzdbg      USE GRID_CONF, ONLY: MYPE

      IMPLICIT NONE

C Arguments

      INTEGER NSP                ! no. of species
      INTEGER NLAYS              ! no. of model layers
      REAL    F( NLAYS )         ! entrainment fraction
      REAL    C( NSP, NLAYS )    ! species concentration
      REAL    SIGMAF( 0:NLAYS )  ! full layer sigma (mono decr)
      REAL    CBELOW( NSP )      ! spec conc in layer below cld base
      INTEGER CLBASE, CLTOP
      REAL    FRAC               ! grid cell fractional cloud cover
      REAL    TCLIFE             ! cloud lifetime (s)
      REAL    DTCLD              ! cloud integration time step

C Parameters

      REAL, PARAMETER :: HALF = 0.5
      REAL, PARAMETER :: CRANKP = 0.5

C Local variables

      INTEGER NLP, K, NL, S          ! index variables
      INTEGER KB

      REAL DTLIM, F1
      REAL TOT1, TOT2
      REAL DTS, DELC, M1UP
      REAL( 8 ) :: AI( NLAYS ), BI( NLAYS ), EI( NLAYS )
      REAL( 8 ) :: DI( NLAYS ), UI( NLAYS )
      REAL( 8 ) :: ALPHA, BETA, GAMA
      REAL VCI( NLAYS,NSP )
      REAL MBARKS( NLAYS ), MDWN( NLAYS )
      REAL DSIGH( NLAYS ), DSIGHI( NLAYS )

! Adjoint variables
      REAL ::  C_ADJ( NSP, NLAYS )
      REAL ::  CBELOW_ADJ( NSP )
      REAL( 8 ) :: UI_ADJ( NLAYS )
      REAL( 8 ) :: DI_ADJ( NLAYS )
      REAL( 8 ) :: BETA_ADJ
      REAL :: VCI_ADJ( NLAYS, NSP )
      REAL :: F1_ADJ
      REAL :: DELC_ADJ

! Checkpointing variables
      REAL( 8 ) :: BI_SAVE( NLAYS, NSP, NLAYS )
      REAL( 8 ) :: AI_SAVE( NLAYS, NSP, NLAYS )
      REAL( 8 ) :: EI_SAVE( NLAYS, NSP, NLAYS )
      REAL( 8 ) :: BI_SAVE1( NLAYS, NSP, NLAYS )
      REAL( 8 ) :: AI_SAVE1( NLAYS, NSP, NLAYS )
      REAL( 8 ) :: ALPHA_SAVE( NLAYS, NSP, (CLTOP - CLBASE) + CLBASE )
      REAL( 8 ) :: GAMA_SAVE( NLAYS, NSP )

! Initialize adjoint variables
      UI_ADJ = 0.0
      DI_ADJ = 0.0
      BETA_ADJ = 0.0
      VCI_ADJ = 0.0
      F1_ADJ = 0.0
      DELC_ADJ = 0.0

C-----------------------------------------------------------------------

      DTLIM = DTCLD
      MDWN ( CLTOP + 1 ) = 0.0
      DSIGH( CLTOP + 1 ) = 1.0
      SIGMAF( 0 ) = 1
      M1UP = 0.0
      KB  = CLBASE - 1
      CLBASE = CLBASE
      DSIGH ( KB ) = SIGMAF( KB ) - 1.0
      DSIGHI( KB ) = 1.0 / DSIGH( KB )

C Compute ACM mixing rate

      DO K = CLTOP, CLBASE, -1
        DSIGH ( K ) = SIGMAF( K ) - SIGMAF( K - 1 )
        DSIGHI( K ) = 1.0 / DSIGH( K )
        MBARKS( K ) = ( 1.0 - F( K ) ) * FRAC / TCLIFE
        MDWN  ( K ) = MBARKS( K )
     &              + MDWN( K + 1 ) * DSIGH( K + 1 ) * DSIGHI( K )
        M1UP  = M1UP + MBARKS( K ) * DSIGH( K )
        DTLIM = MIN( HALF / ( M1UP * DSIGHI( K ) ), DTLIM )
      END DO
      DTLIM = MIN( HALF / ( M1UP * DSIGHI( KB ) ), DTLIM )

      DO S = 1, NSP
        VCI( KB, S ) = CBELOW( S )
        VCI( CLTOP+1,S ) = 9999.0
        DO K = CLBASE, CLTOP
          VCI( K,S ) = C( S,K )
          UI( K )  = 0.0           ! init variable for use below
        END DO
      END DO

      NLP = INT( DTCLD / DTLIM + 1.0 )
      DTS = ( DTCLD / NLP )
      DO NL = 1, NLP      ! loop over sub timestep
        DO S = 1, NSP     ! loop over species
                                                                              
C Compute tendency of CBL concentrations - Semi-Implicit solution

          DO K = CLBASE, CLTOP
            DELC = DTS
     &           * ( MBARKS( K ) * VCI( KB,S )
     &           -   MDWN( K ) * VCI( K,S )
     &           +   DSIGH( K+1 ) * DSIGHI( K ) * MDWN( K+1 ) * VCI( K+1,S ) )
            DI( K ) = VCI( K,S ) + ( 1.0 - CRANKP ) * DELC
            EI( K ) = -CRANKP * MDWN( K ) * DTS * DSIGH( K ) * DSIGHI( K-1 )
            BI( K ) = 1.0 + CRANKP * MDWN( K ) * DTS
            AI( K ) = -CRANKP * MBARKS( K ) * DTS
          END DO

          BI( KB ) = 1.0 + CRANKP * M1UP * DTS * DSIGHI( KB )
          F1 = M1UP * VCI( KB,S )
     &       - MDWN( CLBASE ) * VCI( CLBASE,S ) * DSIGH( CLBASE )
          DI( KB ) = VCI( KB,S ) - ( 1.0 - CRANKP ) * F1 * DSIGHI( KB ) * DTS

C Define arrays A,B,E which make up matrix and D which is RHS

          BETA = DI( KB )
          GAMA = BI( KB )
          ALPHA = 1.0
          DO K = CLBASE, CLTOP
            ALPHA = -ALPHA * EI( K ) / BI( K )

            ! Checkpointing, mdt
            ALPHA_SAVE( NL, S, K ) = ALPHA

            BETA  = ALPHA * DI( K ) + BETA
            GAMA  = ALPHA * AI( K ) + GAMA
          END DO

          ! Checkpointing, mdt
          GAMA_SAVE( NL, S ) = GAMA

          UI( KB )   = BETA / GAMA

          ! Checkpointing, mdt
          BI_SAVE1( NL, S, : ) = BI
          AI_SAVE1( NL, S, : ) = AI

          UI( CLTOP ) = ( DI( CLTOP ) - AI( CLTOP ) * UI( KB ) ) / BI( CLTOP )

          BETA = DI( KB )
          GAMA = BI( KB )
          ALPHA = 1.0

C Back substitution:
          ! Checkpointing, mdt
          BI_SAVE( NL, S, : ) = BI
          AI_SAVE( NL, S, : ) = AI
          EI_SAVE( NL, S, : ) = EI

          DO K = CLTOP - 1, CLBASE, -1
            UI( K ) = ( DI( K ) - AI( K ) * UI( KB ) - EI( K+1 ) * UI( K+1 ) )
     &              / BI( K )
          END DO

C Update concentrations
          DO K = KB, CLTOP
            VCI( K,S ) = UI( K )
          END DO

        END DO                   ! end loop for species
      END DO                 ! end timestep loop

      DO S = 1, NSP
        TOT1 = 0.0
        TOT2 = 0.0

        DO K = 1, CLTOP
          TOT1 = TOT1 + C( S,K ) * ( SIGMAF( K-1 ) - SIGMAF( K ) )
        END DO

        CBELOW( S ) = VCI( KB,S )

        DO K = CLBASE, CLTOP
          C( S,K ) = VCI( K,S )
        END DO

      END DO

!---------------------------------------------------------------------
! Adjoint Code Begins

      DO S = NSP, 1, -1

         DO K = CLTOP, CLBASE, -1

            !------
            ! fwd code:
            ! C( S,K ) = VCI( K,S )
            ! adj code:
            VCI_ADJ( K, S ) = VCI_ADJ( K, S ) + C_ADJ( S, K )
            C_ADJ( S, K ) = 0.0

         END DO

         VCI_ADJ( KB, S ) = VCI_ADJ( K, S ) + CBELOW_ADJ( S )
         CBELOW_ADJ (S) = 0.d0

      END DO

      DO NL = NLP, 1, -1

         DO S = NSP, 1, -1

            DO K = CLTOP, KB, -1

               !------
               ! fwd code:
               ! VCI( K,S ) = UI( K )
               ! adj code:
               UI_ADJ ( K ) = UI_ADJ ( K ) + VCI_ADJ ( K, S )
               VCI_ADJ ( K, S ) = 0.0

            END DO

            ! Checkpointing, mdt
            BI = BI_SAVE( NL, S, : )
            AI = AI_SAVE( NL, S, : )
            EI = EI_SAVE( NL, S, : )

            DO K = CLBASE, CLTOP - 1

               !------
               ! fwd code:
               ! UI( K ) = ( DI( K ) - AI( K ) * UI( KB ) - EI( K+1 ) 
               !           * UI( K+1 ) ) / BI( K ) 
               ! adj code:
               DI_ADJ ( K ) = DI_ADJ ( K ) + UI_ADJ ( K ) / BI( K )
               UI_ADJ ( KB ) = UI_ADJ ( KB ) - AI( K ) * UI_ADJ ( K )
     &                       / BI( K )
               UI_ADJ ( K + 1 ) = UI_ADJ ( K + 1 ) - EI( K + 1 ) * 
     &                            UI_ADJ ( K ) / BI( K )
               UI_ADJ ( K ) = 0.0

            END DO

            !------
            ! fwd code:
            ! BETA = DI( KB )
            ! adj code:
            DI_ADJ ( KB ) = DI_ADJ ( KB ) + BETA_ADJ
            BETA_ADJ = 0.0

            ! Checkpointing, mdt
            BI = BI_SAVE1( NL, S, : )
            AI = AI_SAVE1( NL, S, : )

            !------
            ! fwd code:
            ! UI( CLTOP ) = ( DI( CLTOP ) - AI( CLTOP ) * UI( KB 
            !               ) ) / BI( CLTOP )
            ! adj code:
            DI_ADJ ( CLTOP ) = DI_ADJ ( CLTOP ) + UI_ADJ( CLTOP ) 
     &                       / BI( CLTOP )
            UI_ADJ ( KB ) = UI_ADJ ( KB ) - AI( CLTOP ) * UI_ADJ ( 
     &                      CLTOP ) / BI( CLTOP )
            UI_ADJ ( CLTOP ) = 0.0

            ! Checkpointing, mdt
            GAMA = GAMA_SAVE( NL, S )

            !------
            ! fwd code:
            ! UI( KB ) = BETA / GAMA
            ! adj code:
            BETA_ADJ = BETA_ADJ + UI_ADJ ( KB ) / GAMA 
            UI_ADJ ( KB ) = 0.0

            DO K = CLTOP, CLBASE, -1
 
               ! Checkpointing, mdt
               ALPHA = ALPHA_SAVE( NL, S, K )
 
               !------
               ! fwd code:
               ! BETA = ALPHA * DI( K ) + BETA
               ! adj code:
               DI_ADJ ( K ) = DI_ADJ ( K ) + ALPHA * BETA_ADJ

            END DO

            !------
            ! fwd code:
            ! BETA = DI( KB ) 
            ! adj code:
            DI_ADJ ( KB ) = DI_ADJ ( KB ) + BETA_ADJ
            BETA_ADJ = 0.0

            !------
            ! fwd code:
            ! DI( KB ) = VCI( KB,S ) - ( 1.0 - CRANKP ) * F1 * 
            !            DSIGHI( KB ) * DTS
            ! adj code:
            VCI_ADJ ( KB, S ) = VCI_ADJ ( KB, S ) + DI_ADJ ( KB )
            F1_ADJ = F1_ADJ - ( 1.0 - CRANKP ) * DSIGHI( KB ) * DTS *
     &               DI_ADJ ( KB )
            DI_ADJ ( KB ) = 0.0

            !------
            ! fwd code:
            ! F1 = M1UP * VCI( KB,S ) - MDWN( CLBASE ) * 
            !      VCI( CLBASE,S ) * DSIGH( CLBASE )
            ! adj code:
            VCI_ADJ ( KB, S ) = VCI_ADJ ( KB, S ) + M1UP * F1_ADJ
            VCI_ADJ ( CLBASE, S ) = VCI_ADJ ( CLBASE, S ) - MDWN ( 
     &                              CLBASE ) * DSIGH( CLBASE ) * F1_ADJ
            F1_ADJ = 0.0

            DO K = CLTOP, CLBASE, -1
 
               !------
               ! fwd code:
               ! DI( K ) = VCI( K,S ) + ( 1.0 - CRANKP ) * DELC
               ! adj code:
               VCI_ADJ ( K, S ) = VCI_ADJ ( K, S ) + DI_ADJ ( K )
               DELC_ADJ = DELC_ADJ + ( 1.0 - CRANKP ) * DI_ADJ ( K )
               DI_ADJ ( K ) = 0.0

               !------
               ! fwd code:
               ! DELC = DTS * ( MBARKS( K ) * VCI( KB,S )
               !      -   MDWN( K ) * VCI( K,S ) + DSIGH( K+1 ) 
               !      *   DSIGHI( K ) * MDWN( K+1 ) * VCI( K+1,S ) )
               ! adj code:
               VCI_ADJ ( K + 1, S ) = VCI_ADJ ( K + 1, S ) + DTS * 
     &                 DSIGH( K + 1 ) * DSIGHI( K ) * MDWN ( K + 1 )
     &                 * DELC_ADJ
               VCI_ADJ ( K, S ) = VCI_ADJ ( K, S ) - MDWN ( K ) * DTS
     &                 * DELC_ADJ
               VCI_ADJ ( KB, S ) = VCI_ADJ ( KB, S ) + DTS * MBARKS( K )
     &                 * DELC_ADJ
               DELC_ADJ = 0.0

            END DO

         END DO

      END DO

      DO S = NSP, 1, -1

         DO K = CLTOP, CLBASE, -1

            !------
            ! fwd code:
            ! VCI( K, S ) = C( S, K )
            ! adj code:
            C_ADJ ( S, K ) = C_ADJ ( S, K ) + VCI_ADJ ( K, S )
            VCI_ADJ ( K, S ) = 0.0

         END DO

         !------
         ! fwd code:
         ! VCI( KB, S ) = CBELOW( S )
         ! adj code:
         CBELOW_ADJ ( S ) = CBELOW_ADJ ( S ) + VCI_ADJ ( KB, S )
         VCI_ADJ ( KB, S ) = 0.0

      END DO

      END SUBROUTINE ACMCLD_ADJ
