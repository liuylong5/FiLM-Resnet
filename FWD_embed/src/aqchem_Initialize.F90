! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialization File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : aqchem_Initialize.F90
! Time                 : Tue Feb 26 14:38:30 2013
! Working directory    : /home/kfahey/kpp-2.2.3/aqchem_5p0_kf
! Equation file        : aqchem.kpp
! Output root filename : aqchem
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


MODULE aqchem_Initialize

!  USE aqchem_Parameters, ONLY: dp, NVAR, NFIX
  USE aqchem_Parameters
  USE aqchem_Global
  
  IMPLICIT NONE

CONTAINS

REAL( kind=dp )FUNCTION KMTF ( ACC, DG, MW )

 implicit none
 
 REAL( kind=dp ) ACC  ! accommodation coefficient     ! unitless
 REAL( kind=dp ) DG  ! Gas molecular diffusion coef !m2/s
 REAL( kind=dp ) MW  ! molecular weight
 REAL( kind=dp ) KMT, RHO1, RAD, R, PI, V
 
 RHO1 = 1.D0 !g/mL - density of H2O
 RAD = DDIAM * 0.5D0 ! m
 R = 8.3145d0
 PI = 3.1415926536
 KMT = ( RAD*RAD ) / ( 3.D0 * DG )
 V = SQRT( 8.D0 * R * TEMP_KPP * 1000.d0 / PI / MW ) !V = ( 8.D0 * R * TEMP_KPP * 1000.d0 / PI / MW )**0.5  ! m/s
 KMT = KMT + ( 4 * RAD / ( 3 * V * ACC ) )
 KMT = 1.D0 / KMT  ! from Schwartz, 1986
 KMTF = KMT * LWC_KPP / 1000.D0 ! implied division by RHO1 (=1)
 
 RETURN
 
END FUNCTION KMTF

REAL( kind=dp )FUNCTION KMTB ( HL, DH, ACC, DG, MW )

 implicit none
 
 REAL( kind=dp ) HL, DH
 REAL( kind=dp ) ACC  ! accommodation coefficient
 REAL( kind=dp ) DG  ! Gas molecular diffusion coef
 REAL( kind=dp ) MW  ! molecular weight
 REAL( kind=dp ) HLCONST, RHO1, RAD, R, PI, KMT, V, R2
 
 RAD = DDIAM * 0.5D0
 R = 8.3145d0 ! J / mol-K
 R2 = 0.08206D0 ! L-atm / mol-K  (= R/101325)
 PI = 3.1415926536
 KMT = ( RAD * RAD ) / ( 3.D0 * DG )
 V = SQRT( 8.D0 * R * TEMP_KPP * 1000.d0 / PI / MW ) ! ( 8.D0 * R * TEMP_KPP * 1000.d0 / PI / MW )**0.5  ! m/s
 KMT = KMT + ( 4 * RAD / ( 3 * V * ACC ) )
 KMT = 1.D0 / KMT  ! from Schwartz, 1986
 HLCONST = HL * EXP( DH * ( DELINVT ) )
 KMTB = KMT / ( R2 * TEMP_KPP * HLCONST )
 
 RETURN
 
END FUNCTION KMTB

REAL( kind=dp )FUNCTION ORG ( korg )

          implicit none
       
            REAL( kind=dp ) korg
         
            ORG = korg
            ORG = ORG * PHI2
         
        RETURN
     
END FUNCTION ORG


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialize - function to initialize concentrations
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!kf SUBROUTINE Initialize ( )

 SUBROUTINE Initialize ( TEMP2, PRES_PA, TAUCLD, PRCRATE,  &
   WCAVG, WTAVG, AIRM, ALFA0, ALFA3, GAS, AEROSOL, CTHK1) !, DARK )


  USE aqchem_Global
  USE AQ_DATA
  
!         INCLUDE 'new_AQ_PARAMS.EXT'      ! aqueous chemistry parameters for box

  INTEGER :: i
  REAL(kind=dp) :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!...........Arguments:


      REAL,      INTENT( IN )  :: AIRM      ! total air mass in cloudy layers (mol/m2)
      REAL,      INTENT( IN )  :: ALFA0     ! scav coef for aitken aerosol number
      REAL,      INTENT( IN )  :: ALFA3     ! scav coef for aitken aerosol mass
      REAL,      INTENT( IN )  :: PRCRATE   ! precip rate (mm/hr)
      REAL,      INTENT( IN )  :: PRES_PA   ! pressure (Pa)
      REAL,      INTENT( IN )  :: TAUCLD    ! timestep for cloud (s)
      REAL,      INTENT( IN )  :: TEMP2      ! temperature (K)
      REAL,      INTENT( IN )  :: WCAVG     ! liquid water content (kg/m3)
      REAL,      INTENT( IN )  :: WTAVG     ! total water content (kg/m3)
!      LOGICAL,   INTENT( IN )  :: DARK      ! DARK = TRUE is night,  DARK = FALSE is day     
      REAL( 8 ), INTENT( IN ) :: GAS    ( NGAS )         ! gas phase concentrations (mol/molV)
      REAL( 8 ), INTENT( IN ) :: AEROSOL( NAER, NMODES ) ! aerosol concentrations (mol/molV)
      REAL( 8 ), INTENT( IN ) :: CTHK1


!...........Local Variables:
!	Need to input and declare aqchem argument list
!	also declare CTHK1, TWASH, FE_SOL, FE_III, MN_SOL, MN_II
!	FECOR, MNCOR, NACOR, MGCOR, CACOR, KCOR, SUMPOS,
!	SUMNEG, Kw, CHGBAL

!	REAL( 8 ) :: CTHK1
	REAL( 8 ) :: TWASH
!	REAL( 8 ) :: FE_SOL 
!	REAL( 8 ) :: FE_III 
!	REAL( 8 ) :: MN_SOL 
!	REAL( 8 ) :: MN_II
!	REAL( 8 ) :: FECOR 
!	REAL( 8 ) :: MNCOR 
!	REAL( 8 ) :: NACOR 
!	REAL( 8 ) :: MGCOR 
!	REAL( 8 ) :: CACOR 
!	REAL( 8 ) :: KCOR 
	REAL( 8 ) :: SUMPOS
	REAL( 8 ) :: SUMNEG 
	REAL( 8 ) :: Kw 
!	REAL :: Kw 
	REAL( 8 ) :: CHGBAL 

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  CFACTOR = 1.000000e+00_dp

  x = (0.0)*CFACTOR
  DO i = 1, NVAR
    VAR(i) = x
  END DO

  x = (0.0)*CFACTOR
  DO i = 1, NFIX
    FIX(i) = x
  END DO

! constant rate coefficients
! END constant rate coefficients

! INLINED initializations


        TSTART = 0.0D0
        DT     = TAUCLD
        TEND = TSTART + DT
        RTOL = 1.0D-02
        ATOL = 1.0D-02

        DDIAM = 1.6D-5 ! Droplet diameter = 16 micrometers

        TEMP_KPP = TEMP2
        PRESS = PRES_PA / 101325.D0 ! Pressure (ATM)
        LWC_KPP = WCAVG
        LWCFRAC =  LWC_KPP * 1.D-3 ! L H2O / L AIR
       
        INV_TEMP      = 1.0D0 / TEMP_KPP
        DELINVT = ( 298.d0 - TEMP_KPP ) / ( 298.d0 * TEMP_KPP )
!        CTHK1 = AIRM * TEMP_KPP * 0.08206D0 / ( PRESS * 1000.0D0 ) ! cloud thickness (m)
        TWASH = WTAVG * 1000.0D0 * CTHK1 * 3600.0D0 / &
	(1000.D0 * MAX( 1.0D-20, REAL( PRCRATE, 8 )))
        IF(PRCRATE .GT. 0.0) THEN
        WETFAC_KPP = 1.D0 / TWASH
        ELSE
        WETFAC_KPP = 0.d0
        ENDIF
        ALFA3_KPP = ALFA3
        ALFA0_KPP = ALFA0
!     MGLYH = 3.2E+04  ! should depend on gas phase mech, just testing for now
        PHI2 = 1000.d0 / 6.022d23 / LWCFRAC
	INVPHI2 = 1 / PHI2

! Fraction partitioning to FE(III) and MN(II)

!        IF ( DARK ) THEN
!           FE_III = 0.9D0  ! Night time, GS 01July2011
!        ELSE
!           FE_III = 0.1D0  ! Day time, GS 01July2011
!        END IF
!
!        MN_II = 1.0D0              

! Solubility of Fe and Mn

!        FE_SOL = 0.1D0               
!        MN_SOL = 0.5D0    

! Set initial dynamic concentrations based on input gas and aerosol concentrations 
        
        VAR(ind_G_SO2) = GAS(LSO2)
        VAR(ind_G_HNO3) = GAS(LHNO3) + 2.D0*GAS(LN2O5)
        VAR(ind_G_CO2) = GAS(LCO2)
        VAR(ind_G_NH3) = GAS(LNH3)
        VAR(ind_G_H2O2) = GAS(LH2O2)
        VAR(ind_G_O3) = GAS(LO3)
        VAR(ind_G_HCOOH) = GAS(LFOA)
        VAR(ind_G_MHP) = GAS(LMHP)
        VAR(ind_G_PAA) = GAS(LPAA)
        VAR(ind_G_HCL) = GAS(LHCL)
        VAR(ind_G_GLY) = GAS(LGLY)
        VAR(ind_G_MGLY) = GAS(LMGLY)
!        VAR(ind_G_HO) = GAS(LHO)
!        VAR(ind_G_N2O5) = GAS(LN2O5)   -- added to HNO3
!        VAR(ind_G_H2SO4) = GAS(LH2SO4)  -- added to SO4MIN2
     
     VAR(ind_A_NO3AKN) = AEROSOL(LNO3,AKN)
     VAR(ind_A_NH4AKN) = AEROSOL(LNH4,AKN)
     VAR(ind_A_CLAKN) = AEROSOL(LCL,AKN)
     VAR(ind_A_NAAKN) = AEROSOL(LNA,AKN)
     VAR(ind_A_SO4AKN) = AEROSOL(LSO4,AKN)
     VAR(ind_A_PECAKN) = AEROSOL(LEC,AKN)
!slz     VAR(ind_A_POAAKN) = AEROSOL(LPOA,AKN)
     VAR(ind_A_PRIAKN) = AEROSOL(LPRI,AKN)
!     VAR(ind_A_NUMAKN) = AEROSOL(LNUM,AKN)
     
! Instantaneous droplet activation of ACC and COR modes and dissolution of N2O5 and H2SO4

     VAR(ind_L_SO4MIN2) = AEROSOL(LSO4,ACC) + AEROSOL(LSO4,COR)
     VAR(ind_L_NO3MIN) = AEROSOL(LNO3,ACC) + AEROSOL(LNO3,COR)
     VAR(ind_L_NH4PLUS) = AEROSOL(LNH4,ACC) + AEROSOL(LNH4,COR)
     VAR(ind_L_CLMIN) = AEROSOL(LCL,ACC) + AEROSOL(LCL,COR)
     VAR(ind_L_PRIACC) = AEROSOL(LPRI,ACC)
     VAR(ind_L_NAPLUS) = AEROSOL(LNA,ACC)
!     VAR(ind_L_CAPLUS2) = AEROSOL(LCAACC,ACC)
!     VAR(ind_L_MGPLUS2) = AEROSOL(LMGACC,ACC)
!     VAR(ind_L_KPLUS) = AEROSOL(LKACC,ACC)
     VAR(ind_L_PECACC) = AEROSOL(LEC,ACC)
     VAR(ind_L_ORGC) = AEROSOL(LORGC,ACC)
!slz     VAR(ind_L_POAACC) = AEROSOL(LPOA,ACC)
   
     VAR(ind_L_SO4MIN2) = VAR(ind_L_SO4MIN2) + GAS(LH2SO4)
     
     VAR(ind_L_CO3MIN2) = AEROSOL(LCACO3,COR) + AEROSOL(LMGCO3,COR)
     VAR(ind_L_CAPLUS2) = AEROSOL(LCACO3,COR)
     VAR(ind_L_MGPLUS2) = AEROSOL(LMGCO3,COR)
     VAR(ind_L_KPLUS) = AEROSOL(LK,COR)

! Coarse crustal species from SOILCOR, ANTHCOR, SEASCOR

!        FECOR   = 0.0281D0  * ( 100.0D0 / 55.8D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0467D0  * ( 100.0D0 / 55.8D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )
!        MNCOR   = 0.00078D0 * ( 100.0D0 / 54.9D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0011D0  * ( 100.0D0 / 54.9D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )
!        NACOR   = 0.8373D0  * (  23.0D0 / 23.0D0 ) * AEROSOL(LSEASC,COR) &                    
!                + 0.0652D0  * ( 100.0D0 / 23.0D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0023D0  * ( 100.0D0 / 23.0D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )
!        MGCOR   = 0.0997D0  * (  23.0D0 / 24.3D0 ) * AEROSOL(LSEASC,COR) &                    
!                + 0.0000D0  * ( 100.0D0 / 24.3D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0032D0  * ( 100.0D0 / 24.3D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )
!        CACOR   = 0.0320D0  * (  23.0D0 / 40.1D0 ) * AEROSOL(LSEASC,COR) &                 
!                + 0.0872D0  * ( 100.0D0 / 40.1D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0562D0  * ( 100.0D0 / 40.1D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )
!        KCOR    = 0.0310D0  * (  23.0D0 / 39.1D0 ) * AEROSOL(LSEASC,COR) &                 
!                + 0.0252D0  * ( 100.0D0 / 39.1D0 ) * AEROSOL(LSOILC,COR) / ( 1.0 - 0.04642 ) &
!                + 0.0176D0  * ( 100.0D0 / 39.1D0 ) * AEROSOL(LANTHC,COR) / ( 1.0 - 0.00325 )

! Fe3+ and Mn2+ for SIV oxidation
               
!         VAR(ind_L_FEPLUS3) = FE_SOL * FE_III * ( AEROSOL(LFEACC,ACC) + FECOR )     
!         VAR(ind_L_MNPLUS2) = MN_SOL * MN_II * ( AEROSOL(LMNACC,ACC) + MNCOR )

         VAR(ind_L_FEPLUS3) = AEROSOL(LA3FE,COR)   
         VAR(ind_L_MNPLUS2) = AEROSOL(LB2MN,COR)
      
         VAR(ind_L_NAPLUS) = VAR(ind_L_NAPLUS) + AEROSOL(LNA,COR) !+ NACOR
!         VAR(ind_L_CAPLUS2) = VAR(ind_L_CAPLUS2) + CACOR
!         VAR(ind_L_MGPLUS2) = VAR(ind_L_MGPLUS2) + MGCOR
!         VAR(ind_L_KPLUS) = VAR(ind_L_KPLUS) + KCOR

! Convert concententrations from input units mol/mol air --> molec/cm3 air     

     CFACTOR = PRES_PA / (TEMP_KPP * 8.314510) * PHI

     DO i = 1, NVAR
        VAR(i) = CFACTOR * VAR(i)
     END DO

     DO i = 1, NFIX
        FIX(i) = CFACTOR * FIX(i)
     END DO
    
!kf      FIX(indf_L_H2O) = 1.D0 * INVPHI2   ! set to 1 M -- Kw includes [H2O]
      FIX(indf_L_H2O) = 55.5D0 * INVPHI2   ! 
      FIX(indf_L_HO) = ( GAS(LHO) * PRESS * 30.D0 * EXP( 4.5D3 * &
      DELINVT ) ) * INVPHI2
       
   
!  calculate initial H+ and OH- from electroneutrality and Kw
!
!  Sum of positive ions + H+ = Sum of negative ions + OH-
!
!  [H+] * [OH-] = Kw = 1.00E-14 * EXP( -6.71E+03 * ((1.d0 / TEMP_KPP) - (1.D0 / 298.D0)) )     ! Smith and Martell (1976)
!  [H+] = Kw / [OH-]
!
!  SUMPOS + Kw / [OH-] = SUMNEG + [OH-] -->  [OH-]**2 - (SUMPOS - SUMNEG)*[OH-] - Kw
!  Solve for [OH-] with quadratic formula and plug back into Kw relation to get initial [H+]

!        Kw = 1.00E-14 * EXP( -6.71E+03 * ((1.d0 / TEMP_KPP) - (1.D0 / 298.D0)))
        Kw = 1.00D-14 * EXP( -6.955D+03 * (DELINVT))  ! includes conc of H2O

        SUMPOS = 2.D0 * (VAR(ind_L_CAPLUS2) + VAR(ind_L_MGPLUS2) ) + VAR(ind_L_NAPLUS) &
                 + VAR(ind_L_KPLUS) + VAR(ind_L_NH4PLUS)
        SUMNEG = 2.D0 * (VAR(ind_L_SO4MIN2) + VAR(ind_L_CO3MIN2)) + VAR(ind_L_NO3MIN) &
	         + VAR(ind_L_CLMIN)
    
     SUMPOS = SUMPOS * PHI2
     SUMNEG = SUMNEG * PHI2
         
        CHGBAL = SUMPOS - SUMNEG
    
        VAR(ind_L_OHMIN) = (CHGBAL + SQRT(CHGBAL*CHGBAL + 4.d0 * Kw)) * 0.5D0
        VAR(ind_L_HPLUS) = Kw / VAR(ind_L_OHMIN)

! Does that equal CHGBAL - OHMIN?  
!    write(6,*) VAR(ind_L_HPLUS), VAR(ind_L_OHMIN) - CHGBAL 

!	write(6,*) -DLOG10(VAR(ind_L_HPLUS)) 
    
        IF (VAR(ind_L_OHMIN) .LT. 0.d0) THEN
           print *, 'NEGATIVE INITIAL OHMIN+HPLUS CONC'
           stop
        ENDIF
     
       VAR(ind_L_OHMIN) = VAR(ind_L_OHMIN) * INVPHI2
       VAR(ind_L_HPLUS) = VAR(ind_L_HPLUS) * INVPHI2

  RCONST(1) = ((KMTF(0.11D0,1.28D-5,64.064D0)))
  RCONST(2) = ((KMTF(0.0868D0,1.32D-5,63.013D0)))
  RCONST(3) = ((KMTF(0.00015D0,1.55D-5,44.01D0)))
  RCONST(4) = ((KMTF(0.091D0,2.3D-5,17.031D0)))
  RCONST(5) = ((KMTF(0.1532D0,1.46D-5,34.015D0)))
  RCONST(6) = ((KMTF(0.1D0,1.48D-5,47.998D0)))
  RCONST(7) = ((KMTF(0.0229D0,1.53D-5,46.025D0)))
  RCONST(8) = ((KMTF(0.006758D0,1.31D-5,48.04D0)))
  RCONST(9) = ((KMTF(0.019D0,1.02D-5,76.05D0)))
  RCONST(10) = ((KMTF(0.1158D0,1.89D-5,36.461D0)))
  RCONST(11) = ((KMTF(0.023D0,1.15D-5,58.04D0)))
  RCONST(12) = ((KMTF(0.023D0,1.15D-5,72.06D0)))
  RCONST(13) = ((KMTB(1.4D+00,2.9D+03,0.11D0,1.28D-5,64.064D0)))
  RCONST(14) = ((KMTB(2.1D+05,8.7D+03,0.0868D0,1.32D-5,63.013D0)))
  RCONST(15) = ((KMTB(3.6D-02,2.2D+03,0.00015D0,1.55D-5,44.01D0)))
  RCONST(16) = ((KMTB(6.1D+01,4.2D+03,0.091D0,2.3D-5,17.031D0)))
  RCONST(17) = ((KMTB(8.3D+04,7.4D+03,0.1532D0,1.46D-5,34.015D0)))
  RCONST(18) = ((KMTB(1.14D-02,2.3D+03,0.1D0,1.48D-5,47.998D0)))
  RCONST(19) = ((KMTB(8.9D+03,6.1D+03,0.0229D0,1.53D-5,46.025D0)))
  RCONST(20) = ((KMTB(3.1D+02,5.2D+03,0.006758D0,1.31D-5,48.04D0)))
  RCONST(21) = ((KMTB(8.4D+02,5.3D+03,0.019D0,1.02D-5,76.05D0)))
  RCONST(22) = ((KMTB(1.9D+01,6.0D+02,0.1158D0,1.89D-5,36.461D0)))
  RCONST(23) = ((KMTB(3.6D+05,0.0D+00,0.023D0,1.15D-5,58.04D0)))
  RCONST(24) = ((KMTB(MGLYH,0.0D+0,0.023D0,1.15D-5,72.06D0)))
  RCONST(25) = (ALFA3_KPP)
  RCONST(26) = (ALFA3_KPP)
  RCONST(27) = (ALFA3_KPP)
  RCONST(28) = (ALFA3_KPP)
  RCONST(29) = (ALFA3_KPP)
  RCONST(30) = (ALFA3_KPP)
  RCONST(31) = (ALFA3_KPP)
  RCONST(32) = (ALFA3_KPP)
  RCONST(55) = ((ORG(3.0D10)))
  RCONST(56) = ((ORG(3.0D10)))    
  RCONST(73) = (WETFAC_KPP)
  RCONST(74) = (WETFAC_KPP)
  RCONST(75) = (WETFAC_KPP)
  RCONST(76) = (WETFAC_KPP)
  RCONST(77) = (WETFAC_KPP)
  RCONST(78) = (WETFAC_KPP)
  RCONST(79) = (WETFAC_KPP)
  RCONST(80) = (WETFAC_KPP)
  RCONST(81) = (WETFAC_KPP)
  RCONST(82) = (WETFAC_KPP)
  RCONST(83) = (WETFAC_KPP)
  RCONST(84) = (WETFAC_KPP)
  RCONST(85) = (WETFAC_KPP)
  RCONST(86) = (WETFAC_KPP)
  RCONST(87) = (WETFAC_KPP)
  RCONST(88) = (WETFAC_KPP)
  RCONST(89) = (WETFAC_KPP)
  RCONST(90) = (WETFAC_KPP)
  RCONST(91) = (WETFAC_KPP)
  RCONST(92) = (WETFAC_KPP)
  RCONST(93) = (WETFAC_KPP)
  RCONST(94) = (WETFAC_KPP)
  RCONST(95) = (WETFAC_KPP)
  RCONST(96) = (WETFAC_KPP)
  RCONST(97) = (WETFAC_KPP)
  RCONST(98) = (WETFAC_KPP)
  RCONST(99) = (WETFAC_KPP)
  RCONST(100) = (WETFAC_KPP)
  RCONST(101) = (WETFAC_KPP)
  RCONST(102) = (WETFAC_KPP)
  RCONST(103) = (WETFAC_KPP)
  RCONST(104) = (WETFAC_KPP)
  RCONST(105) = (WETFAC_KPP)
  RCONST(106) = (WETFAC_KPP)
  RCONST(107) = (WETFAC_KPP)


! End INLINED initializations

      
END SUBROUTINE Initialize

! End of Initialize function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE aqchem_Initialize

