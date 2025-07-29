      subroutine aqchem_newton(BB,
     +     RECIPA1,RECIPA2,PSO2F,PNH3F,PHCLF,PFOAF,PHNO3F,PCO2F,
     +     AE, SO4, HSO4, SO3, HSO3, CO3, HCO3, OH, NH4, HCO2, NO3, CL,
     +     INITGAS, SK6TS6, TS6, NA, CA, MG,
     +     MN, FE, K, A, B, 
     +     XLSO2, SO21, SO212, SO212H, SO21H,
     +     XLNH3, NH3DH20, NH31HDH,
     +     XLHCL, HCL1,HCL1H,
     +     XL, FOAH, FOA1H,
     +     XLHNO3, HNO31, HNO31H,
     +     XLCO2, CO21, CO212, CO212H, CO21H, H2OW,
     +     ACT1, ACT2, GM2, SK6,
     +     NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2)
      IMPLICIT NONE
      REAL( 8 ) :: BB              !Ph 
      REAL( 8 ) :: AE              ! guess for H+ conc in cloudwater (mol/liter)
      REAL( 8 ) :: SO4             ! SO4= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HSO4            ! HSO4 concn in cloudwater (mol/liter)
      REAL( 8 ) :: SO3             ! SO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HSO3            ! HSO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: CO3             ! CO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCO3            ! HCO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: OH              ! OH conc in cloudwater (mol/liter)
      REAL( 8 ) :: NH4             ! NH4+ conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCO2            ! HCO2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NO3             ! NO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: CL              ! total Cl-  conc in cloudwater (mol/liter)

      REAL( 8 ) :: MN, FE, K, A, B

      REAL( 8 ) :: INITGAS( NGAS ) ! initial gas partial pressure (atm)      
      REAL( 8 ) :: SK6TS6          !
      REAL( 8 ) :: TS6             ! SO4 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NA              ! Na conc in cloudwater (mol/liter)
      REAL( 8 ) :: CA              ! Calcium conc in cloudwater (mol/liter)
      REAL( 8 ) :: MG              !

      REAL( 8 ) :: RECIPA1         !
      REAL( 8 ) :: RECIPA2         !
      REAL( 8 ) :: PSO2F           ! gas only SO2 partial pressure (atm)
      REAL( 8 ) :: PNH3F           ! gas only NH3 partial pressure (atm)
      REAL( 8 ) :: PHCLF           ! gas only HCL partial pressure (atm)
      REAL( 8 ) :: PFOAF           ! gas only ORGANIC ACID partial press (atm)
      REAL( 8 ) :: PHNO3F          ! gas only HNO3 partial pressure (atm)
      REAL( 8 ) :: PCO2F           ! gas only CO2 partial pressure (atm)

      REAL( 8 ) :: XLSO2           !
      REAL( 8 ) :: XLNH3           !
      REAL( 8 ) :: XLHCL           ! const in calc of HCL final partial pres
      REAL( 8 ) :: XL              !
      REAL( 8 ) :: XLHNO3          !
      REAL( 8 ) :: XLCO2           !

C...dissociation constants
      REAL( 8 ) :: SO21            ! First dissociation constant for SO2
      REAL( 8 ) :: SO21H           ! SO21*SO2H
      REAL( 8 ) :: SO212           ! SO21*SO22
      REAL( 8 ) :: SO212H          ! SO21*SO22*SO2H
      REAL( 8 ) :: NH3DH20         !
      REAL( 8 ) :: NH31HDH         !
      REAL( 8 ) :: HCL1            ! First dissociation constant for HCL
      REAL( 8 ) :: HCL1H           ! HCL1*HCLH
      REAL( 8 ) :: FOAH            ! Henry's Law constant for FOA
      REAL( 8 ) :: FOA1H           ! FOAH*FOA1
      REAL( 8 ) :: HNO31           ! First dissociation constant for HNO3
      REAL( 8 ) :: HNO31H          !
      REAL( 8 ) :: CO21            ! First dissociation constant for CO2
      REAL( 8 ) :: CO21H           ! CO2H*CO21
      REAL( 8 ) :: CO212           ! CO21*CO22
      REAL( 8 ) :: CO212H          ! CO2H*CO21*CO22
      REAL( 8 ) :: SK6             !
      REAL( 8 ) :: H2OW            !

      REAL( 8 ) :: ACT1            ! activity correction factor!single ions
      REAL( 8 ) :: ACT2            ! activity factor correction!double ions
      REAL( 8 ) :: GM2             ! activity correction factor
!
      INTEGER  NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2
!local variables
      REAL( 8 ), PARAMETER :: SMALL = 1.D-6
      REAL( 8 ) :: FB              !
      REAL( 8 ) :: FBP             !
      REAL( 8 ) :: BBP             !
      REAL( 8 ) :: STEP
      REAL( 8 ) :: ALPHA
      REAL( 8 ) :: STION, GM1LOG, GM2LOG, GM1

      BBP=1.D0
      CALL AQCHEM_FUNp(bb, bbp, fb, fbp, ae, so4, hso4, so3, hso3, co3
     +  , hco3, oh, nh4, hco2, no3, cl, initgas, sk6ts6, ts6, na, ca, mg,
     +  xlso2, so21, so212, so212h, so21h, xlnh3, nh3dh20, nh31hdh, xlhcl,
     +  hcl1, hcl1h, xl, foah, foa1h, xlhno3, hno31, hno31h, xlco2, co21,
     +  co212, co212h, co21h, h2ow, act1, act2, gm2, sk6, ngas, lso2, lnh3,
     +  lhcl, lfoa, lhno3, lco2)

!slz limit the Newton step size
      STEP=-FB/FBP
      ALPHA=.1D0
      DO WHILE ( ABS(ALPHA*STEP).GT.SMALL)
       ALPHA=.1D0*ALPHA
      END DO
      STEP=ALPHA*STEP
      BB=BB+STEP

      CALL aqchem_fun(BB,
     +     RECIPA1,RECIPA2,PSO2F,PNH3F,PHCLF,PFOAF,PHNO3F,PCO2F,
     +     AE, SO4, HSO4, SO3, HSO3, CO3, HCO3, OH, NH4, HCO2, NO3, CL,
     +     INITGAS, SK6TS6, TS6, NA, CA, MG,
     +     XLSO2, SO21, SO212, SO212H, SO21H,
     +     XLNH3, NH3DH20, NH31HDH,
     +     XLHCL, HCL1,HCL1H,
     +     XL, FOAH, FOA1H,
     +     XLHNO3, HNO31, HNO31H,
     +     XLCO2, CO21, CO212, CO212H, CO21H, H2OW,
     +     ACT1, ACT2, GM2, SK6,
     +     NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2)

!slz update ACT1, ACT2, GM2
          STION = 0.5D0 * (AE + NH4 + OH + HCO3 + HSO3
     &          + 4.0D0 * (SO4 + CO3 + SO3 + CA + MG + MN)
     &          + NO3 + HSO4 + 9.0D0 * FE + NA + K + CL + A + B + HCO2)
          GM1LOG = -0.509D0 * ( SQRT( STION )
     &           / ( 1.0D0 + SQRT( STION ) ) - 0.2D0 * STION )
          GM2LOG = GM1LOG * 4.0D0
          GM1  = 10.0D0 ** GM1LOG
          GM2  = MAX( 10.0D0 ** GM2LOG, 1.0D-30 )
          ACT1 = MAX( GM1 * GM1, 1.0D-30 )
          ACT2 = MAX( GM1 * GM1 * GM2, 1.0D-30 )
!slz update gas partial pressures and aerosol liquid phase concentrations
      CALL aqchem_fun(BB,
     +     RECIPA1,RECIPA2,PSO2F,PNH3F,PHCLF,PFOAF,PHNO3F,PCO2F,
     +     AE, SO4, HSO4, SO3, HSO3, CO3, HCO3, OH, NH4, HCO2, NO3, CL,
     +     INITGAS, SK6TS6, TS6, NA, CA, MG,
     +     XLSO2, SO21, SO212, SO212H, SO21H,
     +     XLNH3, NH3DH20, NH31HDH,
     +     XLHCL, HCL1,HCL1H,
     +     XL, FOAH, FOA1H,
     +     XLHNO3, HNO31, HNO31H,
     +     XLCO2, CO21, CO212, CO212H, CO21H, H2OW,
     +     ACT1, ACT2, GM2, SK6,
     +     NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2)

!slz update GM2
          STION = 0.5D0 * (AE + NH4 + OH + HCO3 + HSO3
     &          + 4.0D0 * (SO4 + CO3 + SO3 + CA + MG + MN)
     &          + NO3 + HSO4 + 9.0D0 * FE + NA + K + CL + A + B + HCO2)
          GM1LOG = -0.509D0 * ( SQRT( STION )
     &           / ( 1.0D0 + SQRT( STION ) ) - 0.2D0 * STION )
          GM2LOG = GM1LOG * 4.0D0
          GM2  = MAX( 10.0D0 ** GM2LOG, 1.0D-30 )

      end subroutine

!        Generated by TAPENADE     (INRIA, Tropics team)
!slz post-processed by tpn2cmaq
!  Tapenade 3.6 (r4343) - 10 Feb 2012 10:52
!
!  Differentiation of aqchem_fun in forward (tangent) mode:
!   variations   of useful results: fb
!   with respect to varying inputs: bb
!   RW status of diff variables: bb:in fb:out
      SUBROUTINE AQCHEM_FUNp(bb, bbd, fb, fbd, ae, so4, hso4, so3, hso3, co3
     +  , hco3, oh, nh4, hco2, no3, cl, initgas, sk6ts6, ts6, na, ca, mg,
     +  xlso2, so21, so212, so212h, so21h, xlnh3, nh3dh20, nh31hdh, xlhcl,
     +  hcl1, hcl1h, xl, foah, foa1h, xlhno3, hno31, hno31h, xlco2, co21,
     +  co212, co212h, co21h, h2ow, act1, act2, gm2, sk6, ngas, lso2, lnh3,
     +  lhcl, lfoa, lhno3, lco2)
        IMPLICIT NONE
!
! lower limit guess of cloudwater pH
        REAL*8 :: bb
        REAL*8 :: bbd
! functional value
        REAL*8 :: fb
        REAL*8 :: fbd
! guess for H+ conc in cloudwater (mol/liter)
        REAL*8 :: ae
        REAL*8 :: aed
! SO4= conc in cloudwater (mol/liter)
        REAL*8 :: so4
        REAL*8 :: so4d
! HSO4 concn in cloudwater (mol/liter)
        REAL*8 :: hso4
        REAL*8 :: hso4d
! SO3= conc in cloudwater (mol/liter)
        REAL*8 :: so3
        REAL*8 :: so3d
! HSO3 conc in cloudwater (mol/liter)
        REAL*8 :: hso3
        REAL*8 :: hso3d
! CO3= conc in cloudwater (mol/liter)
        REAL*8 :: co3
        REAL*8 :: co3d
! HCO3 conc in cloudwater (mol/liter)
        REAL*8 :: hco3
        REAL*8 :: hco3d
! OH conc in cloudwater (mol/liter)
        REAL*8 :: oh
        REAL*8 :: ohd
! NH4+ conc in cloudwater (mol/liter)
        REAL*8 :: nh4
        REAL*8 :: nh4d
! HCO2 conc in cloudwater (mol/liter)
        REAL*8 :: hco2
        REAL*8 :: hco2d
! NO3 conc in cloudwater (mol/liter)
        REAL*8 :: no3
        REAL*8 :: no3d
! total Cl-  conc in cloudwater (mol/liter)
        REAL*8 :: cl
        REAL*8 :: cld
!
        INTEGER :: ngas, lso2, lnh3, lhcl, lfoa, lhno3, lco2
!
! initial gas partial pressure (atm)
        REAL*8 :: initgas(ngas)
!
        REAL*8 :: sk6ts6
! SO4 conc in cloudwater (mol/liter)
        REAL*8 :: ts6
! Na conc in cloudwater (mol/liter)
        REAL*8 :: na
! Calcium conc in cloudwater (mol/liter)
        REAL*8 :: ca
!
        REAL*8 :: mg
!
!
        REAL*8 :: xlso2
!
        REAL*8 :: xlnh3
! const in calc of HCL final partial pres
        REAL*8 :: xlhcl
!
        REAL*8 :: xl
!
        REAL*8 :: xlhno3
!
        REAL*8 :: xlco2
!
!...dissociation constants
! First dissociation constant for SO2
        REAL*8 :: so21
! SO21*SO2H
        REAL*8 :: so21h
! SO21*SO22
        REAL*8 :: so212
! SO21*SO22*SO2H
        REAL*8 :: so212h
!
        REAL*8 :: nh3dh20
!
        REAL*8 :: nh31hdh
! First dissociation constant for HCL
        REAL*8 :: hcl1
! HCL1*HCLH
        REAL*8 :: hcl1h
! Henry's Law constant for FOA
        REAL*8 :: foah
! FOAH*FOA1
        REAL*8 :: foa1h
! First dissociation constant for HNO3
        REAL*8 :: hno31
!
        REAL*8 :: hno31h
! First dissociation constant for CO2
        REAL*8 :: co21
! CO2H*CO21
        REAL*8 :: co21h
! CO21*CO22
        REAL*8 :: co212
! CO2H*CO21*CO22
        REAL*8 :: co212h
!
        REAL*8 :: sk6
!
        REAL*8 :: h2ow
!
! activity correction factor!single ions
        REAL*8 :: act1
! activity factor correction!double ions
        REAL*8 :: act2
! activity correction factor
        REAL*8 :: gm2
!local variables
!
        REAL*8 :: recipa1
        REAL*8 :: recipa1d
!
        REAL*8 :: recipa2
        REAL*8 :: recipa2d
! gas only SO2 partial pressure (atm)
        REAL*8 :: pso2f
        REAL*8 :: pso2fd
! gas only NH3 partial pressure (atm)
        REAL*8 :: pnh3f
        REAL*8 :: pnh3fd
! gas only HCL partial pressure (atm)
        REAL*8 :: phclf
        REAL*8 :: phclfd
! gas only ORGANIC ACID partial press (atm)
        REAL*8 :: pfoaf
        REAL*8 :: pfoafd
! gas only HNO3 partial pressure (atm)
        REAL*8 :: phno3f
        REAL*8 :: phno3fd
! gas only CO2 partial pressure (atm)
        REAL*8 :: pco2f
        REAL*8 :: pco2fd
        REAL*8 :: pwy1
        REAL*8 :: pwy1d
!
!
        pwy1d = -bbd
        pwy1 = -bb
        aed = 10.0d0**pwy1*LOG(10.0d0)*pwy1d
        ae = 10.0d0**pwy1
        recipa1d = -(act1*aed/(ae*act1)**2)
        recipa1 = 1.0d0/(ae*act1)
        recipa2d = -(act2*(aed*ae+ae*aed)/(ae*ae*act2)**2)
        recipa2 = 1.0d0/(ae*ae*act2)
!
!...calculate final gas phase partial pressure of SO2, NH3, HNO3
!...  HCOOH, and CO2 (atm)
!
        pso2fd = -(initgas(lso2)*xlso2*(so21*recipa1d+so212*recipa2d)/(1.0d0+
     +    xlso2*(1.0d0+so21*recipa1+so212*recipa2))**2)
        pso2f = initgas(lso2)/(1.0d0+xlso2*(1.0d0+so21*recipa1+so212*recipa2))
!
        pnh3fd = -(initgas(lnh3)*xlnh3*nh3dh20*aed/(1.0d0+xlnh3*(1.0d0+nh3dh20
     +    *ae))**2)
        pnh3f = initgas(lnh3)/(1.0d0+xlnh3*(1.0d0+nh3dh20*ae))
!
        phclfd = -(initgas(lhcl)*xlhcl*hcl1*recipa1d/(1.0d0+xlhcl*(1.0d0+hcl1*
     +    recipa1))**2)
        phclf = initgas(lhcl)/(1.0d0+xlhcl*(1.0d0+hcl1*recipa1))
!
        pfoafd = -(initgas(lfoa)*xl*foa1h*recipa1d/(1.0d0+xl*(foah+foa1h*
     +    recipa1))**2)
        pfoaf = initgas(lfoa)/(1.0d0+xl*(foah+foa1h*recipa1))
!
        phno3fd = -(initgas(lhno3)*xlhno3*hno31*recipa1d/(1.0d0+xlhno3*(1.0d0+
     +    hno31*recipa1))**2)
        phno3f = initgas(lhno3)/(1.0d0+xlhno3*(1.0d0+hno31*recipa1))
!
        pco2fd = -(initgas(lco2)*xlco2*(co21*recipa1d+co212*recipa2d)/(1.0d0+
     +    xlco2*(1.0d0+co21*recipa1+co212*recipa2))**2)
        pco2f = initgas(lco2)/(1.0d0+xlco2*(1.0d0+co21*recipa1+co212*recipa2))
!
!...calculate liquid phase concentrations (moles/liter)
!
        so4d = -(sk6ts6*gm2*aed/(ae*gm2+sk6)**2)
        so4 = sk6ts6/(ae*gm2+sk6)
        hso4d = -so4d
        hso4 = ts6 - so4
        so3d = so212h*(pso2fd*recipa2+pso2f*recipa2d)
        so3 = so212h*pso2f*recipa2
        hso3d = so21h*(pso2fd*recipa1+pso2f*recipa1d)
        hso3 = so21h*pso2f*recipa1
        co3d = co212h*(pco2fd*recipa2+pco2f*recipa2d)
        co3 = co212h*pco2f*recipa2
        hco3d = co21h*(pco2fd*recipa1+pco2f*recipa1d)
        hco3 = co21h*pco2f*recipa1
        ohd = h2ow*recipa1d
        oh = h2ow*recipa1
        nh4d = nh31hdh*(pnh3fd*ae+pnh3f*aed)
        nh4 = nh31hdh*pnh3f*ae
        hco2d = foa1h*(pfoafd*recipa1+pfoaf*recipa1d)
        hco2 = foa1h*pfoaf*recipa1
        no3d = hno31h*(phno3fd*recipa1+phno3f*recipa1d)
        no3 = hno31h*phno3f*recipa1
! new for sea salt
        cld = hcl1h*(phclfd*recipa1+phclf*recipa1d)
        cl = hcl1h*phclf*recipa1
!
!...compute functional value
!
        fbd = aed + nh4d + 2.0d0*(-co3d-so3d-so4d) - ohd - hco3d - hso3d -
     +    no3d - hso4d - hco2d - cld
        fb = ae + nh4 + na + 2.0d0*(ca+mg-co3-so3-so4) - oh - hco3 - hso3 -
     +    no3 - hso4 - hco2 - cl
      END SUBROUTINE AQCHEM_FUNp

      subroutine aqchem_fun(BB,
     +     RECIPA1,RECIPA2,PSO2F,PNH3F,PHCLF,PFOAF,PHNO3F,PCO2F,
     +     AE, SO4, HSO4, SO3, HSO3, CO3, HCO3, OH, NH4, HCO2, NO3, CL,
     +     INITGAS, SK6TS6, TS6, NA, CA, MG,
     +     XLSO2, SO21, SO212, SO212H, SO21H,
     +     XLNH3, NH3DH20, NH31HDH,
     +     XLHCL, HCL1,HCL1H,
     +     XL, FOAH, FOA1H,
     +     XLHNO3, HNO31, HNO31H,
     +     XLCO2, CO21, CO212, CO212H, CO21H, H2OW,
     +     ACT1, ACT2, GM2, SK6,
     +     NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2)
      IMPLICIT NONE
    
      REAL( 8 ) :: BB              ! lower limit guess of cloudwater pH
      REAL( 8 ) :: AE              ! guess for H+ conc in cloudwater (mol/liter)
      REAL( 8 ) :: SO4             ! SO4= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HSO4            ! HSO4 concn in cloudwater (mol/liter)
      REAL( 8 ) :: SO3             ! SO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HSO3            ! HSO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: CO3             ! CO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCO3            ! HCO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: OH              ! OH conc in cloudwater (mol/liter)
      REAL( 8 ) :: NH4             ! NH4+ conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCO2            ! HCO2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NO3             ! NO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: CL              ! total Cl-  conc in cloudwater (mol/liter)

      REAL( 8 ) :: INITGAS( NGAS ) ! initial gas partial pressure (atm)      
      REAL( 8 ) :: SK6TS6          !
      REAL( 8 ) :: TS6             ! SO4 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NA              ! Na conc in cloudwater (mol/liter)
      REAL( 8 ) :: CA              ! Calcium conc in cloudwater (mol/liter)
      REAL( 8 ) :: MG              !

      REAL( 8 ) :: XLSO2           !
      REAL( 8 ) :: XLNH3           !
      REAL( 8 ) :: XLHCL           ! const in calc of HCL final partial pres
      REAL( 8 ) :: XL              !
      REAL( 8 ) :: XLHNO3          !
      REAL( 8 ) :: XLCO2           !

C...dissociation constants
      REAL( 8 ) :: SO21            ! First dissociation constant for SO2
      REAL( 8 ) :: SO21H           ! SO21*SO2H
      REAL( 8 ) :: SO212           ! SO21*SO22
      REAL( 8 ) :: SO212H          ! SO21*SO22*SO2H
      REAL( 8 ) :: NH3DH20         !
      REAL( 8 ) :: NH31HDH         !
      REAL( 8 ) :: HCL1            ! First dissociation constant for HCL
      REAL( 8 ) :: HCL1H           ! HCL1*HCLH
      REAL( 8 ) :: FOAH            ! Henry's Law constant for FOA
      REAL( 8 ) :: FOA1H           ! FOAH*FOA1
      REAL( 8 ) :: HNO31           ! First dissociation constant for HNO3
      REAL( 8 ) :: HNO31H          !
      REAL( 8 ) :: CO21            ! First dissociation constant for CO2
      REAL( 8 ) :: CO21H           ! CO2H*CO21
      REAL( 8 ) :: CO212           ! CO21*CO22
      REAL( 8 ) :: CO212H          ! CO2H*CO21*CO22
      REAL( 8 ) :: SK6             !
      REAL( 8 ) :: H2OW            !

      REAL( 8 ) :: ACT1            ! activity correction factor!single ions
      REAL( 8 ) :: ACT2            ! activity factor correction!double ions
      REAL( 8 ) :: GM2             ! activity correction factor
!
      INTEGER  NGAS, LSO2,LNH3,LHCL,LFOA,LHNO3,LCO2
!local variables
      REAL( 8 ) :: RECIPA1         !
      REAL( 8 ) :: RECIPA2         !
      REAL( 8 ) :: PSO2F           ! gas only SO2 partial pressure (atm)
      REAL( 8 ) :: PNH3F           ! gas only NH3 partial pressure (atm)
      REAL( 8 ) :: PHCLF           ! gas only HCL partial pressure (atm)
      REAL( 8 ) :: PFOAF           ! gas only ORGANIC ACID partial press (atm)
      REAL( 8 ) :: PHNO3F          ! gas only HNO3 partial pressure (atm)
      REAL( 8 ) :: PCO2F           ! gas only CO2 partial pressure (atm)


      AE = 10.0D0 ** ( -BB )
      RECIPA1 = 1.0D0 / ( AE * ACT1 )
      RECIPA2 = 1.0D0 / ( AE * AE * ACT2 )

C...calculate final gas phase partial pressure of SO2, NH3, HNO3
C...  HCOOH, and CO2 (atm)

      PSO2F = INITGAS( LSO2 ) / ( 1.0D0 + XLSO2 * ( 1.0D0 + SO21 * RECIPA1
     &        + SO212 * RECIPA2 ) )

      PNH3F = INITGAS( LNH3 ) / ( 1.0D0 + XLNH3 * ( 1.0D0 + NH3DH20 * AE ) )

      PHCLF = INITGAS( LHCL ) / ( 1.0D0 + XLHCL *  ( 1.0D0 + HCL1 * RECIPA1 ) )

      PFOAF = INITGAS( LFOA ) / ( 1.0D0 + XL * ( FOAH + FOA1H * RECIPA1 ) )

      PHNO3F = INITGAS( LHNO3 ) / ( 1.0D0 + XLHNO3 * ( 1.0D0 + HNO31 * RECIPA1 ) )

      PCO2F = INITGAS( LCO2 ) / ( 1.0D0 + XLCO2 * ( 1.0D0 + CO21 * RECIPA1
     &        + CO212 * RECIPA2 ) )

C...calculate liquid phase concentrations (moles/liter)

      SO4  = SK6TS6 / ( AE * GM2 + SK6 )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2F  * RECIPA2
      HSO3 = SO21H   * PSO2F  * RECIPA1
      CO3  = CO212H  * PCO2F  * RECIPA2
      HCO3 = CO21H   * PCO2F  * RECIPA1
      OH   = H2OW    * RECIPA1
      NH4  = NH31HDH * PNH3F  * AE
      HCO2 = FOA1H   * PFOAF  * RECIPA1
      NO3  = HNO31H  * PHNO3F * RECIPA1
      CL   = HCL1H   * PHCLF  * RECIPA1 ! new for sea salt

      end subroutine
