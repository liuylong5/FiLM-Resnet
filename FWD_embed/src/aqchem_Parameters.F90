! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
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
! File                 : aqchem_Parameters.F90
! Time                 : Tue Dec  2 14:12:20 2014
! Working directory    : /home/kfahey/kpp-2.2.3/ASSUMPTION_TEST/CLEAN_INIT/TESTS/CMAQ_BASE_022414/3D/FINAL_120114
! Equation file        : aqchem.kpp
! Output root filename : aqchem
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE aqchem_Parameters

  USE aqchem_Precision
  PUBLIC
  SAVE


! TOTAL_SPECIES - Total Number of species
  INTEGER, PARAMETER :: TOTAL_SPECIES = 85 
! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 85 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 82 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 55 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 3 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 107 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 83 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 242 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 243 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 83 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 85 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_A_NO3AKN = 1 
  INTEGER, PARAMETER :: ind_A_NH4AKN = 2 
  INTEGER, PARAMETER :: ind_A_CLAKN = 3 
  INTEGER, PARAMETER :: ind_A_NAAKN = 4 
  INTEGER, PARAMETER :: ind_A_SO4AKN = 5 
  INTEGER, PARAMETER :: ind_A_PECAKN = 6 
  INTEGER, PARAMETER :: ind_A_POAAKN = 7 
  INTEGER, PARAMETER :: ind_A_PRIAKN = 8 
  INTEGER, PARAMETER :: ind_L_PRIACC = 9 
  INTEGER, PARAMETER :: ind_L_PECACC = 10 
  INTEGER, PARAMETER :: ind_L_POAACC = 11 
  INTEGER, PARAMETER :: ind_L_FEPLUS3 = 12 
  INTEGER, PARAMETER :: ind_L_MNPLUS2 = 13 
  INTEGER, PARAMETER :: ind_L_CAPLUS2 = 14 
  INTEGER, PARAMETER :: ind_L_MGPLUS2 = 15 
  INTEGER, PARAMETER :: ind_L_KPLUS = 16 
  INTEGER, PARAMETER :: ind_L_NAPLUS = 17 
  INTEGER, PARAMETER :: ind_WD_SO2 = 18 
  INTEGER, PARAMETER :: ind_WD_HNO3 = 19 
  INTEGER, PARAMETER :: ind_WD_CO2 = 20 
  INTEGER, PARAMETER :: ind_WD_NH4OH = 21 
  INTEGER, PARAMETER :: ind_WD_H2O2 = 22 
  INTEGER, PARAMETER :: ind_WD_O3 = 23 
  INTEGER, PARAMETER :: ind_WD_HCOOH = 24 
  INTEGER, PARAMETER :: ind_WD_MHP = 25 
  INTEGER, PARAMETER :: ind_WD_PAA = 26 
  INTEGER, PARAMETER :: ind_WD_H2SO4 = 27 
  INTEGER, PARAMETER :: ind_WD_HCL = 28 
  INTEGER, PARAMETER :: ind_WD_GLY = 29 
  INTEGER, PARAMETER :: ind_WD_MGLY = 30 
  INTEGER, PARAMETER :: ind_WD_NO3MIN = 31 
  INTEGER, PARAMETER :: ind_WD_NH4PLUS = 32 
  INTEGER, PARAMETER :: ind_WD_CLMIN = 33 
  INTEGER, PARAMETER :: ind_WD_PRIACC = 34 
  INTEGER, PARAMETER :: ind_WD_FEPLUS3 = 35 
  INTEGER, PARAMETER :: ind_WD_MNPLUS2 = 36 
  INTEGER, PARAMETER :: ind_WD_PECACC = 37 
  INTEGER, PARAMETER :: ind_WD_ORGC = 38 
  INTEGER, PARAMETER :: ind_L_ORGC = 39 
  INTEGER, PARAMETER :: ind_WD_POAACC = 40 
  INTEGER, PARAMETER :: ind_WD_HPLUS = 41 
  INTEGER, PARAMETER :: ind_WD_CAPLUS2 = 42 
  INTEGER, PARAMETER :: ind_WD_MGPLUS2 = 43 
  INTEGER, PARAMETER :: ind_WD_KPLUS = 44 
  INTEGER, PARAMETER :: ind_WD_NAPLUS = 45 
  INTEGER, PARAMETER :: ind_G_HNO3 = 46 
  INTEGER, PARAMETER :: ind_G_CO2 = 47 
  INTEGER, PARAMETER :: ind_G_NH3 = 48 
  INTEGER, PARAMETER :: ind_G_H2O2 = 49 
  INTEGER, PARAMETER :: ind_G_O3 = 50 
  INTEGER, PARAMETER :: ind_G_HCOOH = 51 
  INTEGER, PARAMETER :: ind_G_MHP = 52 
  INTEGER, PARAMETER :: ind_L_GLY = 53 
  INTEGER, PARAMETER :: ind_G_GLY = 54 
  INTEGER, PARAMETER :: ind_L_MGLY = 55 
  INTEGER, PARAMETER :: ind_G_MGLY = 56 
  INTEGER, PARAMETER :: ind_G_PAA = 57 
  INTEGER, PARAMETER :: ind_G_HCL = 58 
  INTEGER, PARAMETER :: ind_G_SO2 = 59 
  INTEGER, PARAMETER :: ind_L_HNO3 = 60 
  INTEGER, PARAMETER :: ind_L_NO3MIN = 61 
  INTEGER, PARAMETER :: ind_L_NH4OH = 62 
  INTEGER, PARAMETER :: ind_L_NH4PLUS = 63 
  INTEGER, PARAMETER :: ind_L_OHMIN = 64 
  INTEGER, PARAMETER :: ind_L_HCOOH = 65 
  INTEGER, PARAMETER :: ind_L_HCOOMIN = 66 
  INTEGER, PARAMETER :: ind_L_H2SO4 = 67 
  INTEGER, PARAMETER :: ind_L_HCL = 68 
  INTEGER, PARAMETER :: ind_L_CLMIN = 69 
  INTEGER, PARAMETER :: ind_L_CO3MIN2 = 70 
  INTEGER, PARAMETER :: ind_L_H2CO3 = 71 
  INTEGER, PARAMETER :: ind_L_HCO3MIN = 72 
  INTEGER, PARAMETER :: ind_L_HSO4MIN = 73 
  INTEGER, PARAMETER :: ind_L_H2O2 = 74 
  INTEGER, PARAMETER :: ind_L_PAA = 75 
  INTEGER, PARAMETER :: ind_L_SO4MIN2 = 76 
  INTEGER, PARAMETER :: ind_L_MHP = 77 
  INTEGER, PARAMETER :: ind_L_SO3MIN2 = 78 
  INTEGER, PARAMETER :: ind_L_HPLUS = 79 
  INTEGER, PARAMETER :: ind_L_O3 = 80 
  INTEGER, PARAMETER :: ind_L_HSO3MIN = 81 
  INTEGER, PARAMETER :: ind_L_SO2 = 82 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_L_H2O = 83 
  INTEGER, PARAMETER :: ind_L_HO = 84 
  INTEGER, PARAMETER :: ind_DUMMY = 85 

! Index declaration for dummy species


! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_L_H2O = 1 
  INTEGER, PARAMETER :: indf_L_HO = 2 
  INTEGER, PARAMETER :: indf_DUMMY = 3 

! Index declaration for mechanism species

  INTEGER, PARAMETER :: spc_ind_G_SO2 = 1 
  INTEGER, PARAMETER :: spc_ind_G_HNO3 = 2 
  INTEGER, PARAMETER :: spc_ind_G_CO2 = 3 
  INTEGER, PARAMETER :: spc_ind_G_NH3 = 4 
  INTEGER, PARAMETER :: spc_ind_G_H2O2 = 5 
  INTEGER, PARAMETER :: spc_ind_G_O3 = 6 
  INTEGER, PARAMETER :: spc_ind_G_HCOOH = 7 
  INTEGER, PARAMETER :: spc_ind_G_MHP = 8 
  INTEGER, PARAMETER :: spc_ind_G_PAA = 9 
  INTEGER, PARAMETER :: spc_ind_G_HCL = 10 
  INTEGER, PARAMETER :: spc_ind_G_GLY = 11 
  INTEGER, PARAMETER :: spc_ind_G_MGLY = 12 
  INTEGER, PARAMETER :: spc_ind_A_NO3AKN = 13 
  INTEGER, PARAMETER :: spc_ind_A_NH4AKN = 14 
  INTEGER, PARAMETER :: spc_ind_A_CLAKN = 15 
  INTEGER, PARAMETER :: spc_ind_A_NAAKN = 16 
  INTEGER, PARAMETER :: spc_ind_A_SO4AKN = 17 
  INTEGER, PARAMETER :: spc_ind_A_PECAKN = 18 
  INTEGER, PARAMETER :: spc_ind_A_POAAKN = 19 
  INTEGER, PARAMETER :: spc_ind_A_PRIAKN = 20 
  INTEGER, PARAMETER :: spc_ind_L_SO2 = 21 
  INTEGER, PARAMETER :: spc_ind_L_HNO3 = 22 
  INTEGER, PARAMETER :: spc_ind_L_H2CO3 = 23 
  INTEGER, PARAMETER :: spc_ind_L_NH4OH = 24 
  INTEGER, PARAMETER :: spc_ind_L_H2O2 = 25 
  INTEGER, PARAMETER :: spc_ind_L_O3 = 26 
  INTEGER, PARAMETER :: spc_ind_L_HCOOH = 27 
  INTEGER, PARAMETER :: spc_ind_L_MHP = 28 
  INTEGER, PARAMETER :: spc_ind_L_PAA = 29 
  INTEGER, PARAMETER :: spc_ind_L_H2SO4 = 30 
  INTEGER, PARAMETER :: spc_ind_L_HCL = 31 
  INTEGER, PARAMETER :: spc_ind_L_GLY = 32 
  INTEGER, PARAMETER :: spc_ind_L_MGLY = 33 
  INTEGER, PARAMETER :: spc_ind_L_SO4MIN2 = 34 
  INTEGER, PARAMETER :: spc_ind_L_NO3MIN = 35 
  INTEGER, PARAMETER :: spc_ind_L_NH4PLUS = 36 
  INTEGER, PARAMETER :: spc_ind_L_CLMIN = 37 
  INTEGER, PARAMETER :: spc_ind_L_PRIACC = 38 
  INTEGER, PARAMETER :: spc_ind_L_PECACC = 39 
  INTEGER, PARAMETER :: spc_ind_L_ORGC = 40 
  INTEGER, PARAMETER :: spc_ind_L_POAACC = 41 
  INTEGER, PARAMETER :: spc_ind_L_HPLUS = 42 
  INTEGER, PARAMETER :: spc_ind_L_OHMIN = 43 
  INTEGER, PARAMETER :: spc_ind_L_FEPLUS3 = 44 
  INTEGER, PARAMETER :: spc_ind_L_MNPLUS2 = 45 
  INTEGER, PARAMETER :: spc_ind_L_HSO3MIN = 46 
  INTEGER, PARAMETER :: spc_ind_L_SO3MIN2 = 47 
  INTEGER, PARAMETER :: spc_ind_L_HCO3MIN = 48 
  INTEGER, PARAMETER :: spc_ind_L_CO3MIN2 = 49 
  INTEGER, PARAMETER :: spc_ind_L_HCOOMIN = 50 
  INTEGER, PARAMETER :: spc_ind_L_HSO4MIN = 51 
  INTEGER, PARAMETER :: spc_ind_L_CAPLUS2 = 52 
  INTEGER, PARAMETER :: spc_ind_L_MGPLUS2 = 53 
  INTEGER, PARAMETER :: spc_ind_L_KPLUS = 54 
  INTEGER, PARAMETER :: spc_ind_L_NAPLUS = 55 
  INTEGER, PARAMETER :: spc_ind_WD_SO2 = 56 
  INTEGER, PARAMETER :: spc_ind_WD_HNO3 = 57 
  INTEGER, PARAMETER :: spc_ind_WD_CO2 = 58 
  INTEGER, PARAMETER :: spc_ind_WD_NH4OH = 59 
  INTEGER, PARAMETER :: spc_ind_WD_H2O2 = 60 
  INTEGER, PARAMETER :: spc_ind_WD_O3 = 61 
  INTEGER, PARAMETER :: spc_ind_WD_HCOOH = 62 
  INTEGER, PARAMETER :: spc_ind_WD_MHP = 63 
  INTEGER, PARAMETER :: spc_ind_WD_PAA = 64 
  INTEGER, PARAMETER :: spc_ind_WD_H2SO4 = 65 
  INTEGER, PARAMETER :: spc_ind_WD_HCL = 66 
  INTEGER, PARAMETER :: spc_ind_WD_GLY = 67 
  INTEGER, PARAMETER :: spc_ind_WD_MGLY = 68 
  INTEGER, PARAMETER :: spc_ind_WD_NO3MIN = 69 
  INTEGER, PARAMETER :: spc_ind_WD_NH4PLUS = 70 
  INTEGER, PARAMETER :: spc_ind_WD_CLMIN = 71 
  INTEGER, PARAMETER :: spc_ind_WD_PRIACC = 72 
  INTEGER, PARAMETER :: spc_ind_WD_FEPLUS3 = 73 
  INTEGER, PARAMETER :: spc_ind_WD_MNPLUS2 = 74 
  INTEGER, PARAMETER :: spc_ind_WD_PECACC = 75 
  INTEGER, PARAMETER :: spc_ind_WD_ORGC = 76 
  INTEGER, PARAMETER :: spc_ind_WD_POAACC = 77 
  INTEGER, PARAMETER :: spc_ind_WD_HPLUS = 78 
  INTEGER, PARAMETER :: spc_ind_WD_CAPLUS2 = 79 
  INTEGER, PARAMETER :: spc_ind_WD_MGPLUS2 = 80 
  INTEGER, PARAMETER :: spc_ind_WD_KPLUS = 81 
  INTEGER, PARAMETER :: spc_ind_WD_NAPLUS = 82 
  INTEGER, PARAMETER :: spc_ind_L_H2O = 83 
  INTEGER, PARAMETER :: spc_ind_L_HO = 84 
  INTEGER, PARAMETER :: spc_ind_DUMMY = 85 

END MODULE aqchem_Parameters

