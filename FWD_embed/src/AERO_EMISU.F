
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_fwd/AERO_EMIS.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  MODULE AERO_EMIS contains emissions code required for the modal
C     aerosol module in CMAQ
C                                 Coded by Dr. Francis S. Binkowski
C                                      and Dr. Jeffrey O. Young
C
C  CONTAINS: SUBROUTINE RDEMIS_AE
C            Variable declarations needed for other subroutines in CMAQ's
C             vertical diffusion module
C
C  DEPENDENT UPON:  NONE
C
C  REVISION HISTORY:
C
C   30 Aug 01 J.Young:  dyn alloc - Use HGRD_DEFN
C   09 Oct 03 J.Gipson: added MW array for AE emis species to module contents
C   31 Jan 05 J.Young:  dyn alloc - establish both horizontal & vertical
C                       domain specifications in one module, GRID_CONF
C   26 Apr 05 P.Bhave:  removed code supporting the "old type" of emission 
C                        files that had unspeciated PM10 and PM2.5 only
C                       removed need for 'AERO_SPC.EXT' by declaring the 
C                        required variables locally
C   13 Jun 05 P.Bhave:  added vars needed for sea-salt emission processing
C                       inherit N_AE_EMIS,AE_EMIS,AE_EMIS_MAP from AE_EMIS.EXT
C                       moved RHO* parameters from RDEMIS_AE to this module
C                        for use by SSEMIS routine
C   24 Aug 07 J.Young:  Modified to enable in-line plume rise calculation for
C                       3D pt source emissions. Distinguish between PM (primary,
C                       unspeciated, file data) and AE (model speciated). Re-
C                       named RDEMIS_AE to GET_AERO_EMIS.
C   11 Apr 08 J.Kelly:  added code to emit coarse surface area
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      MODULE AERO_EMISU

      USE GRID_CONF           ! horizontal & vertical domain specifications
       
      IMPLICIT NONE
      
      INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names and map table
            
C Array dimensions
      INTEGER, PARAMETER :: NSSDIAG = 12   ! number of species in sea-salt
                                           ! diagnostic emission file
      INTEGER, PARAMETER :: NSSSPC = 4     ! number of chemical species in
                                           ! fresh sea-salt aerosol
      INTEGER, PARAMETER :: NSSMOD = 3     ! number of lognormal modes in
                                           ! sea-salt aerosol

C Molar masses (molecular wt.) and reciprocals of each aerosol species
      REAL, SAVE :: AE_EM_MW( N_AE_EMIS )
      REAL, SAVE :: R_AE_EM_MW( N_AE_EMIS )

C     Bulk component densities [ kg/m3 ]
      REAL, PARAMETER :: RHOSO4  = 1.8E3  ! density of sulfate aerosol
      REAL, PARAMETER :: RHONO3  = 1.8E3  ! density of nitrate aerosol
      REAL, PARAMETER :: RHOORG  = 2.0E3  ! density of organic aerosol
      REAL, PARAMETER :: RHOSOIL = 2.6E3  ! density of soil dust
      REAL, PARAMETER :: RHOSEAS = 2.2E3  ! density of marine aerosol
      REAL, PARAMETER :: RHOANTH = 2.2E3  ! density of elemental carbon and
                                          !  unspeciated anthropogenic aerosol
      REAL, PARAMETER :: RHOH2O  = 1.0E3  ! density of water
                                             
C Aerosol species names in sea-salt output file
      CHARACTER( 16 ), SAVE :: WRSS_SPC( NSSDIAG )  ! in sea-salt output file

C Species names in the sea-salt-emissions diagnostic file
      DATA WRSS_SPC(  1 ) / 'ANAJ'  /   ! accumulation mode sodium
      DATA WRSS_SPC(  2 ) / 'ACLJ'  /   ! accumulation mode chloride
      DATA WRSS_SPC(  3 ) / 'ASO4J' /   ! accumulation mode sulfate
      DATA WRSS_SPC(  4 ) / 'AH2OJ' /   ! accumulation mode water
      DATA WRSS_SPC(  5 ) / 'ANAK'  /   ! coarse mode sodium
      DATA WRSS_SPC(  6 ) / 'ACLK'  /   ! coarse mode chloride
      DATA WRSS_SPC(  7 ) / 'ASO4K' /   ! coarse mode sulfate 
      DATA WRSS_SPC(  8 ) / 'AH2OK' /   ! coarse mode water   
      DATA WRSS_SPC(  9 ) / 'ANUMJ' /   ! accumulation mode number
      DATA WRSS_SPC( 10 ) / 'ANUMK' /   ! coarse mode number
      DATA WRSS_SPC( 11 ) / 'ASRFJ' /   ! accumulation mode surface area
      DATA WRSS_SPC( 12 ) / 'ASRFK' /   ! coarse mode surface area

C Indices to sea-salt species in mass-emission arrays
      INTEGER, PARAMETER :: KNA  = 1   ! position of sodium
      INTEGER, PARAMETER :: KCL  = 2   ! position of chloride
      INTEGER, PARAMETER :: KSO4 = 3   ! position of sulfate
      INTEGER, PARAMETER :: KH2O = 4   ! position of water

      CONTAINS

C-----------------------------------------------------------------------

C ///////////////////////////////////////////////////////////////////////////
C  SUBROUTINE GET_AERO_EMIS reads PM emissions from gridded input file and
C  optionally from plume rise calculated point source emissions (via VDEMIS
C  array) and converts into molar-mixing-ratio units, for use in the vertical
C  diffusion routines
C
C  Key Subroutines/Functions Called:  SSEMIS
C
C  Revision History:
C
C   30 Aug 01 J.Young:  dynamic allocation - Use INTERPX
C   29 Jul 03 P.Bhave:  added compatibility with emission files that contain 
C                       PM10, PEC, POA, PNO3, PSO4, and PMF, but do not 
C                       contain PMC
C   20 Aug 03 J.Young:  return aero emissions in molar mixing ratio, ppm units
C   09 Oct 03 J.Gipson: added MW array for AE emis species to module contents
C   01 Sep 04 P.Bhave:  changed MW for primary organics from 120 to 220 g/mol,
C                       to match MWPOA in subroutine ORGAER3.
C   31 Jan 05 J.Young:  dyn alloc - removed HGRD_ID, VGRID_ID, and COORD_ID 
C                       include files because those parameters are now 
C                       inherited from the GRID_CONF module
C   26 Apr 05 P.Bhave:  removed code supporting the "old type" of emission 
C                        files that had unspeciated PM10 and PM2.5 only
C                       removed need for 'AERO_CONST.EXT' by declaring the
C                        required variables locally
C                       simplified the CONVM, CONVN, CONVS calculations
C                       updated and enhanced in-line documentation
C   03 May 05 P.Bhave:  fixed bug in the H2SO4 unit conversion, initially
C                        identified by Jinyou Liang of CARB
C   13 Jun 05 P.Bhave:  calculate sea-salt emissions; execute if MECHNAME = AE4
C                        read input fields from new OCEAN_1 file
C                        read extra input fields from MET_CRO_2D and MET_CRO_3D
C                        write diagnostic sea-salt emission file
C                        added TSTEP to call vector for diagnostic output file
C                       inherit MWs from AE_SPC.EXT instead of hardcoding
C                       find pointers to CGRID indices instead of hardcoding
C   08 Mar 07 P.Bhave&   added capability for emission files that contain 
C             S.Roselle:  POC or POA
C   30 Jan 08 P.Bhave:  added compatibility with AE5 mechanisms
C   23 Mar 08 J.Young:  modifications to allow for in-line point source emissions
C   11 Apr 08 J.Kelly:  added code to emit coarse surface area
C   09 Sep 08 P.Bhave:  backward compatibility with AE4 mechanisms
C   30 Mar 10 D.Wong:   removed code to let non-I/O processor to open CTM_SSEMIS_1
C
C  References:
C    CRC76,        "CRC Handbook of Chemistry and Physics (76th Ed)",
C                   CRC Press, 1995
C    Hobbs, P.V.   "Basic Physical Chemistry for the Atmospheric Sciences",
C                   Cambridge Univ. Press, 206 pp, 1995.
C    Snyder, J.P.  "Map Projections-A Working Manual", U.S. Geological Survey
C                   Paper 1395 U.S.GPO, Washington, DC, 1987.
c    Binkowski & Roselle  Models-3 Community Multiscale Air Quality (CMAQ)
C                   model aerosol component 1: Model Description.  
C                   J. Geophys. Res., Vol 108, No D6, 4183 
C                   doi:10.1029/2001JD001409, 2003
C-----------------------------------------------------------------------

         SUBROUTINE GET_AERO_EMISU ( JDATE, JTIME, TSTEP, EMLYRS, RJACM,
     &                              VDEMIS_AE )

         USE PT3D_DEFN, ONLY: PT3DEMIS, PMEMIS_PT ! from plume rise
         USE PMEM_DEFN                            ! defines PM emissions species
         USE PTMAP, ONLY: N_SPC_PTPM, PTPM_MAP    ! defines emissions file(s) PM mapping

         USE SUBST_MODULES                        ! stenex (for BARRIER)
!        USE SUBST_UTIL_MODULE                    ! stenex
     
         INCLUDE SUBST_RXCMMN    ! to get mech name
         INCLUDE SUBST_AE_SPC    ! aerosol species names and molecular weights
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
         INCLUDE SUBST_CONST     ! physical and mathematical constants
         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
!        INCLUDE SUBST_IODECL    ! I/O definitions and declarations
#include      SUBST_IODECL       ! I/O definitions and declarations
         INCLUDE SUBST_IOFDESC   ! file header data structure
         INCLUDE SUBST_FILES_ID  ! file name parameters

C *** Arguments

         INTEGER, INTENT( IN ) :: JDATE      ! current model date, coded YYYYDDD
         INTEGER, INTENT( IN ) :: JTIME      ! current model time, coded HHMMSS
         INTEGER, INTENT( IN ) :: TSTEP( 2 ) ! time step vector (HHMMSS)
                                             ! TSTEP(1) = local output step
                                             ! TSTEP(2) = sciproc sync. step (chem)
         INTEGER, INTENT( IN ) :: EMLYRS   ! number of vertical layers on emis file
         REAL,    INTENT( IN ) :: RJACM( NCOLS,NROWS,NLAYS ) ! recip of mid-layer
                                                             ! Jacobian [1/m]
         REAL,    INTENT( OUT )   :: VDEMIS_AE( :,:,:,: ) ! aerosol emissions
                                                       ! [ppmv/s] for mass & number spcs
                                                       ! [m2/mol/s] for surface area spcs

C *** Local Variables

C     Geometric Constants
         REAL, PARAMETER :: F6DPI = 6.0 / PI
         REAL, PARAMETER :: F6DPIM9 = 1.0E-9 * F6DPI

C     Aerosol version name
         CHARACTER( 16 ), SAVE :: AE_VRSN

C     Indices in the VDEMIS_AE array
         INTEGER, SAVE :: VSO4AJ    ! accumulation mode sulfate
         INTEGER, SAVE :: VSO4AI    ! Aitken mode sulfate
         INTEGER, SAVE :: VNO3AJ    ! accumulation mode nitrate
         INTEGER, SAVE :: VNO3AI    ! Aitken mode nitrate
         INTEGER, SAVE :: VORGPAJ   ! accumulation mode primary anthropogenic
                                    ! organic aerosol
         INTEGER, SAVE :: VORGPAI   ! Aitken mode primary anthropogenic
                                    ! organic aerosol
         INTEGER, SAVE :: VECJ      ! accumulation mode elemental carbon
         INTEGER, SAVE :: VECI      ! Aitken mode elemental carbon
         INTEGER, SAVE :: VP25AJ    ! accumulation mode unspeciated PM2.5
         INTEGER, SAVE :: VP25AI    ! Aitken mode unspeciated PM2.5
         INTEGER, SAVE :: VANTHA    ! coarse mode anthropogenic aerosol
         INTEGER, SAVE :: VSEAS     ! coarse mode marine aerosol
         INTEGER, SAVE :: VSOILA    ! coarse mode soil-derived aerosol
         INTEGER, SAVE :: VAT0      ! Aitken mode number
         INTEGER, SAVE :: VAC0      ! accumulation mode number
         INTEGER, SAVE :: VCOR0     ! coarse mode number
         INTEGER, SAVE :: VSURFAT   ! Aitken mode surface area
         INTEGER, SAVE :: VSURFAC   ! accumulation mode surface area
         INTEGER, SAVE :: VSURFCO   ! coarse mode surface area
         INTEGER, SAVE :: VH2OAJ    ! accumulation mode water 
         INTEGER, SAVE :: VH2OAI    ! Aitken mode water
         INTEGER, SAVE :: VNAJ      ! accumulation mode sodium
         INTEGER, SAVE :: VNAI      ! Aitken mode sodium
         INTEGER, SAVE :: VCLJ      ! accumulation mode chloride
         INTEGER, SAVE :: VCLI      ! Aitken mode chloride
         INTEGER, SAVE :: VNAK      ! coarse mode sodium
         INTEGER, SAVE :: VCLK      ! coarse mode chloride
         INTEGER, SAVE :: VSO4K     ! coarse mode sulfate
         INTEGER, SAVE :: VH2OK     ! coarse mode water

C     Geometric mean diameter by volume (or mass) of emitted particles in 
C     each mode [ m ].  See paragraph #14 of Binkowski & Roselle (2003)
         REAL, PARAMETER :: DGVEM_AT = 0.03E-6 ! Aitken mode
         REAL, PARAMETER :: DGVEM_AC = 0.3E-6  ! accumulation mode
         REAL, PARAMETER :: DGVEM_CO = 6.0e-6  ! coarse mode

C     Geometric standard deviation of emitted particles in each mode, as
C     described in paragraph #14 of Binkowski & Roselle (2003)
         REAL, PARAMETER :: SGEM_AT = 1.7      ! Aitken mode       
         REAL, PARAMETER :: SGEM_AC = 2.0      ! accumulation mode 
         REAL, PARAMETER :: SGEM_CO = 2.2      ! coarse mode       

C     Variables for converting mass emissions rate to number emissions rate
         REAL, SAVE :: FACTNUMAT               ! Aitken mode
         REAL, SAVE :: FACTNUMAC               ! accumulation mode
         REAL, SAVE :: FACTNUMC                ! Coarse mode

C     Variables for converting mass emissions rate to 2nd moment emissions rate
         REAL, SAVE :: FACTM2AT                ! Aitken mode
         REAL, SAVE :: FACTM2AC                ! acumulation mode
         REAL, SAVE :: FACTM2CO                ! coarse mode

C     Variables for calculating the volume of each grid cell
         REAL  DX1, DX2                          ! grid-cell width and length [m]
         REAL  GRDAREA                           ! grid area [m2]
         REAL, ALLOCATABLE, SAVE :: GRDHGT( : )  ! grid height [sigma]
         REAL, ALLOCATABLE, SAVE :: GRDVOL( : )  ! grid volume [m2*sigma]

C     Emission rate of all aerosol species interpolated to current time
         REAL, ALLOCATABLE, SAVE :: EMBUFF( :,:,:,: )  ! in all grid cells
         REAL, ALLOCATABLE, SAVE :: PM_EM( : )         ! in one grid cell
         
C     Factor for converting aerosol emissions from input units ...
         REAL CONVEM_PM_MASS                       ! into [ug/sec]
         REAL, ALLOCATABLE, SAVE :: CONVEM_PM( : ) ! into [ug/m2/sec]
         REAL GSFAC                                ! into [ug/m3/sec]

C     Variables interpolated from the meteorological input files
         REAL PRES( NCOLS,NROWS )              ! atmospheric pressure [Pa]
         REAL TA  ( NCOLS,NROWS )              ! air temperature [K]
         REAL QV  ( NCOLS,NROWS )              ! H2O mass mixing ratio [kg/kg]
         REAL DENS( NCOLS,NROWS,NLAYS )        ! air density [kg/m3]
         REAL WSPD10( NCOLS,NROWS )            ! wind speed at 10m [m/s]

C     Variables for converting emission rates into molar-mixing-ratio units
         REAL, PARAMETER :: GPKG = 1.0E+03     ! g/kg
         REAL, PARAMETER :: MGPG = 1.0E+06     ! ug/g
         REAL, PARAMETER :: RAVO = 1.0 / AVO   ! reciprocal of Avogadro
!        REAL, PARAMETER :: MWAIR = 28.9628    ! molar mass of dry air [g/mol]
                                               ! assuming 78.06% N2, 21% O2, 
                                               ! and 0.943% Ar 
                                               ! (Source : Hobbs, 1995) pp 69-70
         REAL  CONVM     ! conversion factor for mass emissions [m3/mol]
         REAL  CONVN     ! conversion factor for number emissions [1e6*m3]
         REAL  CONVS     ! conversion factor for surface area emissions [m3/mol]

C     Domain decomposition info from emission and meteorology files
         INTEGER  GXOFF, GYOFF          ! origin offset
         INTEGER, SAVE :: STARTCOL, ENDCOL, STARTROW, ENDROW
         INTEGER, SAVE :: STRTCOL_O1, ENDCOL_O1, STRTROW_O1, ENDROW_O1
         INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C     Aerosol mass-emission rates [ug/m3/s]
         REAL EPOC          ! fine primary organic carbon
         REAL EPEC          ! fine primary elemental carbon
         REAL EPNO3         ! fine primary nitrate
         REAL EPSO4         ! fine primary sulfate
         REAL EPMF          ! fine primary unspeciated pm
         REAL EPMCO         ! coarse anthropogenic aerosol
         REAL ESOILCO       ! coarse soil-derived aerosol
         REAL ESEASCO       ! coarse marine aerosol

C     Variables for handling vapor-phase sulfuric acid emissions
         REAL EMSULF                          ! emission rate [ppmv/s]
         REAL, PARAMETER :: MWH2SO4 = 98.0    ! molar mass [g/mol]
         INTEGER, SAVE :: VSULF               ! index to H2SO4 in VDEMIS array

C     Speciation factors for coarse mode emissions, from paragraph #15 of
C     Binkowski & Roselle (2003)
         REAL, PARAMETER :: FAC_DUST  = 0.90  ! (fugitive dust)/PMC
         REAL, PARAMETER :: FAC_OTHER = 0.10  ! (non-fugitive dust)/PMC

C     Variables for handling sea-salt emissions
         REAL  OCEAN( NCOLS,NROWS )           ! fractional seawater cover
         REAL  SZONE( NCOLS,NROWS )           ! fractional surf-zone cover
         REAL  SSOUTD( NSSDIAG )              ! all emission rates for 
                                              !  diagnostic output file
         REAL  SSOUTM( NSSSPC,NSSMOD )        ! mass emission rates [ug/m3/s]
         REAL  SSOUTN( NSSMOD )               ! number emission rates [1/m3/s]
         REAL  SSOUTS( NSSMOD )               ! surface-area emisrates [m2/m3/s]

C     Grid-specific values for sea-salt calculations
         REAL  OFRAC                          ! fractional seawater cover
         REAL  SFRAC                          ! fractional surf-zone cover
         REAL  BLKPRS                         ! atmospheric pressure [Pa]
         REAL  BLKTA                          ! air temperature [K]
         REAL  BLKQV                          ! H2O mass mixing ratio [kg/kg]
         REAL  U10                            ! wind speed at 10m [m/s]
         REAL  RLAY1HT                        ! reciprocal of layer-1 hgt [1/m]
         REAL  AIRVOL                         ! grid-cell volume [m3]

C     Factors for splitting primary carbon emissions into Aitken and
C     accumulation modes, from paragraph #12 of Binkowski & Roselle (2003)
         REAL, PARAMETER :: FACEM25_ACC  = 0.999 ! accumulation mode
         REAL, PARAMETER :: FACEM25_ATKN = 0.001 ! Aitken mode

C     Mode-specific mass-emission rates [ug/m3/s]
         REAL EPM25AT       ! Aitken mode unspeciated aerosol
         REAL EPM25AC       ! accumulation mode unspeciated aerosol
         REAL EPORGAT       ! Aitken mode primary organic aerosol
         REAL EPORGAC       ! accumulation mode primary organic aerosol
         REAL EPECAT        ! Aitken mode elemental carbon
         REAL EPECAC        ! accumulation mode elemental carbon
         REAL EPSO4AT       ! Aitken mode primary sulfate
         REAL EPSO4AC       ! accumulation mode primary sulfate
         REAL EPSO4CO       ! coarse mode primary sulfate
         REAL EPNO3AT       ! Aitken mode primary nitrate
         REAL EPNO3AC       ! accumulation mode primary nitrate
         REAL EPNAAT        ! Aitken mode sodium
         REAL EPNAAC        ! accumulation mode sodium
         REAL EPNACO        ! coarse mode sodium
         REAL EPCLAT        ! Aitken mode chloride
         REAL EPCLAC        ! accumulation mode chloride
         REAL EPCLCO        ! coarse mode chloride
         REAL EPH2OAT       ! Aitken mode primary water
         REAL EPH2OAC       ! accumulation mode primary water
         REAL EPH2OCO       ! coarse mode primary water

C     Factors for converting aerosol mass concentration [ug/m3] to 3rd
C     moment concentration [m3/m3]
         REAL, PARAMETER :: SO4FAC  = F6DPIM9 / RHOSO4
         REAL, PARAMETER :: NO3FAC  = F6DPIM9 / RHONO3
         REAL, PARAMETER :: ORGFAC  = F6DPIM9 / RHOORG
         REAL, PARAMETER :: SOILFAC = F6DPIM9 / RHOSOIL
         REAL, PARAMETER :: SEASFAC = F6DPIM9 / RHOSEAS
         REAL, PARAMETER :: ANTHFAC = F6DPIM9 / RHOANTH

C     Third moment emissions rates [m3/m3/s]
         REAL EMISM3AT      ! Aitken mode
         REAL EMISM3AC      ! accumulation mode
         REAL EMISM3CO      ! coarse mode

C     Number emissions rates [1/m3/s]
         REAL EM_NUMATKN    ! Aitken mode
         REAL EM_NUMACC     ! accumulation mode
         REAL EM_NUMCOR     ! coarse mode

C     Surface area emission rates [m2/m3/s]
         REAL EM_SRFATKN    ! Aitken mode
         REAL EM_SRFACC     ! accumulation mode
         REAL EM_SRFCOR     ! coarse mode

C     Variables for writing out sea-salt emission rates
         INTEGER, SAVE     :: WSTEP  = 0    ! local write counter
         INTEGER, SAVE     :: SDATE, STIME  ! scenario start date & time
         INTEGER           :: MDATE, MTIME  ! SSEMIS write date&time
         INTEGER, EXTERNAL :: ENVINT        ! fxn not declared in IODECL3.EXT
         INTEGER, EXTERNAL :: TIME2SEC      ! fxn not declared in IODECL3.EXT
         INTEGER, EXTERNAL :: SEC2TIME      ! fxn not declared in IODECL3.EXT
         REAL, ALLOCATABLE, SAVE :: SSBF( :,:,: ) ! seasalt emission accumulator
         REAL              :: WRSS( NCOLS,NROWS ) ! seasalt emission write buffer
         LOGICAL, SAVE     :: SSEMDIAG      ! flag for creating SSEMIS output file
         CHARACTER( 16 ), SAVE :: CTM_SSEMDIAG = 'CTM_SSEMDIAG' ! environment var
                                            ! for SSEMDIAG file
         LOGICAL, EXTERNAL :: ENVYN         ! get environment variable as boolean

C     Miscellaneous variables
         LOGICAL, SAVE     :: FIRSTIME = .TRUE.
         INTEGER, SAVE     :: LOGDEV
         INTEGER, EXTERNAL :: INDEX1, TRIMLEN
         INTEGER INDX
         INTEGER STATUS                   ! ENV..., ALLOCATE status
         CHARACTER( 80 ) :: VARDESC       ! environment variable description
         CHARACTER( 96 ) :: XMSG = ' '
         CHARACTER( 16 ), SAVE :: PNAME = 'GET_AERO_EMIS'
         CHARACTER( 16 ) :: VNAME       ! temp var for species names
         INTEGER C, R, L, N, V              ! Loop indices

         REAL AE_TEMP        !Temporary variable

C ----------------------------------------------------------------------

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            LOGDEV = INIT3()

            IF ( INDEX ( MECHNAME, 'AE3' ) .GT. 0 ) THEN
               AE_VRSN  = 'AE3'
            ELSE IF ( INDEX ( MECHNAME, 'AE4' ) .GT. 0 ) THEN
               AE_VRSN  = 'AE4'
            ELSE IF ( INDEX ( MECHNAME, 'AE5' ) .GT. 0 ) THEN
               AE_VRSN  = 'AE5'
            ELSE
               XMSG = 'This version of the emission processing code '
     &             // 'can only be used with the AE3, AE4, or AE5 '
     &             // 'aerosol mechanisms.'
                     CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF ! check on MECHNAME

C *** Set indices for the VDEMIS_AE array using the AE_EMIS table
C     This is the subset of the AE_SPC table that have emissions

            VNAME = 'ASO4J'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VSO4AJ = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'ASO4I'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VSO4AI = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'ANO3J'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VNO3AJ = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'ANO3I'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VNO3AI = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'AORGPAJ'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VORGPAJ = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'AORGPAI'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VORGPAI = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'AECJ'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VECJ = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'AECI'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VECI = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'A25J'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VP25AJ = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'A25I'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VP25AI = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'ACORS'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VANTHA = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            IF ( AE_VRSN .EQ. 'AE3' ) THEN
               VNAME = 'ASEAS'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VSEAS = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
            END IF

            VNAME = 'ASOIL'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VSOILA = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'NUMATKN'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VAT0 = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'NUMACC'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VAC0 = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'NUMCOR'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VCOR0 = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'SRFATKN'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VSURFAT = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            VNAME = 'SRFACC'
            N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
            IF ( N .NE. 0 ) THEN
               VSURFAC = N
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

            IF ( AE_VRSN .NE. 'AE3' ) THEN

               VNAME = 'AH2OJ'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VH2OAJ = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'AH2OI'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VH2OAI = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ANAJ'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VNAJ = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ANAI'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VNAI = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ACLJ'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VCLJ = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ACLI'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VCLI = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ANAK'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VNAK = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ACLK'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VCLK = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'ASO4K'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VSO4K = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF

               VNAME = 'AH2OK'
               N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
               IF ( N .NE. 0 ) THEN
                  VH2OK = N
               ELSE
                  XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
               END IF
               
               IF ( AE_VRSN .EQ. 'AE5' ) THEN
                  VNAME = 'SRFCOR'
                  N = INDEX1( VNAME, N_AE_EMIS, AE_EMIS )
                  IF ( N .NE. 0 ) THEN
                     VSURFCO = N
                  ELSE
                     XMSG = 'Could not find ' // VNAME // 'in AE_EMIS table'
                     CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
                  END IF
               END IF               

C *** set up for sea-salt emission processing

C *** Open the ocean file, which contains the ocean and surf-zone fractions

               IF ( .NOT. OPEN3( OCEAN_1, FSREAD3, PNAME ) ) THEN
                  XMSG = 'Could not open '// OCEAN_1 // ' file'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
               
               IF ( .NOT. DESC3( OCEAN_1 ) ) THEN
                  XMSG = 'Could not get '// OCEAN_1 // ' file description'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               END IF

               CALL SUBHFILE ( OCEAN_1, GXOFF, GYOFF,
     &                         STRTCOL_O1, ENDCOL_O1, STRTROW_O1, ENDROW_O1 )

C *** Get sea-salt-emission diagnostic output file flag.

               SSEMDIAG = .FALSE.         ! default
               VARDESC = 'Flag for writing the sea-salt-emission diagnostic file'
               SSEMDIAG = ENVYN( CTM_SSEMDIAG, VARDESC, SSEMDIAG, STATUS )
               IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
               IF ( STATUS .EQ. 1 ) THEN
                  XMSG = 'Environment variable improperly formatted'
                  CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               ELSE IF ( STATUS .EQ. -1 ) THEN
                  XMSG =
     &                'Environment variable set, but empty ... Using default:'
                  WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
               ELSE IF ( STATUS .EQ. -2 ) THEN
                  XMSG = 'Environment variable not set ... Using default:'
                  WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
               END IF

               IF ( SSEMDIAG ) THEN    ! Open the sea-salt emission dignostic file

                  SDATE = ENVINT( 'CTM_STDATE', 'Scenario Start (YYYYJJJ)', 0, STATUS )
                  STIME = ENVINT( 'CTM_STTIME', 'Scenario Start (HHMMSS)', 0, STATUS )

                  IF ( MYPE .EQ. 0 ) THEN
                     CALL OPSSEMIS ( SDATE, STIME, TSTEP )
!                    IF ( .NOT. CLOSE3( CTM_SSEMIS_1 ) ) THEN
!                       XMSG = 'Could not close '// CTM_SSEMIS_1 // ' file'
!                       CALL M3EXIT( PNAME, SDATE, STIME, XMSG, XSTAT1 )
!                    END IF
                  END IF

!                 CALL SUBST_BARRIER

C Make sure all processors open file

!                 IF ( .NOT. OPEN3( CTM_SSEMIS_1, FSRDWR3, PNAME ) ) THEN
!                    XMSG = 'Could not open '// CTM_SSEMIS_1 // ' file'
!                    CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!                 END IF

                  ALLOCATE ( SSBF( NSSDIAG,MY_NCOLS,MY_NROWS ), STAT = STATUS )
                  IF ( STATUS .NE. 0 ) THEN
                     XMSG = '*** SSBF memory allocation failed'
                     CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                  END IF

               END IF
               
            END IF   ! AE_VRSN .NE. 'AE3'

C *** Save array of MWs and reciprocals [AE_EM_MW( V ) required for pa_update.F]

            DO V = 1, N_AE_EMIS
               AE_EM_MW( V ) = AE_MOLWT( AE_EMIS_MAP( V ) )
               R_AE_EM_MW( V ) = 1.0 / AE_EM_MW( V )
            END DO

C *** Calculate factors for converting 3rd moment emission rates into number
C     emission rates.  See Equation 7b of Binkowski & Roselle (2003)

            FACTNUMAT = EXP( 4.5 * LOG( SGEM_AT ) ** 2 ) / DGVEM_AT ** 3
            FACTNUMAC = EXP( 4.5 * LOG( SGEM_AC ) ** 2 ) / DGVEM_AC ** 3
            FACTNUMC  = EXP( 4.5 * LOG( SGEM_CO ) ** 2 ) / DGVEM_CO ** 3

C *** Calculate factors for converting 3rd moment emission rates into 2nd 
C     moment emission rates.  See Equation 7c of Binkowski & Roselle (2003)

            FACTM2AT  = EXP( 0.5 * LOG( SGEM_AT ) ** 2 ) / DGVEM_AT
            FACTM2AC  = EXP( 0.5 * LOG( SGEM_AC ) ** 2 ) / DGVEM_AC
            FACTM2CO  = EXP( 0.5 * LOG( SGEM_CO ) ** 2 ) / DGVEM_CO

C *** Allocate memory for PM_EM, EMBUFF, GRDHGT, GRDVOL, CONVEM_PM, and SSBF

            ALLOCATE ( PM_EM( N_EMIS_PM ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** PM_EM memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( EMBUFF( MY_NCOLS,MY_NROWS,EMLYRS,N_EMIS_PM ),
     &                 STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** EMBUFF memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( CONVEM_PM( EMLYRS ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** CONVEM_PM memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( GRDHGT( EMLYRS ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** GRDHGT memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            ALLOCATE ( GRDVOL( EMLYRS ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** GRDVOL memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
C *** Calculate length and width of each grid cell
C     note: crude estimate is made for LAT/LONG coordinate systems

            IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
               DX1 = DG2M * XCELL_GD ! in m
               DX2 = DG2M * YCELL_GD
     &             * COS( PI180*( YORIG_GD + YCELL_GD
     &             * FLOAT( GL_NROWS/2 ) ) ) ! in m
            ELSE
               DX1 = XCELL_GD        ! in m
               DX2 = YCELL_GD        ! in m
            END IF

C *** Calculate height of grid cell in each layer in sigma coordinates
C     Multiply by grid area [m2] to obtain grid volume

            GRDAREA = DX1 * DX2
            DO L = 1, EMLYRS
               GRDHGT( L ) = X3FACE_GD( L ) - X3FACE_GD( L-1 )
               GRDVOL( L ) = GRDHGT( L ) * GRDAREA
            END DO

C *** Calculate scaling factor for converting aerosol emissions from
C     their input units to [ug/s] and then to [ug/m2/s] using layer-
C     specific grid-cell volume

            IF ( UNITSCKPM .EQ. 'G/S' .OR.
     &           UNITSCKPM .EQ. 'g/s' ) THEN
               CONVEM_PM_MASS = MGPG                  ! (g/s) -> (ug/s)
            ELSE IF ( UNITSCKPM .EQ. 'KG/HR' .OR.
     &                UNITSCKPM .EQ. 'kg/hr' ) THEN
               CONVEM_PM_MASS = GPKG * MGPG / 3600.0  ! (kg/hr) -> (ug/s)
            ELSE
               XMSG = 'Units incorrect on ' // EMIS_1
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
            END IF

            DO L = 1, EMLYRS
               CONVEM_PM( L ) = CONVEM_PM_MASS / GRDVOL( L )
            END DO

C *** Find location of sulfuric acid vapor in VDEMIS array

            VNAME = 'SULF'
            INDX = INDEX1 ( VNAME, N_GC_EMIS, GC_EMIS )
            IF ( INDX .NE. 0 ) THEN
               VSULF = INDX  ! index for vapor-phase H2SO4 emissions
            ELSE
               XMSG = 'Could not find ' // VNAME // 'in gas table'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
            END IF

C *** Get domain decomposition info from input files

            CALL SUBHFILE ( EMIS_1, GXOFF, GYOFF,
     &                      STARTCOL, ENDCOL, STARTROW, ENDROW )

            CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                      STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         END IF    ! FIRSTIME

C ----------------------------------------------------------------------

C *** Read air density [kg/m3], atmospheric pressure [Pa], air temperature 
C     [K], specific humidity [kg H2O / kg air], and 10m wind speed [m/s] 
C     from meteorology file.  Interpolate to the current time.  Store results 
C     in DENS, PRES, TA, QV, and WSPD10 arrays.

         IF ( .NOT. INTERPX( MET_CRO_3D, 'DENS', PNAME,
     &                       STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                       1,NLAYS, JDATE, JTIME, DENS ) ) THEN
            XMSG = 'Could not interpolate DENS from ' // MET_CRO_3D
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( AE_VRSN .NE. 'AE3' ) THEN
            IF ( .NOT. INTERPX( MET_CRO_3D, 'PRES', PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                          1,1, JDATE, JTIME, PRES ) ) THEN
               XMSG = 'Could not interpolate PRES from ' // MET_CRO_3D
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            IF ( .NOT. INTERPX( MET_CRO_3D, 'TA', PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                          1,1, JDATE, JTIME, TA ) ) THEN
               XMSG = 'Could not interpolate TA from ' // MET_CRO_3D
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
            
            IF ( .NOT. INTERPX( MET_CRO_3D, 'QV', PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                          1,1, JDATE, JTIME, QV ) ) THEN
               XMSG = 'Could not interpolate QV from ' // MET_CRO_3D
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            IF ( .NOT. INTERPX( MET_CRO_2D, 'WSPD10', PNAME,
     &                          STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                          1, 1, JDATE, JTIME, WSPD10 ) ) THEN
               XMSG = 'Could not find WSPD10 in ' // MET_CRO_2D //
     &                ' search for WIND10.'
               CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
               IF ( .NOT. INTERPX( MET_CRO_2D, 'WIND10', PNAME,
     &                             STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3,
     &                             1, 1, JDATE, JTIME, WSPD10 ) ) THEN
                  XMSG = 'Could not find WIND10 or WSPD10 in ' // MET_CRO_2D
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
            END IF


         END IF   ! AE_VRSN .NE. 'AE3'

C *** Loop over each grid cell

         DO L = 1, EMLYRS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS

C *** Calculate scaling factor for converting mass emissions into [ug/m3/s]
C     note: RJACM converts grid heights from sigma coordinates to meters
C     Also calculate scaling factors for converting to molar-mixing-ratio units

                  GSFAC = CONVEM_PM( L ) * RJACM( C,R,L )
                  CONVM = MWAIR / GPKG / DENS( C,R,L )   ! [m3/mol]
                  CONVN = CONVM * RAVO * 1.0E+06         ! 10^6 * [m3]
                  CONVS = CONVM                          ! [m3/mol]

                  AE_TEMP = GSFAC * 1.E+6 / 86400.

                  VDEMIS_AE( VSO4AI, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VSO4AI )
                  VDEMIS_AE( VSO4AJ, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VSO4AJ )
                  VDEMIS_AE( VNO3AI, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VNO3AI )
                  VDEMIS_AE( VNO3AJ, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VNO3AJ )
                  VDEMIS_AE( VORGPAI,L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VORGPAI )
                  VDEMIS_AE( VORGPAJ,L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VORGPAJ )
                  VDEMIS_AE( VECI,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VECI )
                  VDEMIS_AE( VECJ,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VECJ )
                  VDEMIS_AE( VP25AI, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VP25AI )
                  VDEMIS_AE( VP25AJ, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VP25AJ )
                  VDEMIS_AE( VANTHA, L,C,R ) = AE_TEMP   * CONVM * R_AE_EM_MW( VANTHA )
                  VDEMIS_AE( VSOILA, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VSOILA )
                  VDEMIS_AE( VAT0,   L,C,R ) = 0.
                  VDEMIS_AE( VAC0,   L,C,R ) = 0.
                  VDEMIS_AE( VCOR0,  L,C,R ) = 0.
                  VDEMIS_AE( VSURFAT,L,C,R ) = 0.
                  VDEMIS_AE( VSURFAC,L,C,R ) = 0.

                  IF ( AE_VRSN .EQ. 'AE3' ) THEN
                     VDEMIS_AE( VSEAS,  L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VSEAS )
                  ELSE
                     VDEMIS_AE( VH2OAI, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VH2OAI )
                     VDEMIS_AE( VH2OAJ, L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VH2OAJ )
                     VDEMIS_AE( VNAI,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VNAI )
                     VDEMIS_AE( VNAJ,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VNAJ )
                     VDEMIS_AE( VCLI,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VCLI )
                     VDEMIS_AE( VCLJ,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VCLJ )
                     VDEMIS_AE( VNAK,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VNAK )
                     VDEMIS_AE( VCLK,   L,C,R ) = AE_TEMP  * CONVM * R_AE_EM_MW( VCLK )
                     VDEMIS_AE( VSO4K,  L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VSO4K )
                     VDEMIS_AE( VH2OK,  L,C,R ) = AE_TEMP * CONVM * R_AE_EM_MW( VH2OK )
                  END IF   ! AE_VRSN .EQ. 'AE3'

                  IF ( AE_VRSN .EQ. 'AE5' ) THEN
                     VDEMIS_AE( VSURFCO,L,C,R ) = 0.
                  END IF

               END DO   ! loop on MY_NCOLS
            END DO   ! loop on MY_NROWS
         END DO   ! loop on EMLYRS


         RETURN

         END SUBROUTINE GET_AERO_EMISU

      END MODULE AERO_EMISU
      
