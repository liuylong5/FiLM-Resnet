
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/aero/aero5_adj_fwd/aero_driver.F,v 1.4 2013/05/14 21:06:51 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C >>> 08/04/2000 Changes necessary to be able to read and process
C two different types of emissions files.
C the first type is the existing opperational PM2.5 & PM10 unspeciated
C file. The new file format has speciated emissions. 
C >>> This version uses the FORTRAN 90 feature for runtime memory
C allocation.

C 1/12/99 David Wong at LM: 
C   -- introduce new variable MY_NUMBLKS (eliminate NUMBLKS)
C   -- re-calculate NOXYZ accordingly
C FSB Updated for inclusion of surface area / second moment
C 25 Sep 00 (yoj) various bug fixes, cleanup to coding standards
C   Jeff - Dec 00 - move CGRID_MAP into f90 module
C FSB/Jeff - May 01 - optional emissions processing
C   Jerry Gipson - Jun 01 - added SOA linkages for saprc99
C   Bill Hutzell - Jun 01 - simplified CBLK mapping
C   Jerry Gipson - Jun 03 - modified for new soa treatment
C   Jerry Gipson - Aug 03 - removed SOA prod form alkenes & added 
C       emission adjustment factors for ALK & TOL ( RADM2 & SAPRC99 only)
C   Shawn Roselle - Jan 04
C   - removed SOA from transported aerosol surface area
C   - fixed bug in calculation of wet parameters.  Previously, DRY aerosol
C      parameters were being written to the AERDIAG files and mislabeled
C      as WET.
C   Prakash Bhave - May 04
C   - changed AERODIAG species (added RH; removed M0 & M2dry)
C   Jeff Young - Jan 05 - dyn alloc 
C   - establish both horizontal & vertical domain specifications in one module
c   Uma Shankar and Prakash Bhave - Jun 05
c   - added code to handle the following species: ANAI, ANAJ, ANAK, ACLI, 
c     ACLJ, ACLK, ASO4K, AH2OK, ANO3K, and HCL; removed code for ASEAS
c   - removed obsolete MW variables
C   Prakash Bhave - Jul 05 - added PM25 mass-fraction calculations
C   Jeff Young - Feb 06 - trap fractional humidity above 0.005
C   Prakash Bhave - Apr 06 - added GAMMA_N2O5 to the AEROPROC call vector
C       and the aerosol diagnostic file
C   Prakash Bhave - May 06 - changed units of DG variables from m to um in
C       the aerosol diagnostic file as suggested by Dr. Bill Hutzell
C   Sergey Napelenok - Sep 07 - SOA updates
C   - added code to handle the following species: AALKJ, ATOL1J, ATOL2J, 
C     ATOL3J, AXYL1J, AXYL2J, AXYL3J, ABNZ1J, ABNZ2J, ABNZ3J, AISO1J, AISO2J, 
C     AISO3J, ATRP1J, ATRP2J, ASQTJ, AORGCJ, TOLNRXN, TOLHRXN, XYLNRXN, 
C     XYLHRXN, BNZNRXN, BNZHRXN, ISOPRXN, and SESQRXN
C   - removed following species: AORGAI, AORGAJ, AORGBI, AORGBJ, OLIRXN, 
C     CSLRXN, TOLRXN, XYLRXN
C   Prakash Bhave - Oct 07 - SOA updates
C   - added semi-volatile vapors to the CBLK array; moved ppm -> ug/m3 unit
C     conversion from the ORGAER subroutine to this program
C   - updated definition of DRY aerosol to include nonvolatile SOA species
C   - removed adjustment factors for TOLAER (SPTOL, RDTOL) because benzene is
C     now an explicit species so all of the reacted TOL can produce SOA
C   - removed code to handle TERPSP (obsolete); renamed TERPRXN as TRPRXN
C   David Wong - Jan 08 - rearranged calculation of dry 3rd moments to avoid
C      NaN on some compilers (using the M3SUBT variable)
C   Prakash Bhave - Jan 08 - updated MECHNAME check from AE4 to AE5
C   Golam Sarwar -  Mar 08 - added a heterogeneous reaction producing HONO 
C   J.Kelly - Apr 08 - added code to account for new species (ANH4K &
C                      SRFCOR) and variable coarse std. deviation
C                    - removed MW coding now located in AERO_INFO.f
C                    - added LIMIT_Sg flag for call to GETPAR
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AERO ( CGRID, JDATE, JTIME, TSTEP )

      USE GRID_CONF           ! horizontal & vertical domain specifications
      USE CGRID_SPCS          ! CGRID species number and offsets

      USE aero_data           ! aero variable data
      USE precursor_data      ! gas phase aero precursor data
      USE soa_defn            ! gas soa data
      USE met_data            ! air properties data, eg., T, P, RH

      USE CHK_FILES_DEFS      ! adjoint checkpoint files definition
      USE CHK_FILES           ! adjoint checkpoint files routines

      IMPLICIT NONE
 
C *** includes:
 
      INCLUDE SUBST_GC_SPC    ! gas chemistry species table
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map
                              ! table
      INCLUDE SUBST_AE_SPC    ! aerosol species table
                              ! table
      INCLUDE SUBST_NR_SPC    ! non-reactive species table
      INCLUDE SUBST_GC_G2AE   ! gas chem aerosol species and map table
      INCLUDE SUBST_NR_N2AE   ! non-react aerosol species and map table
      INCLUDE SUBST_RXCMMN    ! to get mech name
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
!     INCLUDE SUBST_IODECL    ! I/O definitions and declarations
#include      SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID  ! file name parameters
                              ! (req IOPARMS)
      INCLUDE SUBST_PACTL_ID  ! process analysis

C *** arguments:

C *** CGRID is conc field (including gas and aerosol variables)
      REAL, POINTER :: CGRID( :,:,:,: )              !  concentrations
      INTEGER      JDATE        ! Current model date , coded YYYYDDD
      INTEGER      JTIME        ! Current model time , coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

      INTEGER, SAVE :: LOGDEV             ! unit number for the log file

C *** local variables:
      CHARACTER( 16 ), SAVE :: PNAME = 'AERO_DRIVER'
      CHARACTER( 16 ) :: VNAME            ! variable name
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER   MDATE, MTIME, MSTEP  ! julian date, time and 
                                     ! timestep in sec
      INTEGER   C, R, L, V, N        ! loop counters
      INTEGER   SPC, SIZE            ! species loop counter
      INTEGER   STRT, FINI           ! loop induction variables
      INTEGER   ALLOCSTAT            ! memory allocation status

      LOGICAL   LERROR               ! Error flag

C *** External Functions not previously declared in IODECL3.EXT:
      INTEGER, EXTERNAL :: SECSDIFF, SEC2TIME, TIME2SEC, INDEX1
      LOGICAL, EXTERNAL :: ENVYN     ! get environment variable as boolean

C *** Grid description
      REAL DX1                 ! Cell x-dimension
      REAL DX2                 ! Cell y-dimension
      INTEGER GXOFF, GYOFF     ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

C *** Variable to set time step for writing visibility file
      INTEGER, SAVE :: WSTEP  = 0   ! local write counter
      LOGICAL, SAVE :: WRITETIME = .FALSE. ! local write flag

C *** meteorological variables
      REAL PRES   ( NCOLS,NROWS,NLAYS )  ! Atmospheric pressure [ Pa ]
      REAL TA     ( NCOLS,NROWS,NLAYS )  ! Air temperature [ K ] 
      REAL DENS   ( NCOLS,NROWS,NLAYS )  ! Air density [ kg/m**-3 ]
      REAL QV     ( NCOLS,NROWS,NLAYS )  ! Water vapor mixing ratio [ kg/kg ]

C *** variables computed and output but not carried in CGRID

C *** visibility variables
      INTEGER, PARAMETER :: N_AE_VIS_SPC = 4
      INTEGER, PARAMETER :: IDCVW1 = 1 ! visual range in deciview (Mie)
      INTEGER, PARAMETER :: IBEXT1 = 2 ! extinction [ 1/km ] (Mie)
      INTEGER, PARAMETER :: IDCVW2 = 3 ! visual range in deciview (Reconst)
      INTEGER, PARAMETER :: IBEXT2 = 4 ! extinction [ 1/km ] (Reconst)
      REAL VIS_SPC( NCOLS,NROWS,N_AE_VIS_SPC ) ! Visual range information

C *** aerosol size distribution variables
      REAL DIAM_SPC( NCOLS,NROWS,NLAYS,23 )

C *** atmospheric properties
      REAL XLM             ! atmospheric mean free path [ m ]
      REAL AMU             ! atmospheric dynamic viscosity [ kg/m s ]

C *** N2O5 heterogeneous reaction probability [ ]
      REAL GAMMA_N2O5      ! computed in SUBROUTINE EQL3

C *** mass fraction of each mode less than 2.5um aerodynamic diameter
      REAL PM25( n_mode )  ! fine fractions

C *** visual range information
      REAL BLKDCV1         ! block deciview (Mie)
      REAL BLKEXT1         ! block extinction [ km**-1 ] (Mie)

      REAL BLKDCV2         ! block deciview (Reconstructed)
      REAL BLKEXT2         ! block extinction [ km**-1 ] (Reconstructed)

C *** other internal aerosol variables
      INTEGER IND                         ! index to be used with INDEX1

C *** synchronization time step [ s ]
      REAL DT

C *** variables to set up for "dry transport "
      REAL M3_WET, M3_DRY   ! third moment with and without water
      REAL M2_WET, M2_DRY   ! second moment with and without water
      REAL M3SUBT           ! temp variable for dry 3rd moment calcs

C *** flag to include water in the 3rd moment calculation
      LOGICAL, PARAMETER :: M3_WET_FLAG = .FALSE.

C *** if LIMIT_Sg = T, atkn & accum std. dev. are not changed by GETPAR
      LOGICAL, PARAMETER :: LIMIT_Sg = .FALSE.

C *** variables aerosol diagnostic file flag
      INTEGER      STATUS            ! ENV... status
      CHARACTER( 80 ) :: VARDESC     ! environment variable description

C *** environment variable for AERDIAG file
      CHARACTER( 16 ), SAVE :: CTM_AERDIAG = 'CTM_AERDIAG'

C *** flag for AERDIAG file [F], default
      LOGICAL, SAVE :: AERDIAG

C *** first pass flag 
      LOGICAL, SAVE :: FIRSTIME = .TRUE.

C *** ratio of molecular weights of water vapor to dry air = 0.622015
      REAL, PARAMETER :: EPSWATER = MWWAT / MWAIR

C *** dry moment factor
      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0

C *** Statement Function **************
      REAL ESATL ! arithmetic statement function for vapor pressure [Pa]
      REAL TT
C *** Coefficients for the equation, ESATL defining saturation vapor pressure
      REAL, PARAMETER :: AL = 610.94
      REAL, PARAMETER :: BL = 17.625
      REAL, PARAMETER :: CL = 243.04

C *** values of AL, BL, and CL are from:
C     Alduchov and Eskridge, "Improved Magnus Form Approximations of
C                            Saturation Vapor Pressure,"
C                            Jour. of Applied Meteorology, vol. 35,
C                            pp 601-609, April, 1996.         

      ESATL( TT ) = AL * EXP( BL * ( TT - 273.15 ) / ( TT - 273.15 + CL ) )

C *** End Statement Function  ********

      logical, save :: now = .true.

C ------------------ begin body of AERO_DRIVER -------------------------

      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

C*** Make sure an ae5 version of the mechanism is being used
         IF ( INDEX ( MECHNAME, 'AE5' ) .LE. 0 ) THEN
            XMSG = 'AERO5 requires an AE5 version of chemical mechanism'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT3 )
         END IF

C *** Get aerosol diagnostic file flag.

         AERDIAG = .FALSE.         ! default
         VARDESC = 'Flag for writing the aerosol diagnostic file'
         AERDIAG = ENVYN( CTM_AERDIAG, VARDESC, AERDIAG, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG = 'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

C *** Open the met files.

         IF ( .NOT. OPEN3( MET_CRO_3D, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open  MET_CRO_3D  file '
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         IF ( .NOT. OPEN3( MET_CRO_2D, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open  MET_CRO_2D file '
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

C *** Set up file structure for visibility file. It has two variables,
C     visual range in deciview units (dimensionless) and extinction in
C     units of (1/km) and is for layer 1 only.

         IF ( MYPE .EQ. 0 ) CALL OPVIS ( JDATE, JTIME, TSTEP( 1 ) )

C *** Open the aerosol parameters file (diameters and standard deviations).

         IF ( AERDIAG .AND.
     &        MYPE .EQ. 0 ) CALL OPDIAM ( JDATE, JTIME, TSTEP( 1 ) )

C Get domain decomp info from the MET_CRO_3D file

         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

C Create new adjoint file and open for update
C         CALL CHK_FILE_OPEN ( ADJ_AERO_CHK, JDATE, JTIME, .TRUE., TSTEP )

      END IF    ! FIRSTIME

      MDATE  = JDATE
      MTIME  = JTIME
      MSTEP = TIME2SEC( TSTEP( 2 ) )
      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      IF ( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) WRITETIME = .TRUE.

C *** Set floating point synchronization time step:
      DT = FLOAT( MSTEP ) ! set time step in seconds

C *** Get Meteorological Variables

C *** pressure [Pa]
      VNAME = 'PRES'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, PRES ) ) THEN
         XMSG = 'Could not interpolate '// TRIM( VNAME ) // ' from MET_CRO_3D '
         CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      END IF

C *** temperature [K]
      VNAME = 'TA'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, TA ) ) THEN
         XMSG = 'Could not interpolate '// TRIM( VNAME ) // ' from MET_CRO_3D '
         CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      END IF

C *** specific humidity [g H2O/g air]
      VNAME = 'QV'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, QV ) ) THEN
         XMSG = 'Could not interpolate '// TRIM( VNAME ) // ' from MET_CRO_3D '
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF 

C *** air density [kg/m3]
      VNAME = 'DENS'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MDATE, MTIME, DENS ) ) THEN
         XMSG = 'Could not interpolate '// TRIM( VNAME ) // ' from MET_CRO_3D '
         CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
      END IF

C Checkpoint concentrations
C      CALL CHK_FILE_WRITE ( ADJ_AERO_CHK, JDATE, JTIME, CGRID )

C --------------------- Begin loops over grid cells --------------------------

C *** initialize conc arrays
      !Write(*,'('' before aeroproc '')')
      !Write(*,'(i5,g16.6)') (spc, cgrid(18, 8, 1, spc), spc=1,137)

      DO L = 1, NLAYS
         DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS

C *** Grid cell meteorological data.
            AIRTEMP  = TA   ( C,R,L )
            AIRPRS   = PRES ( C,R,L )   ! Note pascals
            AIRQV    = QV   ( C,R,L )
            AIRDENS  = DENS ( C,R,L )
            H2OSATVP = ESATL( AIRTEMP ) 
            H2OVP    = AIRPRS * AIRQV / ( EPSWATER  + AIRQV )
            AIRRH    = MAX( 0.005, MIN( 0.99, H2OVP / H2OSATVP ) )

C *** extract grid cell concentrations of aero species from CGRID
C     into aerospc_conc in aero_data module
            CALL EXTRACT_AERO( CGRID( C,R,L,: ) )

C *** extract grid cell concentrations of gas precursors from CGRID (ppm)
C     into precursr_conc in precursor_data
            CALL EXTRACT_PRECURSOR( CGRID( C,R,L,: ) )

C**   Calculate SO4RATE stored in module
            SO4RATE = PRECURSOR_CONC( SULPRD_IDX ) / DT

C *** extract soa concentrations from CGRID
            CALL EXTRACT_SOA( CGRID( C,R,L,: ) )

C *** aerosol process routines
            CALL AEROPROC( DT, C, R, L, GAMMA_N2O5 )

C *** update aerosol variables conc back into CGRID
            CALL UPDATE_AERO( CGRID( C,R,L,: ) )

C *** update precursor variables conc back into CGRID
            CALL UPDATE_PRECURSOR( CGRID( C,R,L,: ) )

C *** update gas soa concentrations back to CGRID
            CALL UPDATE_SOA( CGRID( C,R,L,: ) )

C *** Calculate volume fraction of each mode < 2.5um aerodynamic diameter
            DO N = 1, N_MODE
               CALL INLET25 ( AEROMODE_DIAM( N ), AEROMODE_SDEV( N ),
     &                        AEROMODE_DENS( N ), PM25( N ) )
            END DO

C *** Write aerosol extinction coefficients and deciviews to visibility
C     diagnostic array (lowest vertical layer only)

            IF ( WRITETIME .AND. L .EQ. 1 ) THEN

               CALL GETVISBY( BLKDCV1, BLKEXT1, BLKDCV2, BLKEXT2 )

               VIS_SPC( C,R,IDCVW1 ) = BLKDCV1  ! visual range (deciview)(Mie)

               VIS_SPC( C,R,IBEXT1 ) = BLKEXT1  ! aero extinction (1/km)(Mie)

               VIS_SPC( C,R,IDCVW2 ) = BLKDCV2  ! visual range (deciview)
                                                ! (Reconstructed)

               VIS_SPC( C,R,IBEXT2 ) = BLKEXT2  ! aero extinction (1/k)
                                                ! (Reconstructed)
            END IF

C *** Write wet diameters, 2nd, and 3rd moments to aerosol diagnostic array
C     This assumes that GETPAR was last called with M3_WET_FLAG = .TRUE.

            IF ( WRITETIME .AND. AERDIAG ) THEN

               DIAM_SPC( C,R,L, 7 ) = AEROMODE_DIAM( 1 ) * 1.0E6  ! wet i-mode diameter
               DIAM_SPC( C,R,L, 8 ) = AEROMODE_DIAM( 2 ) * 1.0E6  ! wet j-mode diameter
               DIAM_SPC( C,R,L, 9 ) = AEROMODE_DIAM( 3 ) * 1.0E6  ! wet k-mode diameter
               DIAM_SPC( C,R,L,10 ) = MOMENT2_CONC( 1 )           ! wet i-mode 2nd moment
               DIAM_SPC( C,R,L,11 ) = MOMENT2_CONC( 2 )           ! wet j-mode 2nd moment
               DIAM_SPC( C,R,L,12 ) = MOMENT2_CONC( 3 )           ! wet k-mode 2nd moment
               DIAM_SPC( C,R,L,16 ) = MOMENT3_CONC( 1 )           ! wet i-mode 3rd moment
               DIAM_SPC( C,R,L,17 ) = MOMENT3_CONC( 2 )           ! wet j-mode 3rd moment
               DIAM_SPC( C,R,L,18 ) = MOMENT3_CONC( 3 )           ! wet k-mode 3rd moment
               DIAM_SPC( C,R,L,20 ) = PM25( 1 )                   ! i-mode fine fraction
               DIAM_SPC( C,R,L,21 ) = PM25( 2 )                   ! j-mode fine fraction
               DIAM_SPC( C,R,L,22 ) = PM25( 3 )                   ! coarse-mode fine fraction
               DIAM_SPC( C,R,L,23 ) = GAMMA_N2O5                  ! N2O5 heterorxn probability
            END IF   ! WRITETIME .AND. AERDIAG

C *** Calculate 2nd and 3rd moments of the "dry" aerosol distribution
C     NOTE! "dry" aerosol excludes both H2O and SOA  (Jan 2004 --SJR)
C     EXCEPT!  nonvolatile SOA is part of dry aerosol (Oct 2007 --PVB)

C     Aitken mode.
            M3_WET = MOMENT3_CONC( 1 )
            M3SUBT = ( 1.0E-9 * F6DPI / AEROSPC( AH2O_IDX )%DENSITY )
     &             * AEROSPC_CONC( AH2O_IDX,1 )
            M3_DRY = M3_WET - M3SUBT
            M2_WET = MOMENT2_CONC( 1 )
            M2_DRY = M2_WET * ( M3_DRY / M3_WET ) ** TWOTHIRDS

            MOMENT3_CONC( 1 ) = MAX( CONMIN, M3_DRY )
            MOMENT2_CONC( 1 ) = MAX( CONMIN, M2_DRY )

C     accumulation mode.
            M3_WET = MOMENT3_CONC( 2 )

            M3SUBT = ( 1.0E-9 * F6DPI / AEROSPC(AH2O_IDX)%DENSITY )
     &             * AEROSPC_CONC( AH2O_IDX,2 )
            DO SPC = 1, N_VAPOR
               M3SUBT = M3SUBT
     &                + ( 1.0E-9 * F6DPI / AEROSPC(AORGPA_IDX)%DENSITY )
     &                * AEROSPC_CONC( SOA_AEROMAP( SPC ),2 )
            END DO

            M3_DRY = M3_WET - M3SUBT
            M2_WET = moment2_conc( 2 )
            M2_DRY = M2_WET * ( M3_DRY / M3_WET ) ** TWOTHIRDS

            MOMENT3_CONC( 2 ) = MAX( CONMIN, M3_DRY )
            MOMENT2_CONC( 2 ) = MAX( CONMIN, M2_DRY )

C     coarse mode
            M3_WET = MOMENT3_CONC(3)
            M3SUBT = ( 1.0E-9 * F6DPI / AEROSPC( AH2O_IDX )%DENSITY )
     &             * AEROSPC_CONC( AH2O_IDX,3 )
            M3_DRY = M3_WET - M3SUBT
            M2_WET = MOMENT2_CONC( 3 )
            M2_DRY = M2_WET * ( M3_DRY / M3_WET ) ** TWOTHIRDS

            MOMENT3_CONC( 3 ) = MAX( CONMIN, M3_DRY )
            MOMENT2_CONC( 3 ) = MAX( CONMIN, M2_DRY )

C *** Calculate geometric mean diameters and standard deviations of the
C     "dry" size distribution

            CALL GETPAR( M3_WET_FLAG, LIMIT_Sg )

C *** Write dry aerosol distribution parameters to aerosol diagnostic array

            IF ( WRITETIME .AND. AERDIAG ) THEN

               DIAM_SPC( C,R,L, 1 ) = EXP( AEROMODE_SDEV( 1 ) )
               DIAM_SPC( C,R,L, 2 ) = EXP( AEROMODE_SDEV( 2 ) )
               DIAM_SPC( C,R,L, 3 ) = EXP( AEROMODE_SDEV( 3 ) ) 
               DIAM_SPC( C,R,L, 4 ) = AEROMODE_DIAM( 1 ) * 1.0E6  ! dry i-mode diameter
               DIAM_SPC( C,R,L, 5 ) = AEROMODE_DIAM( 2 ) * 1.0E6  ! dry j-mode diameter
               DIAM_SPC( C,R,L, 6 ) = AEROMODE_DIAM( 3 ) * 1.0e6  ! dry k-mode diameter
               diam_spc( c,r,l,13 ) = MOMENT3_CONC( 1 )           ! dry i-mode 3rd moment
               DIAM_SPC( C,R,L,14 ) = MOMENT3_CONC( 2 )           ! dry j-mode 3rd moment
               DIAM_SPC( C,R,L,15 ) = MOMENT3_CONC( 3 )           ! dry k-mode 3rd moment
               DIAM_SPC( C,R,L,19 ) = AIRRH         ! relative humidity
            END IF   ! WRITETIME .AND. AERDIAG

C *** Calculate aerosol surface area from the dry 2nd moment.  Dry value is
C     used in transport routines.  Put aeroSrf values back to CGRID

            do SIZE = 1, N_MODE
               N = AEROSRF_MAP( SIZE )
               CGRID( C,R,L,N ) =  PI * MOMENT2_CONC( SIZE )
            END DO

            END DO ! loop on MY_COLS
         END DO ! loop on MY_ROWS
      END DO ! loop on NLAYS

      !Write(*,'('' after aeroproc '')')
      !Write(*,'(i5,g16.6)') (spc, cgrid(18, 8, 1, spc), spc=1,137)

C *** If last call this hour, write visibility information.

      IF ( WRITETIME ) THEN
         MDATE = JDATE
         MTIME = JTIME
         CALL NEXTIME ( MDATE, MTIME, TSTEP( 2 ) )
         WSTEP = 0
         WRITETIME = .FALSE.

         IF ( .NOT. WRITE3( CTM_VIS_1, ALLVAR3,
     &                      MDATE, MTIME, VIS_SPC ) ) THEN
            XMSG = 'Could not write ' // CTM_VIS_1 // ' file'
            CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
         END IF

         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &                  'Timestep written to', CTM_VIS_1,
     &                  'for date and time', MDATE, MTIME
  
C *** Write data to the aerosol diagnostic file.

         IF ( AERDIAG ) THEN

            IF ( .NOT. WRITE3( CTM_DIAM_1, ALLVAR3,
     &                         MDATE, MTIME, DIAM_SPC ) ) THEN
               XMSG = 'Could not write ' // CTM_DIAM_1 // ' file'
               CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

            WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &                     'Timestep written to', CTM_DIAM_1,
     &                     'for date and time', MDATE, MTIME

         END IF

      END IF

      RETURN
      END
