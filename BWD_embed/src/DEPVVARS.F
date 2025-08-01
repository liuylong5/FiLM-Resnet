
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_bwd/DEPVVARS.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE DEPVVARS

C-----------------------------------------------------------------------
C Name:     Dry Deposition Variables
C Purpose:  Contains arrays specific to dry deposition species.
C           Initializes dry deposition arrays.
C Revised:  19 Aug 2005  Original version.  (T. Otte and W. Hutzell)
C           25 Nov 2006 J.Young: combine
C           30 Apr 2008  Removed references to RADMDRY.  Added five air
C                        toxic species to LTOTG.  (T. Otte and W. Hutzell)
C-----------------------------------------------------------------------

      INTEGER, PARAMETER :: LDDEP = 4
      INTEGER, PARAMETER :: LTOTG = 54
      CHARACTER( 16 )    :: DEPV_METHOD

C-----------------------------------------------------------------------
C General dry deposition arrays.
C-----------------------------------------------------------------------

      INTEGER, SAVE         :: N_USE_DEPSPC        
      LOGICAL, SAVE         :: USE_DEPSPC( LTOTG )
      CHARACTER( 16 ), SAVE :: DEPSPC( LTOTG )

C----------------------------------------------------------------------
C local deposition velocites [ m/s ]
C----------------------------------------------------------------------

      CONTAINS

         FUNCTION DEPVVARS_INIT( ) RESULT ( SUCCESS )

         USE HGRD_DEFN

         IMPLICIT NONE

         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    ! I/O definitions and declarations

         LOGICAL SUCCESS

         INTEGER ALLOCSTAT
         CHARACTER( 96 ) :: XMSG

         SUCCESS = .TRUE.

C-------------------------------------------------------------------------------
C For M3DRY, set up core species, and include toxic and chlorine compounds.
C-------------------------------------------------------------------------------

         DEPV_METHOD = 'M3DRY'

         DEPSPC(  1 ) = 'SO2'
         DEPSPC(  2 ) = 'SULF'
         DEPSPC(  3 ) = 'NO2'
         DEPSPC(  4 ) = 'NO'
         DEPSPC(  5 ) = 'O3'
         DEPSPC(  6 ) = 'HNO3'
         DEPSPC(  7 ) = 'H2O2'
         DEPSPC(  8 ) = 'ALD'
         DEPSPC(  9 ) = 'HCHO'
         DEPSPC( 10 ) = 'OP'
         DEPSPC( 11 ) = 'PAA'
         DEPSPC( 12 ) = 'ORA'
         DEPSPC( 13 ) = 'NH3'
         DEPSPC( 14 ) = 'PAN'
         DEPSPC( 15 ) = 'HONO'
         DEPSPC( 16 ) = 'CO'
         DEPSPC( 17 ) = 'METHANOL'
         DEPSPC( 18 ) = 'N2O5'
         DEPSPC( 19 ) = 'NO3'
         DEPSPC( 20 ) = 'GEN_ALD'
         DEPSPC( 21 ) = 'CL2'
         DEPSPC( 22 ) = 'HOCL'
         DEPSPC( 23 ) = 'HCL'
         DEPSPC( 24 ) = 'FMCL'
         DEPSPC( 25 ) = 'ICL1'
         DEPSPC( 26 ) = 'ICL2'
         DEPSPC( 27 ) = 'HG'
         DEPSPC( 28 ) = 'HGIIGAS'
         DEPSPC( 29 ) = 'TECDD_2378' 
         DEPSPC( 30 ) = 'PECDD_12378'
         DEPSPC( 31 ) = 'HXCDD_123478' 
         DEPSPC( 32 ) = 'HXCDD_123678'  
         DEPSPC( 33 ) = 'HXCDD_123789'  
         DEPSPC( 34 ) = 'HPCDD_1234678'  
         DEPSPC( 35 ) = 'OTCDD'  
         DEPSPC( 36 ) = 'TECDF_2378'  
         DEPSPC( 37 ) = 'PECDF_12378'  
         DEPSPC( 38 ) = 'PECDF_23478'  
         DEPSPC( 39 ) = 'HXCDF_123478'  
         DEPSPC( 40 ) = 'HXCDF_123678'  
         DEPSPC( 41 ) = 'HXCDF_234678'  
         DEPSPC( 42 ) = 'HXCDF_123789'  
         DEPSPC( 43 ) = 'HPCDF_1234678'  
         DEPSPC( 44 ) = 'HPCDF_1234789'  
         DEPSPC( 45 ) = 'OTCDF'  
         DEPSPC( 46 ) = 'NAPHTHALENE'  
         DEPSPC( 47 ) = '1NITRONAPHTHA'  
         DEPSPC( 48 ) = '2NITRONAPHTHA'  
         DEPSPC( 49 ) = '14NAPHTHOQUIN'  
         DEPSPC( 50 ) = 'HEXMETH_DIIS'
         DEPSPC( 51 ) = 'HYDRAZINE'
         DEPSPC( 52 ) = 'MAL_ANHYDRIDE'
         DEPSPC( 53 ) = 'TOLUENE_DIIS'
         DEPSPC( 54 ) = 'TRIETHYLAMINE'

         N_USE_DEPSPC = 0
         USE_DEPSPC = .FALSE.   ! array

         RETURN

         END FUNCTION DEPVVARS_INIT

      END MODULE DEPVVARS
