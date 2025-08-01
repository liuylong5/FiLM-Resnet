
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

! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/CGRID_SPCS.F,v 1.1 2010/11/17 19:22:51 sergey Exp $ 

! what(1) key, module and SID; SCCS file; date and time of last delta:
! %W% %P% %G% %U%

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE CGRID_SPCS

!-----------------------------------------------------------------------
! Function:
!   Define extent and starting positions of the species variables in the
!   CGRID array

! Preconditions:

! Subroutines and functions called:

! Revision history:
!   Dec 00 - Jeff - Initial implementation
!-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, SAVE :: NSPCSD   ! Number of species in CGRID
      INTEGER, SAVE :: GC_STRT  ! Starting index of gas chemistry species
      INTEGER, SAVE :: AE_STRT  ! Starting index of aerosol species
      INTEGER, SAVE :: NR_STRT  ! Starting index of non-reactive species
      INTEGER, SAVE :: TR_STRT  ! Starting index of tracer species

      CONTAINS

         SUBROUTINE CGRID_MAP ( NSPCSD, GC_STRT, AE_STRT, NR_STRT, TR_STRT )

         INCLUDE SUBST_GC_SPC      ! gas chemistry species table
         INCLUDE SUBST_AE_SPC      ! aerosol species table
         INCLUDE SUBST_NR_SPC      ! non-reactive species table
         INCLUDE SUBST_TR_SPC      ! tracer species table

         INTEGER, INTENT( OUT ) :: NSPCSD
         INTEGER, INTENT( OUT ) :: GC_STRT
         INTEGER, INTENT( OUT ) :: AE_STRT
         INTEGER, INTENT( OUT ) :: NR_STRT
         INTEGER, INTENT( OUT ) :: TR_STRT

         LOGICAL, SAVE :: FIRSTIME = .TRUE.

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.

            NSPCSD = N_GC_SPCD + N_AE_SPC + N_NR_SPC + N_TR_SPC
 
            GC_STRT = 1 ! Always, even if N_GC_SPCS = 0
 
            AE_STRT = N_GC_SPCD + 1
 
            NR_STRT = N_GC_SPCD + N_AE_SPC + 1
 
            TR_STRT = N_GC_SPCD + N_AE_SPC + N_NR_SPC + 1

            END IF

         RETURN
         END SUBROUTINE CGRID_MAP

      END MODULE CGRID_SPCS
