
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
C $Header: /opt/home/shunliu/src/cmaq/4.7.1/models/CCTM/src/par/par/par_init.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PAR_INIT( COLROW, NSPCS, CLOCK, IERROR )
C.......................................................................
C
C Purpose:   Performs initialization necessary for parallel execution
C            of the CCTM program. This consists of initialization of
C            a Message Passing Interface (MPI) library, a parallel
C            Models-3 I/O API library, and a stencil exchange library.
 
C Revision History: 
C      Original version  2/96 by Al Bourgeois for parallel implementation.
C      Modified 6/98 by AJB for PAR_IO library, added error code.
C      Modified 07/08/1998 by AJB to set IO_GRP and MY_IO_PE.
C      Modified 07/29/1998 by AJB to allow setting of NPROW, NPCOL outside.
C               Also added synchronization of return code across processors.
C      Modified 08/24/1998 by AJB to move this from the PAR_IO library to
C               the CCTM.
C      Modified 09/03/1998 by AJB to call SUBDMAP, and do restructuring.
C      Modified 06/14/1999 by AJB to set env variable ALL_PE_SYNC_IO_MODE
C               for pinterp3.
C      Modified 06/16/1999 by Al Bourgeois to remove interprocessor
C         synchronization. This removes the guarantee that all processors
C         return the same error code, and a "hang" state can occur if
C         PM3EXIT is not called on the condition that this function fails.
C      Modified 08/05/1999 by Al Bourgeois to call ALL_PE_IO instead of
C         PUTENV. Also, remove IERROR argument to SUBDMAP.
C      Modified 22 Nov 00 J.Young: To use Dave Wong's f90 stenex
C      Modified 12/27/2000 by Jeff Young to reorder col/row for call to SE_INIT
C          to correspond to DATA_ORI. (PAR_INIT's arguments were also re-
C          ordered, which wasn't really necessary).
C      Modified 07 Feb 01 J.Young: enable dfio and remove pario libraries
C      Modified 14 Feb 01 David Wong: for hgrd dynamic allocation
C      Modified 20 Mar 01 J.Young: move processor-to-subdomain map to HGRD_DEFN
C               16 Apr 03 J.Young: #ifdef for distr run time env
C               31 Jan 05 J.Young: call GRID_INIT for both horz and vert setup
C               16 Jan 06 J.Young: possible PA grid init
C               15 Feb 06 J.Young: rename PIO_INIT to PIO_RE_INIT to preserve
C                                  backward compatiblity
C               25 May 06 J.Young: ck cols/rows arg list order in pio_init
 
C Argument List Description:
C  In:
C     CHARACTER( 2 ) :: COLROW  ! columns/rows arg list check
C     INTEGER NLAYS             ! Number of layers in entire grid
C     INTEGER NSPCS             ! Number of species in CGRID
 
C  Out:
C     REAL    CLOCK             ! Wall-clock time (sec) at MPI initialization
C     INTEGER IERROR            ! Error code: 0 for ok, 1 for fail
 
C Local Variable Description:
 
C  Other Local Variable Descriptions:  see below

C Calls: M3WARN, SE_INIT, MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_WTIME
 
C Notes: Error flag IERROR is set to 1 if an error occurs in any of
C        the initialization routines called fails.
C........................................................................

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE PAGRD_DEFN            ! PA horiz domain specs
      USE SUBST_MODULES         ! stenex
!     USE SUBST_INIT_MODULE     ! stenex

      USE SWAP_SANDIA

      IMPLICIT  NONE

C Include Files

      INCLUDE SUBST_MPICH       ! MPI definitions and parameters
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_PACTL_ID    ! PA control parameters

C Arguments

      CHARACTER( 2 ), INTENT( IN ) :: COLROW ! columns/rows arg list check
      INTEGER, INTENT( IN )        :: NSPCS  ! Number of species in CGRID
      REAL, INTENT( OUT ) :: CLOCK ! Wall-clock time (sec) at MPI initialization
      INTEGER, INTENT( OUT )       :: IERROR ! Error code

C Local Variables 

      INTEGER        NPROCS      ! Total number of processors
      CHARACTER( 80 ) :: MSG     ! Message issued from M3WARN routine
      CHARACTER( 16 ) :: PNAME = 'PAR_INIT'

      CHARACTER( 2 ) :: DATA_ORI
      INTEGER     GEO_ORI        ! 0 (Cartesian), 1 (Matrix)

      LOGICAL, EXTERNAL :: PIO_RE_INIT
C........................................................................

C Initialize error code
      IERROR = 0

C Set orientation
      CALL UPCASE ( COLROW )
      IF ( COLROW .NE. 'CR' ) THEN
         MSG = 'Argument list must have cols/rows orientation.'
         CALL M3WARN ( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF
      DATA_ORI = COLROW
      GEO_ORI  = 0

C Start up MPI
      CALL MPI_INIT( IERROR )
      IF ( IERROR .NE. 0 ) THEN
         MSG = 'Error in MPI_INIT.'
         CALL M3WARN ( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF

C Get wall-clock time
      CLOCK = REAL( MPI_WTIME() )

C Get number of processors
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NPROCS, IERROR )
      IF ( IERROR .NE. 0 ) THEN
         MSG = 'Error in MPI_COMM_SIZE.'
         CALL M3WARN ( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF

C Get my processor rank (i.e., get ID of this processor)
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYPE, IERROR )
      IF ( IERROR .NE. 0 ) THEN
         MSG = 'Error in MPI_COMM_WORLD.'
         CALL M3WARN ( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF

#ifdef cluster
C Distribute run time environment variables to other machines - by bo wang
      write( *,* ) ' MYPE -> NPROCS: ', MYPE, NPROCS
      CALL distr_env ( MYPE, NPROCS )
#endif

C Set up horizontal domain and calculate processor-to-subdomain maps
C and define vertical layer structure (in module GRID_CONF)
      IF ( .NOT. GRID_INIT ( NPROCS, MYPE ) ) THEN
         MSG = '*** Failure defining domain configuration'
         CALL M3WARN ( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF

C Initialize parallel I/O
      IF ( .NOT. PIO_RE_INIT( DATA_ORI, GL_NCOLS, GL_NROWS, NLAYS, NTHIK,
     &                        NCOLS, NROWS, NPCOL, NPROW, NPROCS,
     &                        MYPE, .TRUE. ) ) THEN
         MSG = 'Failed to initialize parallel I/O library.'
         CALL M3WARN( PNAME, 0, 0, MSG )
         IERROR = 1; RETURN
         END IF

      CALL SWAP_INIT_SNL( NPCOL, NPROW )

C Initialize stencil exchange
      CALL SE_INIT( NPROCS, NPCOL, NPROW, GL_NCOLS, GL_NROWS, NLAYS,
     &              NSPCS, MYPE, MNDIS, MEDIS, MSDIS, MWDIS,
     &              DATA_ORI, GEO_ORI )

C Set up horizontal domain and calculate processor-to-subdomain maps for
C process analysis, if required

      IF ( LIPR .OR. LIRR ) THEN
         IF ( .NOT. PAGRD_INIT( NPROCS, MYPE ) ) THEN
            MSG = '*** Failure defining PA domain configuration'
            CALL M3WARN ( PNAME, 0, 0, MSG )
            IERROR = 1; RETURN
            END IF
         END IF

D     write( *,* ) ' MYPE -> NPROCS: ',   MYPE, NPROCS
D     write( *,* ) ' MYPE -> NPCOL: ',    MYPE, NPCOL
D     write( *,* ) ' MYPE -> NPROW: ',    MYPE, NPROW
D     write( *,* ) ' MYPE -> GL_NCOLS: ', MYPE, GL_NCOLS
D     write( *,* ) ' MYPE -> GL_NROWS: ', MYPE, GL_NROWS
D     write( *,* ) ' MYPE -> NLAYS: ',    MYPE, NLAYS
D     write( *,* ) ' MYPE -> NSPCS: ',    MYPE, NSPCS
D     write( *,* ) ' MYPE -> MNDIS: ',    MYPE, MNDIS
D     write( *,* ) ' MYPE -> MEDIS: ',    MYPE, MEDIS
D     write( *,* ) ' MYPE -> MSDIS: ',    MYPE, MSDIS
D     write( *,* ) ' MYPE -> MWDIS: ',    MYPE, MWDIS
D     write( *,* ) ' MYPE -> DATA_ORI: ', MYPE, DATA_ORI
D     write( *,* ) ' MYPE -> GEO_ORI: ',  MYPE, GEO_ORI

      RETURN
      END
