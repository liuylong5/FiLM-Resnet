
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/par/par/par_term.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE PAR_TERM ( CALLER, JDATE, JTIME, BEGTIME )
C ....................................................................
 
C  PURPOSE:    Performs parallel related shutdown and cleanup
C              necessary for the MPI, I/O, and stencil exchange
C              libraries. Writes out elapsed time.
 
C  REVISION HISTORY: 
C       Original version  2/96 by Al Bourgeois for parallel implementaion.
C       Modified 07/29/1998 by AJB to shut down IOAPI. 
C       Modified 08/24/1998 by AJB to remove from PAR_IO and put in CCTM. 
C       Modified 03/09/2001 by David Wong
C         -- re-arrange the call sequence of MPI_FINALIZE so SHUT3 can be
C            called correctly
 
C  ARGUMENT LIST DESCRIPTION: 
C  IN:
C      CHARACTER*16    CALLER    ! Program name.
C      INTEGER         JDATE     ! Current model date, coded YYYYDDD.
C      INTEGER         JTIME     ! Current model time, coded HHMMSS.
C      REAL            BEGTIME   ! Wall-clock time (sec) at initialization.
 
C  OUT: none
 
C  LOCAL VARIABLE DESCRIPTION:  see below
 
C  CALLS:  PIO_TERM, ELAPTIM, M3WARN, MPI_FINALIZE 
 
C .......................................................................

      IMPLICIT  NONE

C INCLUDE FILES

#include      SUBST_IODECL       ! I/O definitions and declarations
      INCLUDE SUBST_MPICH        ! MPI definitions and parameters.

C ARGUMENTS:

      CHARACTER( 16 ) :: CALLER  ! Program name.
      INTEGER         JDATE      ! Current model date, coded YYYYDDD.
      INTEGER         JTIME      ! Current model time, coded HHMMSS.
      REAL            BEGTIME    ! Wall-clock time (sec) at initialization.

C LOCAL VARIABLES: 

      INTEGER         MYPE       ! Local processor id (range is 0 to NPROCS-1)
      INTEGER         IERR       ! Return code from MPI_ calls
      REAL            ENDTIME    ! Wall-clock time (sec) at termination
      INTEGER         STRLEN     ! string length
      CHARACTER( 80 ) :: MSG = ' '   ! Message issued from M3WARN routine
      CHARACTER(  3 ) :: CMYPE   ! Processor ID string
      CHARACTER(  7 ) :: PE_STR  ! String suffix to go with processor ID
      CHARACTER( 16)  :: CALL16  ! First 16 characters of CALLER
      CHARACTER( 26 ) :: PCALLER ! Program/processor name

      INTEGER         ELAPTIME
      INTEGER         LOGDEV
      INTEGER, EXTERNAL :: TRIMLEN ! M3IO lib
      CHARACTER( 24 ), EXTERNAL :: DT2STR ! Construct HH:MM:SS string - M3IO lib

C .......................................................................
C     begin function PAR_TERM( CALLER, JDATE, JTIME, BEGTIME )

C Get final wall-clock time

      ENDTIME = MPI_WTIME()

C Write elapsed time to log 

!     CALL ELAPTIM( BEGTIME, ENDTIME )
      ELAPTIME = NINT ( ENDTIME - BEGTIME )
      LOGDEV = INIT3()
      WRITE( LOGDEV,92001 ) ELAPTIME

C Get my processor rank (i.e., get ID of this processor)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, MYPE, IERR )
      IF ( IERR .NE. 0 ) THEN
         MSG = 'Error in MPI_COMM_WORLD.'
         CALL M3EXIT ( 'PAR_TERM', 0, 0, MSG, 1 )
         END IF

C Shut down IOAPI

!     CALL PIO_TERM( CALLER, JDATE, JTIME )

C Shut down program

!     CALL PM3EXIT( CALLER, JDATE, JTIME, MSG, 0 )

C Construct PCALLER string

      WRITE ( PE_STR,'( A7 )' ) ' on PE '
      WRITE( CMYPE,'( I3.3 )' ) MYPE
      STRLEN = MIN ( 16, TRIMLEN( CALLER ) )
      CALL16 = CALLER( 1:STRLEN )
      PCALLER = CALL16( 1:STRLEN ) // PE_STR // CMYPE
      WRITE( LOGDEV,92010 ) PCALLER

      IF ( JDATE .GT. 0 .OR. JTIME .GT. 0 )
     &   WRITE( LOGDEV,92020 ) DT2STR( JDATE, JTIME ), JDATE, JTIME

C Shut down IOAPI

      IF ( .NOT. SHUT3() ) THEN
         WRITE( LOGDEV,91000 )
     &         'ERROR ABORT in program *** ', PCALLER, ' ***',
     &         'Could not shut down M3IO correctly'
         CALL M3EXIT ( 'PAR_TERM', 0, 0, MSG, 1 )
      END IF

C Clean up communications

      CALL MPI_FINALIZE ( IERR )
      IF ( IERR .NE. 0 ) THEN
         MSG = 'Error terminating MPI'
         CALL M3EXIT ( 'PAR_TERM', 0, 0, MSG, 1 )
         END IF

      RETURN

91000 FORMAT( // 5X , A, A, A, / 5X , A // )
92001 FORMAT( / 'The elapsed time for this job was', I7.6, ' seconds.' )
92010 FORMAT( // 5X , '--->> Normal Completion of program ', A, ' <<---' )
92020 FORMAT( 5X , 'Date and time ', A, ' (', I7, ':', I6.6, ')' )

      END
