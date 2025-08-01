
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/procan/pa/pa_output.F,v 1.1.1.1 2010/06/14 16:03:06 sjr Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)pa_output.F	1.1 /project/mod3/CMAQ/src/procan/pa/SCCS/s.pa_output.F 07 Jul 1997 10:20:55

       SUBROUTINE PA_OUTPUT( CGRID, JDATE, JTIME )

C**********************************************************************
C
C  FUNCTION: To output the Process Analysis and Integrated Reaction
C            Rate data           
C
C  PRECONDITIONS: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: INIT3, ENVSTR, M3EXIT
C
C  REVISION HISTORY: Prototype created by Jerry Gipson, August, 1996
C                    allow env var for file names Jeff, Dec, 1996
C                    Modified May, 1997 by Jerry Gipson to be consistent
C                      with beta CTM
C                    Modified Sept, 1997 by Jerry Gipson to be consistent
C                      with targeted CTM
C                    Modified 1/19/99 by David Wong at LM:
C                      -- add four include files because of new PA_CMN.EXT
C                      -- add DATA_COPY function call to redistribute PA grid
C                    Modified 2/26/99 by David Wong at LM:
C                      -- remove SUBST_AE_SPC, SUBST_NR_SPC, SUBST_TR_SPC,
C                         three .EXT files
C                      -- replaced DATA_COPY function with dimension specific
C                         DATA_COPY function and modified its argument list
C                      -- used ifdef statement to distinguish parallel
C                         implementation of IRR calculation which does not
C                         start at the origin
C                    Modified 22 Nov 00 by J.Young: Dave Wong's f90 stenex
C                         DATA_COPY
C                    30 Mar 01 J.Young: dyn alloc - Use PAGRD_DEFN, which Uses
C                         HGRD_DEFN; requires DBUFF for WRITE3
C                    10 Oct 01 David Wong
C                      -- use DBUFF to hold the output data for IRR data since
C                         parallel WRITE3 can't output subsection of data
C                      -- removed the usage of DBUFF for PA output
C                      -- used a new formula to compute number of indexes in
C                         each IPR and IRR output file rather than file header
C                         to remove the call of DESC3 and use IPRFNAM and 
C                         IRRNAME to determine variable name
C                    10 Sep 02 J.Young: fix bug in IPRNAME, IRRNAME index
C                    31 Jan 05 J.Young: dyn alloc - establish both horizontal &
C                       vertical domain specifications in one module (GRID_CONF)
C                    23 Aug 05 David Wong
C                      -- inserted PIO_INIT subroutine calls before and after
C                         outputing PA data which is smaller than the model
C                         domain size, to compute data indices correctly in each
C                         processor for the entire PARIO library functions
C                    15 Feb 06 J.Young: new pario PIO_RE_INIT from PIO_INIT to
C                         add flag to to control printing decomposition map

C**********************************************************************

      USE GRID_CONF             ! horizontal & vertical domain specifications
      USE PAGRD_DEFN            ! PA horiz domain specs
      USE SUBST_MODULES              ! stenex
!     USE SUBST_DATA_COPY_MODULE     ! stenex

      IMPLICIT NONE 

C..INCLUDES:
!     INCLUDE SUBST_HGRD_ID     ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID     ! vertical dimensioning parameters
      INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC     ! file header data structurer
#include      SUBST_IODECL      ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID    ! file name parameters
!     INCLUDE SUBST_BLKPRM      ! cell blocking parameters

      INCLUDE SUBST_PACTL_ID    ! PA control parameters
      INCLUDE SUBST_PACMN_ID    ! PA common block definitions

C..ARGUMENTS:
!     REAL :: CGRID( :,:,:,: )
      REAL, POINTER :: CGRID( :,:,:,: )
      INTEGER JDATE    ! current model date, format YYYYDDD
      INTEGER JTIME    ! current model time, format HHMMSS

C..PARAMETERS: None 

C..EXTERNAL FUNCTIONS: None 

C..SAVED LOCAL VARIABLES: 

      CHARACTER( 16 ), SAVE :: PNAME = 'PA_OUTPUT'
      CHARACTER( 16 ), SAVE :: IPRFNAM( MXFILE3 )  ! IPR output file names
      CHARACTER( 16 ), SAVE :: IRRFNAM( MXFILE3 )  ! IRR output file names
      INTEGER, SAVE :: NIPRFLS          ! Number of IPR output files
      INTEGER, SAVE :: NIRRFLS          ! Number of IRR output files
      INTEGER, SAVE :: LOGDEV           ! Unit number of output log
      LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call

C..SCRATCH LOCAL VARIABLES:
      CHARACTER( 80 ) ::  XMSG = ' '    ! Error message buffer

      INTEGER C         ! Loop index for columns
      INTEGER ICG       ! Index for species in cgrid array
      INTEGER L         ! Loop index for layers
      INTEGER NIRR      ! Loop index for IRR output variables
      INTEGER NFL       ! Loop index for no. of files
      INTEGER NOUT      ! Loop index for PA output variables
      INTEGER NS        ! Loop index for number cgrid concs saved
      INTEGER PC        ! Index for PA output column
      INTEGER PL        ! Index for PA output level
      INTEGER PR        ! Index for PA output row
      INTEGER R         ! Loop index for rows
      INTEGER VARINDX   ! Index for output variable
      INTEGER NVARINDX  ! number of output variable
     
#ifdef parallel
      REAL TIRROUT( NCOLS,NROWS,NLAYS )
      CHARACTER( 2 ) :: COLROW = 'CR'  ! col/row arg list order for pio_re_init
      LOGICAL, EXTERNAL :: PIO_RE_INIT
      LOGICAL, SAVE :: WFLG = .TRUE.   ! pio_init - print decomposition
#endif
      REAL, ALLOCATABLE :: DBUFF( :,:,: )        ! write buffer
      INTEGER ALLOCSTAT

C**********************************************************************

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  On first call, get and save  the output file names 
      IF( LFIRST ) THEN

         LFIRST = .FALSE.
         LOGDEV = INIT3( )

         IF( NIPRVAR. GT. 0 ) THEN

            IF( MOD( NIPRVAR, MXVARS3 ) .EQ. 0 ) THEN
               NIPRFLS = NIPRVAR / MXVARS3
            ELSE
               NIPRFLS = ( NIPRVAR / MXVARS3 ) + 1
            ENDIF

            DO NFL = 1, NIPRFLS

               IF( NFL .EQ. 1 ) THEN
                  IPRFNAM( NFL ) = CTM_IPR_1
               ELSEIF( NFL .EQ. 2 ) THEN
                  IPRFNAM( NFL ) = CTM_IPR_2
               ELSEIF( NFL .EQ. 3 ) THEN
                  IPRFNAM( NFL ) = CTM_IPR_3
               ELSEIF( NFL .GT. 3 ) THEN
                  XMSG = 'Maximum number of IPR output files exceeded'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               ENDIF

            ENDDO

         ENDIF

         IF( NIRRVAR. GT. 0 ) THEN

            IF( MOD( NIRRVAR, MXVARS3 ) .EQ. 0 ) THEN
               NIRRFLS = NIRRVAR / MXVARS3
            ELSE
               NIRRFLS = ( NIRRVAR / MXVARS3 ) + 1
            ENDIF

            DO NFL = 1, NIRRFLS
  
               IF( NFL .EQ. 1 ) THEN
                  IRRFNAM( NFL ) = CTM_IRR_1
               ELSEIF( NFL .EQ. 2 ) THEN
                  IRRFNAM( NFL ) = CTM_IRR_2
               ELSEIF( NFL .EQ. 3 ) THEN
                  IRRFNAM( NFL ) = CTM_IRR_3
               ELSEIF( NFL .GT. 3 ) THEN
                  XMSG = 'Maximum number of IRR output files exceeded'
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
               ENDIF

            ENDDO

         ENDIF

      ENDIF    ! End LFIRST

#ifdef parallel
      IF ( .NOT. PIO_RE_INIT( COLROW, PACOLS, PAROWS, PALEVS, NTHIK,
     &                        MY_PACOLS, MY_PAROWS, NPCOL, NPROW, NPCOL*NPROW,
     &                        MYPE, WFLG ) ) THEN
         XMSG = 'Failed to re-initialize parallel I/O library.'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         RETURN
      END IF
#endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Output Process analysis data and re-initialize arrays
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( LIPR ) THEN

         DO NFL = 1, NIPRFLS

            NVARINDX = MIN ( MXVARS3, NIPRVAR - (NFL - 1) * MXVARS3 )

            DO NOUT = 1, NVARINDX
               VARINDX = ( NFL - 1 ) * MXVARS3 + NOUT


               IF ( .NOT. WRITE3( IPRFNAM( NFL ), IPRNAME( VARINDX ),
     &                 JDATE, JTIME, DELC( :,:,:,VARINDX ) ) ) THEN 
                  XMSG = 'Could not write ' // IPRFNAM( NFL )
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF

            ENDDO

            WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &            'Timestep written to', IPRFNAM( NFL ),
     &            'for date and time', JDATE, JTIME
         ENDDO

         DO  NOUT = 1, NIPRVAR
            DO L = 1, PALEVS
               DO R = 1, MY_PAROWS
                  DO C = 1, MY_PACOLS
                     DELC( C,R,L,NOUT ) = 0.0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         DO NS = 1, NCSAVE
            ICG = SV2GRD( NS )
#ifdef parallel
            CALL SUBST_DATA_COPY ( CGRID, CSAV, ICG, NS )
#else
            DO L = PA_BEGLEV, PA_ENDLEV
               PL = L - PA_BEGLEV + 1
               DO R = MY_BEGROW, MY_ENDROW
                  PR = R - MY_BEGROW + 1
                  DO C = MY_BEGCOL, MY_ENDCOL
                     PC = C - MY_BEGCOL + 1
                     CSAV( PC,PR,PL,NS ) = CGRID( C,R,L,ICG ) 
                  ENDDO
               ENDDO
            ENDDO
#endif
         ENDDO

      ENDIF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Output the IRR data and re-initialize arrays
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF ( LIRR ) THEN

         ALLOCATE ( DBUFF( MY_PACOLS,MY_PAROWS,PALEVS ), STAT = ALLOCSTAT )

         DO NFL = 1, NIRRFLS

            NVARINDX = MIN (MXVARS3, NIRRVAR - (NFL - 1) * MXVARS3)

            DO NOUT = 1, NVARINDX
               VARINDX = ( NFL - 1 ) * MXVARS3 + NOUT

#ifdef parallel
               CALL SUBST_DATA_COPY ( IRROUT, TIRROUT, VARINDX )
               DBUFF = TIRROUT( 1:MY_PACOLS,1:MY_PAROWS,1:PALEVS )
#else
               DBUFF = IRROUT( 1:MY_PACOLS,1:MY_PAROWS,1:PALEVS,VARINDX )
#endif
               IF ( .NOT. WRITE3( IRRFNAM( NFL ), IRRNAME( VARINDX ),
     &              JDATE, JTIME, DBUFF ) ) THEN 
                  XMSG = 'Could not write ' // IRRFNAM( NFL )
                  CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               ENDIF

            ENDDO

            WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &            'Timestep written to', IRRFNAM( NFL ),
     &            'for date and time', JDATE, JTIME

         ENDDO

         DO NIRR = 1, NIRRVAR
#ifdef parallel
         DO L = MY_IRR_BEGLEV, MY_IRR_ENDLEV
            DO R = MY_IRR_BEGROW, MY_IRR_ENDROW
               DO C = MY_IRR_BEGCOL, MY_IRR_ENDCOL
#else
         DO L = 1, PALEVS
            DO R = 1, MY_PAROWS
               DO C = 1, MY_PACOLS
#endif
                      IRROUT( C,R,L,NIRR ) = 0.0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

      IF ( LIRR ) DEALLOCATE ( DBUFF )

      ENDIF

#ifdef parallel
      WFLG = .FALSE.
      IF ( .NOT. PIO_RE_INIT( COLROW, GL_NCOLS, GL_NROWS, NLAYS, NTHIK,
     &                        NCOLS, NROWS, NPCOL, NPROW, NPCOL*NPROW,
     &                        MYPE, WFLG ) ) THEN
         XMSG = 'Failed to re-initialize parallel I/O library.'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         RETURN
      END IF
#endif

      RETURN
      END
