       SUBROUTINE CHEM_BWD( CGRID, JDATE, JTIME, TSTEP )
      
      USE iso_c_binding
      USE HGRD_DEFN             ! horizontal domain specifications
      USE VGRD_DEFN             ! vertical layer specifications
      USE HRDATA
      USE KPP_Integrator_adj
      USE CHK_FILES_DEFS        ! adjoint checkpoint files definition
      USE CHK_FILES             ! adjoint checkpoint files routines

      IMPLICIT NONE 

C..Includes:

      INCLUDE SUBST_IOPARMS   ! Io/api parameters
      INCLUDE SUBST_IOFDESC   ! Io/api file descriptions
      INCLUDE SUBST_IODECL    ! Io/api declarations
      INCLUDE SUBST_FILES_ID  ! CMAQ files
      INCLUDE SUBST_CONST     ! CMAQ constants
      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs
      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block
      INCLUDE SUBST_GC_EMIS   ! Gas chem emissions name and mapping tables

#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! Emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! Emissions processing in vdif
#endif

      INCLUDE SUBST_PACTL_ID  ! Process analysis control parameters

C..Arguments:
      REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations

      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      INTEGER TSTEP( 2 )      ! Time step vector (HHMMSS)

C..Parameters:
C Pascal to atm conversion factor
      REAL, PARAMETER :: PA2ATM = 1.0 / STDATMPA

C..External Functions:
      INTEGER INDEX1          ! Looks up name in a list
      INTEGER JUNIT           ! Gets logical device number
      INTEGER SEC2TIME        ! Returns time interval from seconds
      INTEGER TIME2SEC        ! Returns seconds in time interval

      INTEGER, PARAMETER :: batch_size = 12288 * 16
      INTEGER :: batch_start, batch_end, batch_count, NG, total_grids

C..Saved Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'KPPDRIVER'     ! Program name

      INTEGER, SAVE :: ISTFL            ! Unit no. of iteration stat output file
      LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call to this subroutine

      REAL, SAVE :: MAOMV              ! Mol Wt of air over Mol Wt of water

C..Scratch Local Variables:
      CHARACTER( 132 ) :: MSG       ! Message text
      CHARACTER(  16 ) :: VNAME     ! Name of I/O API data variable
      
      INTEGER C, E, L, R, S  ! Loop indices
      INTEGER :: MY_COL1, MY_COL2, MY_ROW1, MY_ROW2
      INTEGER :: GLOBAL_R, GLOBAL_C
      INTEGER :: LOCAL_R, LOCAL_C

      INTEGER ITMSTEP         ! Chemistry integration interval (sec)   
      INTEGER MIDDATE         ! Date at time step midpoint
      INTEGER MIDTIME         ! Time at time step midpoint
      INTEGER NDARK           ! Number of layer 1 cells in darkness
      INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER SPC             ! Species loop index
      INTEGER STATUS          ! Status code
      INTEGER MAXCOUNT

  
      LOGICAL LSUNLIGHT       ! Flag for sunlight

      REAL ATMPRES            ! Cell pressure
      REAL CHEMSTEP           ! Chemistry integration interval (min)
      REAL H2O                ! Cell H2O mixing ratio (ppmV)
      REAL TEMP               ! Cell Temperature

      REAL PRES(    NCOLS, NROWS, NLAYS )        ! Cell pressure (Pa)
      REAL QV(      NCOLS, NROWS, NLAYS )        ! Cell water vapor (Kg/Kg air)
      REAL TA(      NCOLS, NROWS, NLAYS )        ! Cell temperature (K)
      REAL RJIN( NPHOTAB )                       ! J-values for a cell
      REAL RJ( NCOLS, NROWS, NLAYS, NPHOTAB )    ! J-values for each cell  

      INTEGER     GXOFF, GYOFF          ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

!slz      REAL*8, ALLOCATABLE, SAVE :: CHEMLAST( :, :, : )
      REAL, ALLOCATABLE :: CHKGRID (:,:,:,:)
      INTEGER ALLOCSTAT

      real batch_value_in(batch_size, N_GC_SPC+NPHOTAB-11)
      real batch_value_in_tr(N_GC_SPC+NPHOTAB-11, batch_size)

      real batch_lambda_in(batch_size, N_GC_SPC-16)
      real batch_lambda_in_tr(N_GC_SPC-16, batch_size)

      real batch_lambda_out(batch_size, N_GC_SPC-16)
      real batch_lambda_out_tr(N_GC_SPC-16, batch_size)

      ! real batch_value_out(batch_size, N_GC_SPC)
      ! real batch_value_out_tr(N_GC_SPC, batch_size)

      REAL forward_value(N_GC_SPC+NPHOTAB-11)
      REAL forward_lambda(N_GC_SPC-16), backward_lambda(N_GC_SPC-16)
      integer :: start_count, end_count, count_rate
      integer :: start_count_s, end_count_s
      real :: elapsed_time

      INTEGER IN_SPE
      PARAMETER (IN_SPE = 56)
      INTEGER IN_INDICES(IN_SPE)
      DATA IN_INDICES /1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
     & 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     & 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     & 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     & 41, 42, 43, 44, 45, 46, 47, 48, 50, 52,
     & 53, 55, 56, 67, 68, 71/

      INTEGER NO_SPE
      PARAMETER (NO_SPE = 16)
      INTEGER NO_INDICES(NO_SPE)
      DATA NO_INDICES /49, 51, 54, 57, 58, 59,
     & 60, 61, 62, 63, 64, 65, 66, 69, 70, 72/


      INTERFACE
      SUBROUTINE UP_LAMBDA (CHKGRID_IN,lambda_in,lambda_out,
     &         N_GC_SPC,NPHOTAB,NCOUNT,MYPE,NG) bind(C,name="UP_LAMBDA")
         REAL(4), dimension(N_GC_SPC+NPHOTAB-11, NCOUNT), INTENT(IN) :: CHKGRID_IN
         REAL(4), dimension(N_GC_SPC-16, NCOUNT), INTENT(IN) :: lambda_in
         REAL(4), dimension(N_GC_SPC-16, NCOUNT), INTENT(OUT) :: lambda_out
         INTEGER, value, INTENT(IN) :: N_GC_SPC, NPHOTAB, NCOUNT, MYPE, NG
      END SUBROUTINE UP_LAMBDA
      END INTERFACE

      call system_clock(count_rate=count_rate)
      call system_clock(count=start_count)  
      IF( N_GC_SPC .EQ. 0 ) RETURN
      
      MY_COL1 = COLSX_PE(1, MYPE + 1)
      MY_COL2 = COLSX_PE(2, MYPE + 1)
   
      MY_ROW1 = ROWSX_PE(1, MYPE + 1)
      MY_ROW2 = ROWSX_PE(2, MYPE + 1)

      IF( LFIRST ) THEN

         LOGDEV = INIT3( )

         IF( MECHNAME .NE. 'CB05CL_AE5_AQ' .AND.
     &       MECHNAME .NE. 'CB05CL_AE5CA_AQ' .AND.
     &       MECHNAME .NE. 'CB05CL_AE5ST_AQ' ) THEN
             MSG = 'This version of the EBI solver can only be used with'
     &            // ' the CB05CL_AE5 chemical mechanisms'
             CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
         ENDIF

         IF( INDEX( MECHNAME, 'AE' ) .NE. 0 ) THEN
           L_AE_VRSN = .TRUE.
         ELSE
           L_AE_VRSN = .FALSE.
         ENDIF

         IF( INDEX( MECHNAME, 'AQ' ) .NE. 0 ) THEN
           L_AQ_VRSN = .TRUE.
         ELSE
           L_AQ_VRSN = .FALSE.
         ENDIF

         IF( LIRR ) THEN
            MSG = 'IRR Analysis not allowed with EBI solver'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
         ENDIF 

         CALL HRINIT

         ITMSTEP = TIME2SEC( TSTEP( 2 ) )
         CHEMSTEP = FLOAT( ITMSTEP ) / 60.0
         WRITE( LOGDEV, 92000 ) CHEMSTEP

         WRITE( LOGDEV, 92020 )
         DO SPC = 1, N_GC_SPC
            WRITE( LOGDEV, 92040 ) GC_SPC( SPC ), RTOL( SPC )
         ENDDO

         MAOMV =  MWAIR / MWWAT

         total_grids = NLAYS * MY_NROWS * MY_NCOLS
c..If emissions processing requested stop
         IF( EMISCH ) THEN 

            MSG = 'ERROR: EBI solver not configured to '//
     &            'process emissions in chemistry'
            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )

         ENDIF   ! End if doing emissions


         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         LFIRST = .FALSE.
	 
	 
C Open checkpoint file as read-only
         CALL CHK_FILE_OPEN ( ADJ_CHEM_CHK, JDATE, JTIME, .FALSE. )

      ENDIF      ! First time
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] LFIRST', elapsed_time, ' seconds'

      call system_clock(count=start_count)  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Set date and time to center of time step, get necessary physical 
C  data, and get photolysis rates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      MIDDATE = JDATE
      MIDTIME = JTIME
      ITMSTEP = TIME2SEC( TSTEP( 2 ) )
      CHEMSTEP = FLOAT( ITMSTEP ) / 60.0D+00
      CALL NEXTIME( MIDDATE, MIDTIME, SEC2TIME( ITMSTEP / 2 ) )

C.. Get ambient temperature in K

      VNAME = 'TA' 
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, TA ) ) THEN
         MSG = 'Could not read TA from MET_CRO_3D'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF
      
C.. Get specific humidity in Kg H2O / Kg air
      VNAME = 'QV'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, QV ) ) THEN
         MSG = 'Could not read QV from MET_CRO_3D'
         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF 
      
C.. Get pressure in Pascals
      VNAME = 'PRES'
      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
     &                    MIDDATE, MIDTIME, PRES ) ) THEN
         MSG = 'Could not read PRES from MET_CRO_3D'
         CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT1 )
      ENDIF

      ALLOCATE ( CHKGRID( MY_NCOLS,MY_NROWS,NLAYS,N_GC_SPC+NPHOTAB+3 ), STAT = ALLOCSTAT )
!slz      ALLOCATE ( CHKGRID( MY_NCOLS,MY_NROWS,NLAYS,N_GC_SPC+1 ), STAT = ALLOCSTAT )

      IF ( ALLOCSTAT .NE. 0 ) THEN
         CALL M3EXIT ( PNAME, JDATE, JTIME, 'CHKGRID. allocation failed in kppdriver', XSTAT1 )
      END IF

C Read concentration checkpoint at sync step
      CALL CHK_FILE_READ ( ADJ_CHEM_CHK, JDATE, JTIME, CHKGRID )
      RJ(:,:,:,1:NPHOTAB) = CHKGRID(:,:,:,N_GC_SPC+1:N_GC_SPC+NPHOTAB)
      
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] GET CHKGRID VALUE', elapsed_time, ' seconds'   




      DO L = 1, NLAYS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                    DO S = 1, N_GC_SPC
                  if (CGRID( C, R, L, S ) .NE. 0.0) THEN
                        write(logdev, *) 'C,R,L,S',C,R,L,S
                        write(logdev,*) 'value:', CGRID( C, R, L, S )
                  endif
                    ENDDO
               ENDDO
            ENDDO
      ENDDO
      DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
      write( logdev, * ) 'before LAMBDA( 4 ) in first cell. ROW', R, 'C', C,':', CGRID( C, R, 1, 4 ) 
            enddo
      ENDDO 

      call system_clock(count=start_count)
      NCOUNT = 0
      NCOUNT2 = 0

      DO L = 1, NLAYS
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                    NCOUNT = NCOUNT + 1
               ENDDO
            ENDDO
      ENDDO

      write(logdev,*) 'the number of total grids are:', NCOUNT
      write(logdev,*) 'NLAYS, MY_NROWS, MY_NCOLS', NLAYS, MY_NROWS, MY_NCOLS     
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] COUNT CELLS', elapsed_time, ' seconds'

      call system_clock(count=start_count_s) 
      batch_count = (total_grids + batch_size - 1) / batch_size
      write(logdev,*) 'BATCH_SIZE and BATCH_COUNT', batch_size, batch_count

      DO NG = 1, batch_count
            batch_start = (NG - 1) * batch_size + 1
            batch_end = MIN(NG * batch_size, total_grids)

            NCOUNT = 0
            NCOUNT2 = 0

            call system_clock(count=start_count) 

            DO L = 1, NLAYS
                  DO R = 1, MY_NROWS
                      DO C = 1, MY_NCOLS
                        NCOUNT = NCOUNT + 1
                        LSUNLIGHT = .FALSE.
                        DO NPH = 1, NPHOTAB
                              RJIN( NPH ) = RJ( C, R, L, NPH )
                              IF( RJ( C, R, L, NPH ) .GT. 0.0 ) LSUNLIGHT = .TRUE.
                        ENDDO 
                          IF (NCOUNT >= batch_start .AND. NCOUNT <= batch_end) THEN
                  DO i = 1, IN_SPE
                  forward_value(i) = CHKGRID(C, R, L, in_indices(i))
                  forward_lambda(i) = CGRID( C, R, L, in_indices(i) )
                  END DO
                  forward_value(IN_SPE+1:IN_SPE+NPHOTAB+3) = CHKGRID(C,R,L,N_GC_SPC+1:N_GC_SPC+NPHOTAB+3)
                  forward_value(IN_SPE+NPHOTAB+4) = MERGE(1, 0, LSUNLIGHT)
                  forward_value(IN_SPE+NPHOTAB+5) = CHEMSTEP
                  batch_value_in(NCOUNT - batch_start + 1, :) = forward_value(:)
                  batch_lambda_in(NCOUNT - batch_start + 1, :) = forward_lambda(:)
                          END IF
                      ENDDO
                  ENDDO
            ENDDO

            call system_clock(count=end_count)
            elapsed_time = real(end_count - start_count) / real(count_rate)
            write(logdev,*) '[time check] LOAD ARRAY batch_value_in', elapsed_time, ' seconds' 
  
            call system_clock(count=start_count) 
            ! call cpp interace function to update lambda!
            batch_value_in_tr = TRANSPOSE(batch_value_in)
            batch_lambda_in_tr = TRANSPOSE(batch_lambda_in)
            write(logdev,*) 'before call update lambda'
            CALL UP_LAMBDA(batch_value_in_tr,batch_lambda_in_tr, batch_lambda_out_tr,
     &        N_GC_SPC, NPHOTAB, batch_end - batch_start + 1, MYPE, NG)
            write(logdev,*) 'after call update lambda'
            batch_lambda_out = TRANSPOSE(batch_lambda_out_tr)

            DO L = 1, NLAYS
                  DO R = 1, MY_NROWS
                      DO C = 1, MY_NCOLS
                        NCOUNT2 = NCOUNT2 + 1
                        IF (NCOUNT2 >= batch_start .AND. NCOUNT2 <= batch_end) THEN
                              backward_lambda(:) = batch_lambda_out(NCOUNT2 - batch_start + 1, :)

                  DO i = 1, IN_SPE
                  CGRID(C, R, L, in_indices(i)) = backward_lambda(i)
                  END DO

                  DO i = 1, NO_SPE
                  CGRID(C, R, L, NO_INDICES(i)) = 0.0
                  END DO

                        END IF
                  ENDDO
              ENDDO
            ENDDO
            call system_clock(count=end_count)
            elapsed_time = real(end_count - start_count) / real(count_rate)
            write(logdev,*) '[time check] Fortran-C binding', elapsed_time, ' seconds' 
        ENDDO
        call system_clock(count=end_count_s)
        elapsed_time = real(end_count_s - start_count_s) / real(count_rate)
        write(logdev,*) '[time check] UPDATE the CONC', elapsed_time, ' seconds'

	       
      DEALLOCATE (CHKGRID)

      RETURN

C*********************** FORMAT STATEMENTS ****************************

92000 FORMAT( / 10X, 'Euler Backward Iterative Parameters -'
     &        / 10X, 'Chemistry Integration Time Interval (min):', F12.4,
     &        / 10X, 'EBI maximum time step (min):              ', F12.4 )

92020 FORMAT( //10X, 'Species convergence tolerances:' )

92040 FORMAT(   10X, A16, 2X, 1PE12.2 )

      END