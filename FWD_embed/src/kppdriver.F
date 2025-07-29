       SUBROUTINE CHEM( CGRID, JDATE, JTIME, TSTEP )

C**********************************************************************
C
C     KPP implementaion based on MEBI backbone (Amir: 11/2010)
C
C**********************************************************************

      USE iso_c_binding
      USE HGRD_DEFN             ! horizontal domain specifications
      USE VGRD_DEFN             ! vertical layer specifications
      USE HRDATA
      USE KPP_Integrator
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
      INCLUDE SUBST_AE_SPC    ! Gas chem species names and MWs


#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! Emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! Emissions processing in vdif
#endif

      INCLUDE SUBST_PACTL_ID  ! Process analysis control parameters

C..Arguments:
!     REAL CGRID( NCOLS, NROWS, NLAYS, * )  ! Species concentrations
      REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations
!      REAL, POINTER :: CHKGRID( :, :, :, : )

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
      LOGICAL, SAVE :: ONEFIRST = .TRUE.
      LOGICAL, SAVE :: TWOFIRST = .TRUE.

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
!      INTEGER VAR             ! Variable number on I/O API file
  
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
      
      real batch_value_in(batch_size, N_GC_SPC+NPHOTAB-11)
      real batch_value_in_tr(N_GC_SPC+NPHOTAB-11, batch_size)

      real batch_value_out(batch_size, N_GC_SPC-6)
      real batch_value_out_tr(N_GC_SPC-6, batch_size)

      INTEGER     GXOFF, GYOFF          ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3

      REAL*8, ALLOCATABLE, SAVE :: CHEMLAST( :, :, : )
      REAL, ALLOCATABLE :: CHKGRID (:,:,:,:)
      INTEGER ALLOCSTAT

CCCCCCCC
      REAL forward_value(N_GC_SPC+NPHOTAB-11), backward_value(N_GC_SPC-6)
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

      INTEGER OUT_SPE
      PARAMETER (OUT_SPE = 66)
      INTEGER OUT_INDICES(OUT_SPE)
      DATA OUT_INDICES /1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
     & 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     & 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     & 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     & 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
     & 51, 52, 53, 54, 55, 56, 63, 64, 65, 66,
     & 67, 68, 69, 70, 71, 72/

      INTERFACE
      SUBROUTINE DNN_CONC (CHKGRID_IN, CHKGRID_OUT, N_GC_SPC, NPHOTAB, NCOUNT, MYPE, NG) bind(C, name="DNN_CONC")
         REAL(4), dimension(N_GC_SPC+NPHOTAB-11, NCOUNT), INTENT(IN) :: CHKGRID_IN
         REAL(4), dimension(N_GC_SPC-6, NCOUNT), INTENT(OUT) :: CHKGRID_OUT
         INTEGER, value, INTENT(IN) :: N_GC_SPC, NPHOTAB, NCOUNT, MYPE, NG
      END SUBROUTINE DNN_CONC
      END INTERFACE

C**********************************************************************
      call system_clock(count_rate=count_rate)
      call system_clock(count=start_count)  
      IF( N_GC_SPC .EQ. 0 ) RETURN

      MY_COL1 = COLSX_PE(1, MYPE + 1)
      MY_COL2 = COLSX_PE(2, MYPE + 1)
   
      MY_ROW1 = ROWSX_PE(1, MYPE + 1)
      MY_ROW2 = ROWSX_PE(2, MYPE + 1)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  On first call, call routines to set-up for EBI solver and 
c  set-up to do emissions here if that option is invoked
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
	 
         ALLOCATE ( CHEMLAST( MY_NCOLS,MY_NROWS,NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            CALL M3EXIT ( PNAME, JDATE, JTIME, 'CHEMLAST 
     &                     allocation failed in kppdriver', XSTAT1 )
            END IF

         CHEMLAST = 1.0d0

      CALL CHK_FILE_OPEN ( ADJ_CHEM_START, JDATE, JTIME, .TRUE., TSTEP )
      ! CALL CHK_FILE_OPEN ( ADJ_CHEM_END, JDATE, JTIME, .TRUE., TSTEP )

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
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] get necessary physical data', elapsed_time, ' seconds'


      call system_clock(count=start_count)   
C.. Get photolysis rates in /min
      CALL PHOT ( MIDDATE, MIDTIME, JDATE, JTIME, NDARK, RJ )   
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] CALL PHOT( NOW NOT COMMENT OUT PHOT ) ', elapsed_time, ' seconds'

      call system_clock(count=start_count)      

C -------------->>
      ALLOCATE ( CHKGRID( MY_NCOLS,MY_NROWS,NLAYS,N_GC_SPC+NPHOTAB+3 ), STAT = ALLOCSTAT )

      IF ( ALLOCSTAT .NE. 0 ) THEN
         CALL M3EXIT ( PNAME, JDATE, JTIME, 'CHEMGRID. allocation failed in kppdriver', XSTAT1 )
      END IF
      CHKGRID(:,:,:,1:N_GC_SPC) = MAX (CGRID(:,:,:,1:N_GC_SPC), 1.0E-19)
      CHKGRID(:,:,:,N_GC_SPC+1:N_GC_SPC+NPHOTAB) = RJ(:,:,:,1:NPHOTAB)
      CHKGRID(:,:,:,N_GC_SPC+NPHOTAB+1) =  TA( :,:,: )
      CHKGRID(:,:,:,N_GC_SPC+NPHOTAB+2) =  PA2ATM * PRES( :,:,: )
      CHKGRID(:,:,:,N_GC_SPC+NPHOTAB+3) = QV( :,:,: ) * MAOMV *  1.0E+06
      ! IF( CALL_FILE_OPS ) THEN
            CALL CHK_FILE_WRITE ( ADJ_CHEM_START, JDATE, JTIME, 
     &            CHKGRID(:, :, :, 1:N_GC_SPC+NPHOTAB+3) )
      ! END IF
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] ALLOCATE  CHKGRID', elapsed_time, ' seconds'    

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
                  END DO
                  forward_value(IN_SPE+1:IN_SPE+NPHOTAB+3) = CHKGRID(C,R,L,N_GC_SPC+1:N_GC_SPC+NPHOTAB+3)
                  forward_value(IN_SPE+NPHOTAB+4) = MERGE(1, 0, LSUNLIGHT)
                  forward_value(IN_SPE+NPHOTAB+5) = CHEMSTEP
                  batch_value_in(NCOUNT - batch_start + 1, :) = forward_value(:)
                              IF ((ONEFIRST) .and. (NG == 1) .and. (NCOUNT == 1)) THEN
                                    WRITE(logdev, *) 'MYPE', MYPE
                                    WRITE(logdev, *) '[value check] JTIME and NCOUNT', JTIME, NCOUNT
                                    WRITE(logdev, *) '[value check] init_conc(4)', CHKGRID(1,1,1,4)
                                    WRITE(logdev, *) '[value check] batch_value_in(1, 4)', batch_value_in(1, 4)
                              END IF
                              IF ((ONEFIRST) .and. (NG == 1) .and. (NCOUNT == 2)) THEN
                                    WRITE(logdev, *) 'MYPE', MYPE
                                    WRITE(logdev, *) '[value check] JTIME and NCOUNT', JTIME, NCOUNT
                                    WRITE(logdev, *) '[value check] init_conc(4)', CHKGRID(2,1,1,4)
                                    WRITE(logdev, *) '[value check] batch_value_in(2, 4)', batch_value_in(2, 4)
                              END IF
                        END IF
                      ENDDO
                  ENDDO
              ENDDO

          call system_clock(count=end_count)
          elapsed_time = real(end_count - start_count) / real(count_rate)
          write(logdev,*) '[time check] LOAD ARRAY batch_value_in', elapsed_time, ' seconds' 

          call system_clock(count=start_count) 
          ! interface call
          batch_value_in_tr = TRANSPOSE(batch_value_in)
          CALL DNN_CONC(batch_value_in_tr, batch_value_out_tr, N_GC_SPC, NPHOTAB, batch_end - batch_start + 1, MYPE, NG)
          batch_value_out = TRANSPOSE(batch_value_out_tr)


          DO L = 1, NLAYS
            DO R = 1, MY_NROWS
                DO C = 1, MY_NCOLS
                  NCOUNT2 = NCOUNT2 + 1
                  IF (NCOUNT2 >= batch_start .AND. NCOUNT2 <= batch_end) THEN
                        backward_value(:) = batch_value_out(NCOUNT2 - batch_start + 1, :)
                        
                        DO S = 1, N_GC_SPC 
                              IF (S .GE. 57 .AND. S .LE. 62) THEN
                                    CGRID( C, R, L, S ) = 1.0E-19
                              ELSE
                                    DO i = 1, OUT_SPE
                                    CGRID( C, R, L, out_indices(i) ) = MAX( backward_value(i), 1.0E-19 )
                                    END DO
                              ENDIF
                        ENDDO   
                      IF ((TWOFIRST) .and. (NG == 1) .and. (NCOUNT2 == 1)) THEN
                        WRITE(logdev, *) 'MYPE', MYPE
                        WRITE(logdev, *) '[value check] JTIME and NCOUNT2', JTIME, NCOUNT2
                        WRITE(logdev, *) '[value check] init_conc(4)', CHKGRID(1,1,1,4)
                        WRITE(logdev, *) '[value check] batch_value_out(1, 4)', batch_value_out(1, 4)
                        WRITE(logdev, *) '[value check] CGRID(1,1,1,4)', CGRID(1,1,1,4)
                        write(logdev, *) '[value check] batch_value_in(1, 4)', batch_value_in(4, 1)
                      END IF
                      IF ((TWOFIRST) .and. (NG == 1) .and. (NCOUNT2 == 2)) THEN
                        WRITE(logdev, *) 'MYPE', MYPE
                        WRITE(logdev, *) '[value check] JTIME and NCOUNT2', JTIME, NCOUNT2
                        WRITE(logdev, *) '[value check] init_conc(4)', CHKGRID(2,1,1,4)
                        WRITE(logdev, *) '[value check] batch_value_out(2, 4)', batch_value_out(2, 4), batch_value_out(4, 2)
                        WRITE(logdev, *) '[value check] CGRID(2,1,1,4)', CGRID(2,1,1,4)
                        write(logdev, *) '[value check] batch_value_in(2, 4)',  batch_value_in(4, 2)
                      END IF
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



      CHKGRID(:,:,:,1:N_GC_SPC)=CGRID(:,:,:,1:N_GC_SPC)
      call system_clock(count=start_count) 
!       CALL CHK_FILE_WRITE ( ADJ_CHEM_END, JDATE, JTIME, 
!      &     CHKGRID(:, :, :, 1:N_GC_SPC) )

      DEALLOCATE (CHKGRID)
      call system_clock(count=end_count)
      elapsed_time = real(end_count - start_count) / real(count_rate)
      write(logdev,*) '[time check] CALL CHK_FILE_WRITE', elapsed_time, ' seconds'

      RETURN

C*********************** FORMAT STATEMENTS ****************************

92000 FORMAT( / 10X, 'Euler Backward Iterative Parameters -'
     &        / 10X, 'Chemistry Integration Time Interval (min):', F12.4,
     &        / 10X, 'EBI maximum time step (min):              ', F12.4 )

92020 FORMAT( //10X, 'Species convergence tolerances:' )

92040 FORMAT(   10X, A16, 2X, 1PE12.2 )

      END
