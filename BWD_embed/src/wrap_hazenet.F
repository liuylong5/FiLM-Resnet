C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C Function:
C   DESIGNED FOR THE PASS OF ARRAY IN JULIA
C   CREATED BY HAZENET@SJTU
      MODULE WRAP_SUBROUTINES
C-----------------------------------------------------------------------
      USE CGRID_SPCS
      USE VGRD_DEFN
      USE HGRD_DEFN
      USE PCGRID_DEFN
      USE ADJOINT_FILES

      INTEGER, ALLOCATABLE, SAVE :: GLOBAL_ASTEP(:) ! advection time steps
      LOGICAL, ALLOCATABLE, SAVE :: GLOBAL_XFIRST(:)
	  
	CONTAINS
C-----------------------------------------------------------------------
      SUBROUTINE GET_LGRID (JULIA_LGRID)
         REAL, INTENT(INOUT) :: JULIA_LGRID(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD) 
         INTEGER :: i, j, k, l
         DO l = 1, NSPCSD
            DO k = 1, NLAYS
               DO j = 1, NROWS
                  DO i = 1, NCOLS
                      JULIA_LGRID(i, j, k, l) = PLGRID(i, j, k, l)
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE GET_LGRID

C-----------------------------------------------------------------------
      SUBROUTINE SET_LGRID (JULIA_LGRID)
         REAL, INTENT(IN) :: JULIA_LGRID(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD)
         INTEGER :: i, j, k, l
         DO l = 1, NSPCSD
            DO k = 1, NLAYS
               DO j = 1, NROWS
                  DO i = 1, NCOLS
                      PLGRID(i, j, k, l) = JULIA_LGRID(i, j, k, l)
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE SET_LGRID

C-----------------------------------------------------------------------
      SUBROUTINE GET_LGRID_EM (JULIA_LGRID_EM)
         REAL, INTENT(INOUT) :: JULIA_LGRID_EM(1:NCOLS, 1:NROWS, 1:NLAYS, 1:N_SPC_PEMIS) 
         INTEGER :: i, j, k, l
         DO l = 1, N_SPC_PEMIS
            DO k = 1, NLAYS
               DO j = 1, NROWS
                  DO i = 1, NCOLS
                      JULIA_LGRID_EM(i, j, k, l) = PLGRID_EM(i, j, k, l)
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE GET_LGRID_EM

C-----------------------------------------------------------------------
      SUBROUTINE SET_LGRID_EM (JULIA_LGRID_EM)
         REAL, INTENT(IN) :: JULIA_LGRID_EM(1:NCOLS, 1:NROWS, 1:NLAYS, 1:N_SPC_PEMIS)
         INTEGER :: i, j, k, l
         DO l = 1, N_SPC_PEMIS
            DO k = 1, NLAYS
               DO j = 1, NROWS
                  DO i = 1, NCOLS
                      PLGRID_EM(i, j, k, l) = JULIA_LGRID_EM(i, j, k, l)
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE SET_LGRID_EM

C-----------------------------------------------------------------------
      SUBROUTINE JULIA_INIT (LGRID, LGRID_EM)
         REAL, POINTER :: LGRID(:, :, :, :)
         REAL, POINTER :: LGRID_EM(:, :, :, :)

         LGRID => PLGRID(1:NCOLS, 1:NROWS, :, :)
         LGRID_EM => PLGRID_EM(1:NCOLS, 1:NROWS, :, :)
         LGRID(:, :, :, :) = 0.0
      END SUBROUTINE JULIA_INIT
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_CGRID_MAP ()            
          INCLUDE SUBST_GC_SPC      ! gas chemistry species table
          INCLUDE SUBST_AE_SPC      ! aerosol species table
          INCLUDE SUBST_NR_SPC      ! non-reactive species table
          INCLUDE SUBST_TR_SPC      ! tracer species table

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
      END SUBROUTINE WRAP_CGRID_MAP
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_RD_FORCE_FILE(ARRAY, JDATE, JTIME, FRCFAC_CURRENT)
C Function:   DESIGNED FOR THE PASS OF ARRAY IN RD_FORCE_FILE
C Arguments:
          INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
          REAL, INTENT(IN) :: FRCFAC_CURRENT
          INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
          REAL :: ARRAY(:,:,:,:)

          FRCFAC = FRCFAC_CURRENT

          CALL RD_FORCE_FILE(ARRAY, JDATE, JTIME)    

      END SUBROUTINE WRAP_RD_FORCE_FILE
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_ADVSTEP(JDATE, JTIME, TSTEP)
C Function:   DESIGNED FOR THE PASS OF ARRAY IN ADVSTEP
C Arguments:
          INTEGER, INTENT(IN) :: JDATE ! current model simulation date (YYYYDDD)
          INTEGER, INTENT(IN) :: JTIME ! current model simulation time (HHMMSS)
          INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                              !   TSTEP(1) = local output step
                              !   TSTEP(2) = sciproc sync. step (chem)
          INTEGER :: ALLOCSTAT
          LOGICAL, SAVE :: FIRSTIME = .TRUE.

          INTERFACE
            SUBROUTINE ADVSTEP(JDATE, JTIME, TSTEP, ASTEP)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: JDATE ! current model simulation date (YYYYDDD)
            INTEGER, INTENT(IN) :: JTIME ! current model simulation time (HHMMSS)
            INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                    !   TSTEP(1) = local output step
                                    !   TSTEP(2) = sciproc sync. step (chem)
            INTEGER, INTENT(OUT) :: ASTEP(:) ! advection time steps
            END SUBROUTINE ADVSTEP
          END INTERFACE   

          IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            ALLOCATE(GLOBAL_ASTEP(1:NLAYS), STAT = ALLOCSTAT)
            IF ( ALLOCSTAT /= 0 ) THEN
                  WRITE(*,*) 'Failure allocating ASTEP'
            END IF
          END IF

          CALL ADVSTEP(JDATE, JTIME, TSTEP, GLOBAL_ASTEP)

      END SUBROUTINE WRAP_ADVSTEP
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_HADV_BWD(LGRID, JDATE, JTIME, TSTEP)
C Function:   DESIGNED FOR THE PASS OF ARRAY IN ADVSTEP
C Arguments:
          REAL, POINTER :: LGRID(:,:,:,:) ! concentrations
          INTEGER, INTENT(IN) :: JDATE ! current model simulation date (YYYYDDD)
          INTEGER, INTENT(IN) :: JTIME ! current model simulation time (HHMMSS)
          INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                              !   TSTEP(1) = local output step
                              !   TSTEP(2) = sciproc sync. step (chem)
          INTEGER :: ALLOCSTAT
          LOGICAL, SAVE :: FIRSTIME = .TRUE.

          INTERFACE
            SUBROUTINE HADV_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)
            IMPLICIT NONE
C Arguments:
            REAL, POINTER :: LGRID(:,:,:,:) ! concentrations
            INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD
            INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
            INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                    ! TSTEP(1) = local output step
                                    ! TSTEP(2) = sciproc sync. step (chem)
            INTEGER, INTENT(IN) :: ASTEP(:) ! advection time steps for layers (HHMMSS)
            LOGICAL, INTENT(IN OUT) :: XFIRST(:)
            END SUBROUTINE HADV_BWD
          END INTERFACE

          IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            ALLOCATE(GLOBAL_XFIRST(1:NLAYS), STAT = ALLOCSTAT)
            IF ( ALLOCSTAT /= 0 ) THEN
                  WRITE(*,*) 'Failure allocating XFIRST'
            END IF
          END IF

          CALL HADV_BWD(LGRID, JDATE, JTIME, TSTEP, GLOBAL_ASTEP, GLOBAL_XFIRST)

      END SUBROUTINE WRAP_HADV_BWD
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_WR_ADJ_FILE(JDATE, JTIME, ADJ_FILE)
C Arguments:
          INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD      
          INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
          CHARACTER(16), INTENT(IN) :: ADJ_FILE
          REAL, POINTER :: LGRID(:, :, :, :)
          REAL, POINTER :: LGRID_EM(:, :, :, :)                

          IF ( ADJ_FILE == ADJ_LGRID ) THEN
             LGRID => PLGRID(1:NCOLS, 1:NROWS, :, :)
             CALL WR_ADJ_FILE(LGRID, JDATE, JTIME, ADJ_FILE)    
          END IF

          IF ( ADJ_FILE == ADJ_LGRID_EM ) THEN
             LGRID_EM => PLGRID_EM(1:NCOLS, 1:NROWS, :, :)
             CALL WR_ADJ_FILE(LGRID_EM, JDATE, JTIME, ADJ_FILE)    
          END IF

      END SUBROUTINE WRAP_WR_ADJ_FILE
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_WR_LGRID_AVG(JULIA_ARRAY, JDATE, JTIME, ADJ_FILE)
C Arguments:
          INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD      
          INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
          CHARACTER(16), INTENT(IN) :: ADJ_FILE
          REAL, INTENT(IN) :: JULIA_ARRAY(1:NCOLS, 1:NROWS, 1:NLAYS, 1:NSPCSD)

          CALL WR_ADJ_FILE(JULIA_ARRAY, JDATE, JTIME, ADJ_FILE)  

      END SUBROUTINE WRAP_WR_LGRID_AVG
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_WR_LGRID_EM_AVG(JULIA_ARRAY, JDATE, JTIME, ADJ_FILE)
C Arguments:
          INTEGER, INTENT(IN) :: JDATE ! current model date, coded YYYYDDD      
          INTEGER, INTENT(IN) :: JTIME ! current model time, coded HHMMSS
          CHARACTER(16), INTENT(IN) :: ADJ_FILE
          REAL, INTENT(IN) :: JULIA_ARRAY(1:NCOLS, 1:NROWS, 1:NLAYS, 1:N_SPC_PEMIS) 

         CALL WR_ADJ_FILE(JULIA_ARRAY, JDATE, JTIME, ADJ_FILE)  

      END SUBROUTINE WRAP_WR_LGRID_EM_AVG
C-----------------------------------------------------------------------
	END MODULE WRAP_SUBROUTINES