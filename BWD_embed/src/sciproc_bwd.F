
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/driver/yamo_adj_bwd/sciproc_bwd.F,v 1.3 2013/05/15 20:14:11 ppercell Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCIPROC_BWD(LGRID, LGRID_EM, 
     &                       JDATE, JTIME, TSTEP, ASTEP, XFIRST)

C-----------------------------------------------------------------------
C Function:
C   Science process adjoint for sensitivity analysis

C Revision History:
C   Aug 2010 by Peter Percell at UH-IMAQS: created for adjoint/4dvar
C   Feb 2011 by Peter Percell at UH-IMAQS: added calls to DECOUPLE/COUPLE
C   May 2013 by Peter Percell at UH-IMAQS:
C   - call HADV_BWD instead of directly managing calls to XADV_BWD and YADV_BWD
C   - changed the calling argument NADVSS back to the traditional ASTEP
C   - changed XFIRST to an array indexed by layer
C-----------------------------------------------------------------------
      USE ADJOINT_FILES
      USE GRID_CONF, ONLY: MYPE
      IMPLICIT NONE

C Include Files:
      INCLUDE SUBST_IODECL      ! I/O definitions and declarations

C Arguments:
      REAL, POINTER :: LGRID(:, :, :, :) ! adjoint variable
      REAL, POINTER :: LGRID_EM(:, :, :, :) ! adjoint variable
      INTEGER, INTENT(IN OUT) :: JDATE ! current model date, coded YYYYDDD
      INTEGER, INTENT(IN OUT) :: JTIME ! current model time, coded HHMMSS
      INTEGER, INTENT(IN) :: TSTEP(2) ! time step vector (HHMMSS)
                                      !   TSTEP(1) = local output step
                                      !   TSTEP(2) = sciproc sync. step (chem)
      INTEGER, INTENT(IN) :: ASTEP(:) ! layer advection time step
      LOGICAL, INTENT(IN OUT) :: XFIRST(:) ! layer XFIRST

C Local Variables:
      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      CHARACTER(36), SAVE :: NMSG = 'After NEXTIME: returned JDATE, JTIME'

      ! added by Hazenet@SJTU
      integer:: time_init, time_begin, time_end, time_final,count_rate,count_max
      character(120) :: LOGFILE
      logical :: exist
      

      INTEGER :: SDATE          ! current science process date, coded YYYYDDD
      INTEGER :: STIME          ! current science process time, coded HHMMSS

      INTERFACE
         SUBROUTINE AERO_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE AERO_BWD
         SUBROUTINE CHEM_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE CHEM_BWD
         SUBROUTINE CLDPROC_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE CLDPROC_BWD
         SUBROUTINE DECOUPLE(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE DECOUPLE
         SUBROUTINE COUPLE_CADJ(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE COUPLE_CADJ
         SUBROUTINE DECOUPLE_CADJ(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE DECOUPLE_CADJ
         SUBROUTINE HDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE HDIFF_BWD
         SUBROUTINE ZADV_BWD(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE ZADV_BWD
         SUBROUTINE HADV_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)
            IMPLICIT NONE
            REAL, POINTER           :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN)     :: JDATE, JTIME
            INTEGER, INTENT(IN)     :: TSTEP(2)
            INTEGER, INTENT(IN)     :: ASTEP(:)
            LOGICAL, INTENT(IN OUT) :: XFIRST(:)
         END SUBROUTINE HADV_BWD
         SUBROUTINE COUPLE(LGRID, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE COUPLE
         SUBROUTINE VDIFF_BWD(LGRID, LGRID_EM, JDATE, JTIME, TSTEP)
            IMPLICIT NONE
            REAL, POINTER       :: LGRID(:, :, :, :)
            REAL, POINTER       :: LGRID_EM(:, :, :, :) ! adjoint variable
            INTEGER, INTENT(IN) :: JDATE, JTIME
            INTEGER, INTENT(IN) :: TSTEP(2)
         END SUBROUTINE VDIFF_BWD
      END INTERFACE

C-----------------------------------------------------------------------
      IF ( FIRSTIME ) THEN
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()
      END IF

      SDATE = JDATE
      STIME = JTIME

      CALL NEXTIME(JDATE, JTIME, -TSTEP(2))
      WRITE(LOGDEV, '(/ 5X, A, I8, I7.6)') NMSG, JDATE, JTIME

      IF ( MYPE == 0 ) THEN
         call getenv( 'MY_LOG_FILE', LOGFILE )
         inquire(file=LOGFILE, exist=exist)
         if (exist) then
            open(10, file=LOGFILE, status="old", position="append", action="write")
         else
            open(10, file=LOGFILE, status="new", action="write")
         end if
         write(10,*) ''
         CALL system_clock ( time_init,count_rate,count_max )
         write(10,*) 'The current sync timestep was', JDATE, JTIME
      END IF

      CALL RD_FORCE_FILE(LGRID, JDATE, JTIME)      
      
      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'RD_FORCE_FILE module used:', (time_end - time_init) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max  )
      ! END IF

      ! CALL AERO_BWD(LGRID, JDATE, JTIME, TSTEP)

      IF ( MYPE == 0 ) THEN
          CALL system_clock ( time_end,count_rate,count_max )
          write(10,*) 'AERO_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
          CALL system_clock ( time_begin,count_rate,count_max   )
      END IF

      CALL CHEM_BWD(LGRID, JDATE, JTIME, TSTEP)

      IF ( MYPE == 0 ) THEN
          CALL system_clock ( time_end,count_rate,count_max )
          write(10,*) 'CHEM_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
          CALL system_clock ( time_begin,count_rate,count_max   )
      END IF

      ! CALL CLDPROC_BWD(LGRID, JDATE, JTIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'CLDPROC_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL DECOUPLE(LGRID, SDATE, STIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'DECOUPLE module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL HDIFF_BWD(LGRID, JDATE, JTIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'HDIFF_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL COUPLE(LGRID, JDATE, JTIME, TSTEP)
      
      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'COUPLE module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL DECOUPLE_CADJ(LGRID, JDATE, JTIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'DECOUPLE_CADJ module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL HADV_BWD(LGRID, JDATE, JTIME, TSTEP, ASTEP, XFIRST)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'HADV_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF


      ! CALL ZADV_BWD(LGRID, JDATE, JTIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'ZADV_BWD module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF


      ! CALL COUPLE_CADJ(LGRID, JDATE, JTIME, TSTEP)

      ! IF ( MYPE == 0 ) THEN
      !     CALL system_clock ( time_end,count_rate,count_max )
      !     write(10,*) 'COUPLE_CADJ module used:', (time_end - time_begin) / count_rate, ' seconds'
      !     CALL system_clock ( time_begin,count_rate,count_max   )
      ! END IF

      ! CALL VDIFF_BWD(LGRID, LGRID_EM, JDATE, JTIME, TSTEP)

      IF ( MYPE == 0 ) THEN
          CALL system_clock ( time_final,count_rate,count_max )
          write(10,*) 'VDIFF_BWD module used:', (time_final - time_begin) / count_rate, ' seconds'
          write(10,*) 'The whole sysnc timestep used:', (time_final - time_init) / count_rate, ' seconds'
          write(10,*) ''
          close(10)
      END IF

      RETURN
      END






