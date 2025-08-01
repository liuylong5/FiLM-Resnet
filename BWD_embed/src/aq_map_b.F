!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.6 (r4343) - 10 Feb 2012 10:52
!
!  Differentiation of aq_map in reverse (adjoint) mode:
!   gradient     of useful results: cend
!   with respect to varying inputs: polc alfa0 alfa2 alfa3 cend
!   RW status of diff variables: polc:out alfa0:out alfa2:out alfa3:out
!                cend:in-out
!
!***********************************************************************
!   Portions of Models-3/CMAQ software were developed or based on      *
!   information from various groups: Federal Government employees,     *
!   contractors working on a United States Government contract, and    *
!   non-Federal sources (including research institutions).  These      *
!   research institutions have given the Government permission to      *
!   use, prepare derivative works, and distribute copies of their      *
!   work in Models-3/CMAQ to the public and to permit others to do     *
!   so.  EPA therefore grants similar permissions for use of the       *
!   Models-3/CMAQ software, but users are requested to provide copies  *
!   of derivative works to the Government without restrictions as to   *
!   use by others.  Users are responsible for acquiring their own      *
!   copies of commercial software associated with Models-3/CMAQ and    *
!   for complying with vendor requirements.  Software copyrights by    *
!   the MCNC Environmental Modeling Center are used with their         *
!   permissions subject to the above restrictions.                     *
!***********************************************************************
!
! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header$
!
! what(1) key, module and SID; SCCS file; date and time of last delta:
! %W% %P% %G% %U%
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AQ_MAP_B(jdate, jtime, wtbar, wcbar, tbarc, pbarc, cthk1,
     +  airm, prate1, taucld, polc, polcb, cend, cendb, remov, removac, alfa0
     +  , alfa0b, alfa2, alfa2b, alfa3, alfa3b)
        USE AERO_DATA_B
        USE CGRID_SPCS
        USE AQ_DATA
        IMPLICIT NONE
#ifdef Incl
      INCLUDE SUBST_RXCMMN         ! Mechanism reaction common block (kf)
        INCLUDE SUBST_GC_SPC
        INCLUDE SUBST_AE_SPC
        INCLUDE SUBST_NR_SPC
        INCLUDE SUBST_TR_SPC
        INCLUDE SUBST_GC_G2AQ
        INCLUDE SUBST_AE_A2AQ
        INCLUDE SUBST_NR_N2AQ
        INCLUDE SUBST_TR_T2AQ 
#endif
        INCLUDE SUBST_IOPARMS
!
! Exit status message
! Number of species in CGRID
        CHARACTER(len=120) :: xmsg
!
! Parameters:
!
!
        INTEGER, PARAMETER :: mxspcs=n_gc_spcd+n_ae_spc+n_nr_spc+n_tr_spc
!
!
        REAL, PARAMETER :: onethird=1.0/3.0
        REAL, PARAMETER :: twothirds=2.0/3.0
!
! Arguments:
!
! current model date, coded YYYYDDD
        INTEGER, INTENT(IN) :: jdate
! current model time, coded HHMMSS
        INTEGER, INTENT(IN) :: jtime
!
! total air mass (moles/m2) in cloudy air
        REAL, INTENT(IN) :: airm
! scav coef for aitken aerosol number
        REAL, INTENT(IN) :: alfa0
        REAL :: alfa0b
! scav coef for aitken aerosol sfc area
        REAL, INTENT(IN) :: alfa2
        REAL :: alfa2b
! scav coef for aitken aerosol mass
        REAL, INTENT(IN) :: alfa3
        REAL :: alfa3b
! aq chem calc cloud thickness
        REAL, INTENT(IN) :: cthk1
! mean cloud pressure (Pa)
        REAL, INTENT(IN) :: pbarc
! storm rainfall rate (mm/hr)
        REAL, INTENT(IN) :: prate1
! variable storing H+ deposition
        REAL, INTENT(INOUT) :: removac
! cloud lifetime
        REAL, INTENT(IN) :: taucld
! mean cloud temp (K)
        REAL, INTENT(IN) :: tbarc
! liq water content of cloud (kg/m3)
        REAL, INTENT(IN) :: wcbar
! total wat cont (kg/m2) int. thru cld
        REAL, INTENT(IN) :: wtbar
!
!
! ave vert conc incloud (moles/mole)
        REAL, INTENT(IN) :: polc(mxspcs)
        REAL :: polcb(mxspcs)
! ending incloud conc (moles/mole)
        REAL, INTENT(INOUT) :: cend(mxspcs)
        REAL, INTENT(INOUT) :: cendb(mxspcs)
! moles/m2 or mm*mol/lit scavenged
        REAL, INTENT(INOUT) :: remov(mxspcs)
!
!
! Local Variables (scalars):
!
! flag for first pass thru
        LOGICAL, SAVE :: firstime=.true.

        CHARACTER( 32 ) :: MECHCHK
!
! program name
        CHARACTER(len=16), SAVE :: pname='aq_map_b'
!
! aerosol loop counter
        INTEGER :: iaer
! aerosol mode loop counter
        INTEGER :: imode
! gas loop counter
        INTEGER :: igas
! surrogate loop counter
        INTEGER :: isrg
! relative pointer variable
        INTEGER :: pntr
! liquid species loop counter
        INTEGER :: spc
!
        REAL*8 :: betaso4
        REAL*8 :: betaso4b
! EXP( -ALFA2 * TAUCLD )
        REAL :: ealfa2t
        REAL :: ealfa2tb
! modal mass at time t
        REAL*8 :: m3new(nmodes)
        REAL*8 :: m3newb(nmodes)
! modal mass at time 0
        REAL*8 :: m3old(nmodes)
        REAL*8 :: m3oldb(nmodes)
!
! hydrogen wet dep (mm mol/liter)
        REAL :: hpwdep
! gas phase conc (mol/mol)
        REAL*8 :: gas(ngas)
        REAL*8 :: gasb(ngas)
! gas phase wet dep array (mm mol/liter)
        REAL*8 :: gaswdep(ngas)
! aerosol conc (mol/mol)
        REAL*8 :: aerosol(naer, nmodes)
        REAL*8 :: aerosolb(naer, nmodes)
! aerosol wet dep array (mm mol/liter)
        REAL*8 :: aerwdep(naer, nmodes)
! weights for surrogate
        REAL*8 :: wsrggas(ngas, mxsrg)
        REAL*8 :: wsrggasb(ngas, mxsrg)
! weights for surrogate
        REAL*8 :: wsrgaer(naer, nmodes, mxsrg)
        REAL*8 :: wsrgaerb(naer, nmodes, mxsrg)
!
! External Functions:
!
! external func to get species pointers
        INTEGER, EXTERNAL :: indexn
!slz        EXTERNAL SIZE
        INTEGER :: res
        REAL :: tmp
        INTEGER :: ad_to
        INTEGER :: ad_to0
        INTEGER :: ad_to1
        INTEGER :: branch
        INTEGER :: ad_to2
        INTEGER :: ad_to3
        INTEGER :: ad_to4
        INTEGER :: ad_to5
        INTEGER :: ad_to6
        REAL*8 :: temp3
        REAL :: temp2
        REAL*8 :: temp1
        REAL*8 :: temp0
        INTRINSIC EXP
        REAL*8 :: tempb2
        REAL*8 :: tempb1
        REAL*8 :: tempb0
        INTRINSIC MAX
        REAL :: tmpb
        REAL*8 :: temp3b
!slz        INTEGER :: SIZE
        REAL*8 :: tempb
        REAL :: temp2b
        REAL*8 :: temp3b0
        REAL*8 :: temp
        REAL*8 :: temp4
!
!-----------------------------------------------------------------------
!
!...Initialization
!...  event-statistics variables.
!
        IF (firstime) then
         firstime=.false.
!slzdbg         CALL AQ_DATA_INIT()
!slzdbg aq_data_init called in aq_map; might better add firstime in routine aq_data_init
        end if
!
!
!
!
!...for subsequent calls, check to make sure some surrogates were
!...  specified, otherwise there is no need to perform aqueous chemistry
!
        res = SIZE(cgrid2aq_map)
        IF (res .NE. 0) THEN
!
!...load gas-phase concentrations
!
          DO igas=1,ngas
!
            gas(igas) = 0.0d0
!
            DO isrg=1,nsrggas(igas)
              pntr = cgrid2aq_map(lsrggas(igas, isrg))
              gas(igas) = gas(igas) + REAL(polc(pntr), 8)
            END DO
            CALL PUSHINTEGER4(isrg - 1)
!
            IF (gas(igas) .GT. 0.0d0) THEN
              DO isrg=1,nsrggas(igas)
                pntr = cgrid2aq_map(lsrggas(igas, isrg))
                wsrggas(igas, isrg) = REAL(polc(pntr), 8)/gas(igas)
              END DO
              CALL PUSHINTEGER4(isrg - 1)
              CALL PUSHCONTROL1B(0)
            ELSE
              DO isrg=1,nsrggas(igas)
                wsrggas(igas, isrg) = 1.0d0/REAL(nsrggas(igas), 8)
              END DO
              CALL PUSHINTEGER4(isrg - 1)
              CALL PUSHCONTROL1B(1)
            END IF
!
!...set background values for gases if no surrogates were specified
!
            IF (nsrggas(igas) .EQ. 0) THEN
              CALL PUSHREAL8(gas(igas))
              gas(igas) = REAL(srggas(igas)%backgnd, 8)*1.0d-6
              CALL PUSHCONTROL1B(1)
            ELSE
              CALL PUSHCONTROL1B(0)
            END IF
          END DO
!
!
!...load aerosol concentrations
!
          DO iaer=1,naer
            DO imode=1,nmodes
!
              IF (srgaer(iaer)%name(imode) .NE. ' ') THEN
                aerosol(iaer, imode) = 0.0d0
!
                DO isrg=1,nsrgaer(iaer, imode)
                  pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                  aerosol(iaer, imode) = aerosol(iaer, imode) + REAL(polc(pntr
     +             ), 8)
                END DO
                CALL PUSHINTEGER4(isrg - 1)
!
                IF (aerosol(iaer, imode) .GT. 0.0d0) THEN
                  DO isrg=1,nsrgaer(iaer, imode)
                    pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                    wsrgaer(iaer, imode, isrg) = REAL(polc(pntr), 8)/REAL(
     +               aerosol(iaer, imode), 8)
                  END DO
                  CALL PUSHINTEGER4(isrg - 1)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  DO isrg=1,nsrgaer(iaer, imode)
                    wsrgaer(iaer, imode, isrg) = 1.0d0/REAL(nsrgaer(iaer,
     +               imode), 8)
                  END DO
                  CALL PUSHINTEGER4(isrg - 1)
                  CALL PUSHCONTROL1B(1)
                END IF
!
!...set background values for aerosols if no surrogates were specified
!
                IF (nsrgaer(iaer, imode) .EQ. 0) THEN
                  IF (srgaer(iaer)%molwt .GT. 0.0) THEN
                    CALL PUSHREAL8(aerosol(iaer, imode))
                    aerosol(iaer, imode) = REAL(srgaer(iaer)%backgnd*1.0e-6*
     +               cthk1, 8)/REAL(srgaer(iaer)%molwt*airm, 8)
                    CALL PUSHCONTROL2B(3)
                  ELSE
                    CALL PUSHREAL8(aerosol(iaer, imode))
                    aerosol(iaer, imode) = REAL(srgaer(iaer)%backgnd*cthk1/
     +               airm, 8)
                    CALL PUSHCONTROL2B(2)
                  END IF
                ELSE
                  CALL PUSHCONTROL2B(1)
                END IF
              ELSE
                CALL PUSHCONTROL2B(0)
              END IF
            END DO
          END DO
          CALL PUSHINTEGER4ARRAY(aerospc_map, n_aerospc*n_mode)
!slztpn          CALL PUSHUNKNOWNTYPE(mapped)
          CALL PUSHINTEGER4(ae_strt)
          CALL EXTRACT_AERO(polc)
!
! *** extract grid cell concentrations of aero species from CGRID
!     into aerospc_conc in aero_data module
!
!
! *** Calculate pseudo aerosol 3rd moment (ignore factors that cancel in the division)
!
          m3old = 0.0
          DO imode=2,n_mode
!
            DO spc=1,n_aerospc
              IF (aerospc(spc)%name(imode) .NE. ' ' .AND. (.NOT.aerospc(spc)%
     +            iswet)) THEN
                m3old(imode) = m3old(imode) + aerospc_conc(spc, imode)*
     +            aerospc_mw(spc)/aerospc(spc)%density
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          CALL PUSHREAL8ARRAY(aerosol, naer*3)
          CALL PUSHREAL8ARRAY(gas, 15)
!
!...perform aqueous-phase chemistry calculations
          MECHCHK = MECHNAME  !kf
!
          CALL AQCHEM(jdate, jtime, tbarc, pbarc, taucld, prate1, wcbar,
     +             wtbar, airm, alfa0, alfa2, alfa3, gas, aerosol, gaswdep,
     +             aerwdep, hpwdep, betaso4, mechchk)
!
!...  compute the scavenging coefficient
          ealfa2t = EXP(-(alfa2*taucld))
!
!...store the amount of hydrogen deposition
!
!
!...Now, re-apportion mass back into cend/remov (cgrid-type) array
!
          DO igas=1,ngas
            DO isrg=1,nsrggas(igas)
              pntr = cgrid2aq_map(lsrggas(igas, isrg))
              cend(pntr) = gas(igas)*wsrggas(igas, isrg)
            END DO
            CALL PUSHINTEGER4(isrg - 1)
          END DO
!
          DO iaer=1,naer
            DO imode=1,nmodes
              IF (srgaer(iaer)%name(imode) .NE. ' ') THEN
                DO isrg=1,nsrgaer(iaer, imode)
                  pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                  cend(pntr) = aerosol(iaer, imode)*wsrgaer(iaer, imode, isrg)
                END DO
                CALL PUSHINTEGER4(isrg - 1)
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
          CALL PUSHINTEGER4ARRAY(aerospc_map, n_aerospc*n_mode)
!slztpn          CALL PUSHUNKNOWNTYPE(mapped)
          CALL PUSHREAL4ARRAY(aerospc_mw, n_aerospc)
          CALL PUSHINTEGER4(ae_strt)
          CALL EXTRACT_AERO(cend)
!
! *** extract grid cell concentrations of aero species from CGRID
!     into aerospc_conc in aero_data module
!
!
! *** Calculate pseudo aerosol 3rd moment (ignore factors that cancel in the division)
!
          m3new = 0.0
          DO imode=2,n_mode
!
            DO spc=1,n_aerospc
              IF (aerospc(spc)%name(imode) .NE. ' ' .AND. (.NOT.aerospc(spc)%
     +            iswet)) THEN
                m3new(imode) = m3new(imode) + aerospc_conc(spc, imode)*
     +            aerospc_mw(spc)/aerospc(spc)%density
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
          END DO
!
!...check for minimums
!
          DO imode=2,n_mode
            IF (m3old(imode) .LT. conmin) THEN
              m3old(imode) = conmin
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
              m3old(imode) = m3old(imode)
            END IF
            IF (m3new(imode) .LT. conmin) THEN
              m3new(imode) = conmin
              CALL PUSHCONTROL1B(0)
            ELSE
              CALL PUSHCONTROL1B(1)
              m3new(imode) = m3new(imode)
            END IF
          END DO
          CALL PUSHREAL4(cend(aerosrf_map(1)))
!
!...Update surface area
!
          cend(aerosrf_map(1)) = polc(aerosrf_map(1))*ealfa2t
          CALL PUSHREAL4(cend(aerosrf_map(2)))
          cend(aerosrf_map(2)) = polc(aerosrf_map(2))*EXP(-(betaso4*taucld*
     +      onethird))*(m3new(2)/m3old(2))**twothirds
!
!
!slzdbg          polcb = 0.0
          m3oldb = 0.0_8
          m3newb = 0.0_8
          tmpb = cendb(aerosrf_map(3))
          cendb(aerosrf_map(3)) = 0.0
!slzdbg          polcb = 0.0
          m3oldb = 0.0_8
          m3newb = 0.0_8
          temp3 = m3new(3)/m3old(3)
          temp4 = temp3**twothirds
          temp2 = cend(aeronum_map(3))/polc(aeronum_map(3))
          IF (temp2 .LE. 0.0 .AND. (onethird .EQ. 0.0 .OR. onethird .NE. INT(
     +        onethird))) THEN
            temp2b = 0.0
          ELSE
            temp2b = onethird*temp2**(onethird-1)*polc(aerosrf_map(3))*temp4*
     +        tmpb/polc(aeronum_map(3))
          END IF
          temp3b = temp2**onethird*tmpb
          IF (temp3 .LE. 0.0 .AND. (twothirds .EQ. 0.0 .OR. twothirds .NE. INT
     +        (twothirds))) THEN
            temp3b0 = 0.0
          ELSE
            temp3b0 = twothirds*temp3**(twothirds-1)*polc(aerosrf_map(3))*
     +        temp3b/m3old(3)
          END IF
          cendb(aeronum_map(3)) = cendb(aeronum_map(3)) + temp2b
          polcb(aeronum_map(3)) = polcb(aeronum_map(3)) - temp2*temp2b
          polcb(aerosrf_map(3)) = polcb(aerosrf_map(3)) + temp4*temp3b
          m3newb(3) = m3newb(3) + temp3b0
          m3oldb(3) = m3oldb(3) - temp3*temp3b0
          CALL POPREAL4(cend(aerosrf_map(2)))
          temp = m3new(2)/m3old(2)
          temp1 = temp**twothirds
          temp0 = -(taucld*onethird*betaso4)
          tempb1 = EXP(temp0)*cendb(aerosrf_map(2))
          IF (temp .LE. 0.0 .AND. (twothirds .EQ. 0.0 .OR. twothirds .NE. INT(
     +        twothirds))) THEN
            tempb2 = 0.0
          ELSE
            tempb2 = twothirds*temp**(twothirds-1)*polc(aerosrf_map(2))*tempb1
     +        /m3old(2)
          END IF
          betaso4b = -(EXP(temp0)*polc(aerosrf_map(2))*temp1*taucld*onethird*
     +      cendb(aerosrf_map(2)))
          polcb(aerosrf_map(2)) = polcb(aerosrf_map(2)) + temp1*tempb1
          m3newb(2) = m3newb(2) + tempb2
          m3oldb(2) = m3oldb(2) - temp*tempb2
          cendb(aerosrf_map(2)) = 0.0
          CALL POPREAL4(cend(aerosrf_map(1)))
          polcb(aerosrf_map(1)) = polcb(aerosrf_map(1)) + ealfa2t*cendb(
     +      aerosrf_map(1))
          ealfa2tb = polc(aerosrf_map(1))*cendb(aerosrf_map(1))
          cendb(aerosrf_map(1)) = 0.0
          DO imode=n_mode,2,-1
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) m3newb(imode) = 0.0_8
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) m3oldb(imode) = 0.0_8
          END DO
          aerospc_concb = 0.0
          DO imode=n_mode,2,-1
            DO spc=n_aerospc,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) aerospc_concb(spc, imode) = aerospc_concb(spc
     +            , imode) + aerospc_mw(spc)*m3newb(imode)/aerospc(spc)%
     +            density
            END DO
          END DO
          CALL POPINTEGER4(ae_strt)
          CALL POPREAL4ARRAY(aerospc_mw, n_aerospc)
!slztpn          CALL POPUNKNOWNTYPE(mapped)
          CALL POPINTEGER4ARRAY(aerospc_map, n_aerospc*n_mode)
          moment2_concb = 0.0
          moment0_concb = 0.0
          CALL EXTRACT_AERO_B(cend, cendb)
          wsrgaerb = 0.0_8
          aerosolb = 0.0_8
          DO iaer=naer,1,-1
            DO imode=nmodes,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) THEN
                CALL POPINTEGER4(ad_to6)
                DO isrg=ad_to6,1,-1
                  pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                  aerosolb(iaer, imode) = aerosolb(iaer, imode) + wsrgaer(iaer
     +              , imode, isrg)*cendb(pntr)
                  wsrgaerb(iaer, imode, isrg) = wsrgaerb(iaer, imode, isrg) +
     +              aerosol(iaer, imode)*cendb(pntr)
                  cendb(pntr) = 0.0
                END DO
              END IF
            END DO
          END DO
          gasb = 0.0_8
          wsrggasb = 0.0_8
          DO igas=ngas,1,-1
            CALL POPINTEGER4(ad_to5)
            DO isrg=ad_to5,1,-1
              pntr = cgrid2aq_map(lsrggas(igas, isrg))
              gasb(igas) = gasb(igas) + wsrggas(igas, isrg)*cendb(pntr)
              wsrggasb(igas, isrg) = wsrggasb(igas, isrg) + gas(igas)*cendb(
     +          pntr)
              cendb(pntr) = 0.0
            END DO
          END DO
          alfa2b = alfa2b - EXP(-(taucld*alfa2))*taucld*ealfa2tb
          CALL POPREAL8ARRAY(gas, 15)
          CALL POPREAL8ARRAY(aerosol, naer*3)
          CALL AQCHEM_B(jdate, jtime, tbarc, pbarc, taucld, prate1, wcbar,
     +            wtbar, airm, alfa0, alfa0b, alfa2, alfa3, alfa3b, gas, gasb
     +            , aerosol, aerosolb, gaswdep, aerwdep, hpwdep, betaso4,
     +            betaso4b, mechchk)
          aerospc_concb = 0.0
          DO imode=n_mode,2,-1
            DO spc=n_aerospc,1,-1
              CALL POPCONTROL1B(branch)
              IF (branch .NE. 0) aerospc_concb(spc, imode) = aerospc_concb(spc
     +            , imode) + aerospc_mw(spc)*m3oldb(imode)/aerospc(spc)%
     +            density
            END DO
          END DO
          CALL POPINTEGER4(ae_strt)
!slztpn          CALL POPUNKNOWNTYPE(mapped)
          CALL POPINTEGER4ARRAY(aerospc_map, n_aerospc*n_mode)
          moment2_concb = 0.0
          moment0_concb = 0.0
          CALL EXTRACT_AERO_B(polc, polcb)
          DO iaer=naer,1,-1
            DO imode=nmodes,1,-1
              CALL POPCONTROL2B(branch)
              IF (branch .LT. 2) THEN
                IF (branch .EQ. 0) GOTO 100
              ELSE IF (branch .EQ. 2) THEN
                CALL POPREAL8(aerosol(iaer, imode))
                aerosolb(iaer, imode) = 0.0_8
              ELSE
                CALL POPREAL8(aerosol(iaer, imode))
                aerosolb(iaer, imode) = 0.0_8
              END IF
              CALL POPCONTROL1B(branch)
              IF (branch .EQ. 0) THEN
                CALL POPINTEGER4(ad_to3)
                DO isrg=ad_to3,1,-1
                  pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                  tempb0 = REAL(aerosol(iaer, imode), 8)
                  polcb(pntr) = polcb(pntr) + wsrgaerb(iaer, imode, isrg)/tempb0
                  aerosolb(iaer, imode) = aerosolb(iaer, imode) - REAL(polc(
     +             pntr), 8)*wsrgaerb(iaer, imode, isrg)/tempb0**2
                  wsrgaerb(iaer, imode, isrg) = 0.0_8
                END DO
              ELSE
                CALL POPINTEGER4(ad_to4)
                DO isrg=ad_to4,1,-1
                  wsrgaerb(iaer, imode, isrg) = 0.0_8
                END DO
              END IF
              CALL POPINTEGER4(ad_to2)
              DO isrg=ad_to2,1,-1
                pntr = cgrid2aq_map(lsrgaer(iaer, imode, isrg))
                polcb(pntr) = polcb(pntr) + aerosolb(iaer, imode)
              END DO
              aerosolb(iaer, imode) = 0.0_8
100         CONTINUE
            END DO
          END DO
          DO igas=ngas,1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) THEN
              CALL POPREAL8(gas(igas))
              gasb(igas) = 0.0_8
            END IF
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPINTEGER4(ad_to0)
              DO isrg=ad_to0,1,-1
                pntr = cgrid2aq_map(lsrggas(igas, isrg))
                tempb = wsrggasb(igas, isrg)/gas(igas)
                polcb(pntr) = polcb(pntr) + tempb
                gasb(igas) = gasb(igas) - REAL(polc(pntr), 8)*tempb/gas(igas)
                wsrggasb(igas, isrg) = 0.0_8
              END DO
            ELSE
              CALL POPINTEGER4(ad_to1)
              DO isrg=ad_to1,1,-1
                wsrggasb(igas, isrg) = 0.0_8
              END DO
            END IF
            CALL POPINTEGER4(ad_to)
            DO isrg=ad_to,1,-1
              pntr = cgrid2aq_map(lsrggas(igas, isrg))
              polcb(pntr) = polcb(pntr) + gasb(igas)
            END DO
            gasb(igas) = 0.0_8
          END DO
        END IF
      END SUBROUTINE AQ_MAP_B
