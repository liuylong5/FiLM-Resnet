C        Generated by TAPENADE     (INRIA, Tropics team)
C  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
C
C  Differentiation of tri in reverse (adjoint) mode:
C   gradient     of useful results: x
C   with respect to varying inputs: x b
C   RW status of diff variables: x:in-out b:out
C
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
C
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/vdiff/acm2_inline_adj_bwd/tri_bwd.F,v 1.1 2011/02/22 01:10:45 ppercell Exp $ 
C
C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)tridiag.F	1.1 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.tridiag.F 03 Jun 1997 11:45:09
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE TRI_BWD(l, d, u, bb, xb)

C Revision History:
C   Oct 2010 by Peter Percell at UH-IMAQS: created for adjoint/4dvar

      USE VGRD_DEFN           ! vertical layer specifications

      IMPLICIT NONE
C
C Includes:
C
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table

      INTEGER     N_SPC_DIFF               ! global diffusion species
      PARAMETER ( N_SPC_DIFF = N_GC_DIFF   ! 
     &                       + N_AE_DIFF
     &                       + N_NR_DIFF
     &                       + N_TR_DIFF)
C
C Arguments:
C
C subdiagonal
      REAL l(nlays)
C diagonal
      REAL d(nlays)
C superdiagonal
      REAL u(nlays)
C R.H. side
      REAL bb(n_spc_diff, nlays)
C solution
      REAL xb(n_spc_diff, nlays)
C
C Local Variables:
C
      REAL gam(nlays)
      REAL bet
      INTEGER v, k
      INTEGER ii2
      INTEGER ii1
C
C Decomposition and forward substitution:
      bet = 1.0/d(1)
C
      DO k=2,nlays
        gam(k) = bet*u(k-1)
        CALL PUSHREAL4(bet)
        bet = 1.0/(d(k)-l(k)*gam(k))
      ENDDO
      DO k=1,nlays-1,1
        DO v=n_spc_diff,1,-1
          xb(v, k+1) = xb(v, k+1) - gam(k+1)*xb(v, k)
        ENDDO
      ENDDO
      DO ii1=1,nlays
        DO ii2=1,n_spc_diff
          bb(ii2, ii1) = 0.0
        ENDDO
      ENDDO
      DO k=nlays,2,-1
        DO v=n_spc_diff,1,-1
          bb(v, k) = bb(v, k) + bet*xb(v, k)
          xb(v, k-1) = xb(v, k-1) - bet*l(k)*xb(v, k)
          xb(v, k) = 0.0
        ENDDO
        CALL POPREAL4(bet)
      ENDDO
      DO v=n_spc_diff,1,-1
        bb(v, 1) = bb(v, 1) + bet*xb(v, 1)
        xb(v, 1) = 0.0
      ENDDO
      END
