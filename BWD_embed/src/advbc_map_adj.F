
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
C $Header: /Volumes/Data/CVS/CMAQ_CVSrepos/CCTM/src/hadv/yamo_cadj_fwd/advbc_map.F,v 1.1 2013/05/14 15:53:42 ppercell Exp $
 
C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%
 
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE ADVBC_MAP_ADJ ( CONCMIN, BCNAME, BCFAC, FZERO)
 
C-----------------------------------------------------------------------
C Function:
C   Set up advected species mapping to BCON data - read variable BCNAME
C   from file - if BCNAME is blank, skip read and zero out the BCON array
C   for the variable corresponding to the blank BCNAME
C   BCNAME returned is either the file variable name, the file name
C   corresponding to the ICBC surrogate nate, or blank

C Preconditions:

C Subroutines and functions called:

C   OPEN3, DESC3, INDEX1, TRIMLEN, FINDEX, M3MESG, M3EXIT

C Revision history:
C   Jeff - Aug 1997 Based on beta version

C-----------------------------------------------------------------------

      IMPLICIT NONE
 
C Includes:

      INCLUDE SUBST_GC_ADV    ! gas chem advection species and map table
      INCLUDE SUBST_AE_ADV    ! aerosol advection species and map table
      INCLUDE SUBST_NR_ADV    ! non-react advection species and map table
      INCLUDE SUBST_TR_ADV    ! tracer advection species and map table
      INCLUDE SUBST_GC_ICBC   ! gas chem ic/bc surrogate names and map table
      INCLUDE SUBST_AE_ICBC   ! aerosol ic/bc surrogate names and map table
      INCLUDE SUBST_NR_ICBC   ! non-react ic/bc surrogate names and map table
      INCLUDE SUBST_TR_ICBC   ! tracer ic/bc surrogate names and map table
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
      INCLUDE SUBST_IOFDESC   ! file header data structure
      INCLUDE SUBST_IODECL    ! I/O definitions and declarations
      INCLUDE SUBST_FILES_ID  ! file name parameters

C Arguments:

      CHARACTER( 16 ) :: CONCMIN         ! min value allowed for BCs
      CHARACTER( 16 ) :: BCNAME( * )     ! BC name for adv species
      REAL               BCFAC ( * )     ! Scale factor for BCs
      LOGICAL         :: FZERO           !slz if true, use zeros

C External Functions not previously declared in IODECL3.EXT:
 
      INTEGER, EXTERNAL :: FINDEX, INDEX1, TRIMLEN

C Local variables:
 
      CHARACTER( 16 ) :: PNAME = 'ADVBC_MAP_ADJ'
      CHARACTER( 16 ) :: BLNK = ' '
      CHARACTER( 96 ) :: XMSG = ' '
      CHARACTER( 16 ) :: WRVAR( N_GC_ADV + N_AE_ADV + N_NR_ADV + N_TR_ADV )

      INTEGER    WRCNT                ! no. of species not on file
      INTEGER :: XDATE = 999999       ! dummy date
      INTEGER :: XTIME = 999999       ! dummy time
      INTEGER    INDX                 ! Species index
      INTEGER    ISUR                 ! Surrogate index
      INTEGER    SPC, VAR             ! loop copunters
      INTEGER    LOGDEV

C Map advected gas chemistry species to bc's

      LOGDEV = INIT3 ()
 
      SPC = 0
 
      IF ( N_GC_ADV .GT. 0 ) THEN
         IF ( .NOT. OPEN3( BNDY_GASC_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // trim(BNDY_GASC_1) // ' file; use zero values'
!slzx            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            CALL M3WARN(PNAME, XDATE, XTIME, XMSG) !slz
            FZERO=.TRUE. !slz
            RETURN
            END IF

         WRITE( LOGDEV,* ) "Adjoint boundary conditions applied"
         WRITE( LOGDEV,* ) "Adjoint boundary conditions file: "//trim(BNDY_GASC_1)
 
         IF ( .NOT. DESC3( BNDY_GASC_1 ) ) THEN
            XMSG = 'Could not get ' // trim(BNDY_GASC_1) // ' file description'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF
 
         WRCNT = 0
         DO VAR = 1, N_GC_ADV
            SPC = SPC + 1
            BCFAC( SPC ) = 1.0
            INDX = INDEX1( GC_ADV( VAR ), NVARS3D, VNAME3D )
            IF ( INDX .EQ. 0 ) THEN      ! name not on file, look in ICBC table
                ISUR = FINDEX( GC_ADV_MAP( VAR ), N_GC_ICBC, GC_ICBC_MAP )
                IF ( ISUR .NE. 0 ) THEN  ! name found in ICBC table
                   INDX = INDEX1( GC_ICBC( ISUR), NVARS3D, VNAME3D )
                   BCFAC( SPC ) = GC_ICBC_FAC( ISUR )
                   END IF
               END IF
            IF ( INDX .NE. 0 ) THEN      ! name either on file or in table
               BCNAME( SPC ) = VNAME3D( INDX )
               ELSE
               BCNAME( SPC ) = BLNK
               WRCNT = WRCNT + 1
               WRVAR( WRCNT ) = GC_ADV( VAR )
               END IF
            END DO
         IF ( WRCNT .GT. 0 ) THEN
            WRITE( LOGDEV,1009 ) BNDY_GASC_1( 1:TRIMLEN( BNDY_GASC_1 ) ),
     &                           CONCMIN( 1:TRIMLEN( CONCMIN ) )
            DO VAR = 1, WRCNT
               WRITE( LOGDEV,1015 ) WRVAR( VAR )( 1:TRIMLEN( WRVAR( VAR ) ) )
               END DO
            END IF 
         END IF

C Map advected aerosol species to bc's

      IF ( N_AE_ADV .GT. 0 ) THEN

         IF ( .NOT. OPEN3( BNDY_AERO_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // BNDY_AERO_1 // ' file'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         IF ( .NOT. DESC3( BNDY_AERO_1 ) ) THEN
            XMSG = 'Could not get ' // BNDY_AERO_1 // ' file description'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         WRCNT = 0
         DO VAR = 1, N_AE_ADV
            SPC = SPC + 1
            BCFAC( SPC ) = 1.0
            INDX = INDEX1( AE_ADV( VAR ), NVARS3D, VNAME3D )
            IF ( INDX .EQ. 0 ) THEN
                ISUR = FINDEX( AE_ADV_MAP( VAR ), N_AE_ICBC, AE_ICBC_MAP )
                IF ( ISUR .NE. 0 ) THEN
                   INDX = INDEX1( AE_ICBC( ISUR ), NVARS3D, VNAME3D )
                   BCFAC( SPC ) = AE_ICBC_FAC( ISUR )
                   END IF
               END IF
            IF ( INDX .NE. 0 ) THEN
               BCNAME( SPC ) = VNAME3D( INDX ) 
               ELSE
               BCNAME( SPC ) = BLNK
               WRCNT = WRCNT + 1
               WRVAR( WRCNT ) = AE_ADV( VAR )
               END IF
            END DO
         IF ( WRCNT .GT. 0 ) THEN
            WRITE( LOGDEV,1009 ) BNDY_AERO_1( 1:TRIMLEN( BNDY_AERO_1 ) ),
     &                           CONCMIN( 1:TRIMLEN( CONCMIN ) )
            DO VAR = 1, WRCNT
               WRITE( LOGDEV,1015 ) WRVAR( VAR )( 1:TRIMLEN( WRVAR( VAR ) ) )
               END DO
            END IF 
         END IF
  
C Map advected non-reactive species to bc's

      IF ( N_NR_ADV .GT. 0 ) THEN

         IF ( .NOT. OPEN3( BNDY_NONR_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // BNDY_NONR_1 // ' file'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         IF ( .NOT. DESC3( BNDY_NONR_1 ) ) THEN
            XMSG = 'Could not get ' // BNDY_NONR_1 // ' file description'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         WRCNT = 0
         DO VAR = 1, N_NR_ADV
            SPC = SPC + 1
            BCFAC( SPC ) = 1.0
            INDX = INDEX1( NR_ADV( VAR ), NVARS3D, VNAME3D )
            IF ( INDX .EQ. 0 ) THEN
               ISUR = FINDEX( NR_ADV_MAP( VAR ), N_NR_ICBC, NR_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  INDX = INDEX1( NR_ICBC( ISUR ), NVARS3D, VNAME3D )
                  BCFAC( SPC ) = NR_ICBC_FAC( ISUR )
                  END IF
               END IF
            IF ( INDX .NE. 0 ) THEN
               BCNAME( SPC ) = VNAME3D( INDX ) 
               ELSE
               BCNAME( SPC ) = BLNK
               WRCNT = WRCNT + 1
               WRVAR( WRCNT ) = NR_ADV( VAR )
               END IF
            END DO
         IF ( WRCNT .GT. 0 ) THEN
            WRITE( LOGDEV,1009 ) BNDY_NONR_1( 1:TRIMLEN( BNDY_NONR_1 ) ),
     &                           CONCMIN( 1:TRIMLEN( CONCMIN ) )
            DO VAR = 1, WRCNT
               WRITE( LOGDEV,1015 ) WRVAR( VAR )( 1:TRIMLEN( WRVAR( VAR ) ) )
               END DO
            END IF 
         END IF
  
C Map advected tracer species to bc's

      IF ( N_TR_ADV .GT. 0 ) THEN

         IF ( .NOT. OPEN3( BNDY_TRAC_1, FSREAD3, PNAME ) ) THEN
            XMSG = 'Could not open ' // BNDY_TRAC_1 // ' file'
           CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         IF ( .NOT. DESC3( BNDY_TRAC_1 ) ) THEN
            XMSG = 'Could not get ' // BNDY_TRAC_1 // ' file description'
            CALL M3EXIT( PNAME, XDATE, XTIME, XMSG, XSTAT1 )
            END IF

         WRCNT = 0
         DO VAR = 1, N_TR_ADV
            SPC = SPC + 1
            BCFAC( SPC ) = 1.0
            INDX = INDEX1( TR_ADV ( VAR ) , NVARS3D, VNAME3D )
            IF ( INDX .EQ. 0 ) THEN
               ISUR = FINDEX( TR_ADV_MAP( VAR ), N_TR_ICBC, TR_ICBC_MAP )
               IF ( ISUR .NE. 0 ) THEN
                  INDX = INDEX1( TR_ICBC( ISUR ) , NVARS3D, VNAME3D )
                  BCFAC( SPC ) = TR_ICBC_FAC( ISUR )
                  END IF
               END IF
            IF ( INDX .NE. 0 ) THEN
               BCNAME( SPC ) = VNAME3D( INDX ) 
               ELSE
               BCNAME( SPC ) = BLNK
               WRCNT = WRCNT + 1
               WRVAR( WRCNT ) = TR_ADV( VAR )
               END IF
            END DO
         IF ( WRCNT .GT. 0 ) THEN
            WRITE( LOGDEV,1009 ) BNDY_TRAC_1( 1:TRIMLEN( BNDY_TRAC_1 ) ),
     &                           CONCMIN( 1:TRIMLEN( CONCMIN ) )
            DO VAR = 1, WRCNT
               WRITE( LOGDEV,1015 ) WRVAR( VAR )( 1:TRIMLEN( WRVAR( VAR ) ) )
               END DO
            END IF 
         END IF

      RETURN

1009  FORMAT( / 5X, 'No BC''s in file ', A, ' for the following adv species:',
     &          1X, 'Set to ', A )

1015  FORMAT( 10X, A )

      END
