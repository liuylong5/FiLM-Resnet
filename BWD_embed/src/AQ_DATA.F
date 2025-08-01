
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header$

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Module aq_data

C  Defines gas and aerosol surrogate data.
C  this module replace the AQ_PARAMS.EXT include file
C
C  Revision History:
C     First version was coded in May 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C
C-----------------------------------------------------------------------
      Use aero_data, Only: n_aerospc, conmin

      USE HGRD_DEFN

      Implicit None

      Integer, Parameter :: ngas = 15   ! number of gas-phase species for AQCHEM
      Integer, Parameter :: naddaer = 7 ! number of additional aerosol surrogate species for AQCHEM
      Integer, Parameter :: naer = n_aerospc + naddaer  ! number of aerosol surrogate species for AQCHEM
      Integer, Parameter :: nliqs = 44  ! number of liquid-phase species in AQCHEM
      Integer, Parameter :: nmodes = 3  ! number of aerosol modes
      Integer, Parameter :: mxsrg = 20  ! maximum number of species for surrogate

      integer, parameter :: akn = 1     ! number for Aitken mode
      Integer, Parameter :: acc = 2     ! number for Accumulation mode
      Integer, Parameter :: cor = 3     ! number for Coarse mode

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Gas surrogates species 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Type sur_gas_type
         Character( 16 ) :: name               ! surrogate name
         Real            :: backgnd            ! background value
      End Type sur_gas_type

      Type ( sur_gas_type ), Parameter :: srggas( ngas ) = (/
C                          Name       Background
C                       ----------    ----------
     & sur_gas_type( 'SO2       ',        0.0 ),
     & sur_gas_type( 'HNO3      ',        0.0 ),
     & sur_gas_type( 'N2O5      ',        0.0 ),
     & sur_gas_type( 'CO2       ',      340.0 ),
     & sur_gas_type( 'NH3       ',        0.0 ),
     & sur_gas_type( 'H2O2      ',        0.0 ),
     & sur_gas_type( 'O3        ',        0.0 ),
     & sur_gas_type( 'FOA       ',        0.0 ),
     & sur_gas_type( 'MHP       ',        0.0 ),
     & sur_gas_type( 'PAA       ',        0.0 ),
     & sur_gas_type( 'H2SO4     ',        0.0 ),
     & sur_gas_type( 'HCL       ',        0.0 ),
     & sur_gas_type( 'GLY       ',        0.0 ),
     & sur_gas_type( 'MGLY      ',        0.0 ),
     & sur_gas_type( 'HO        ',        0.0 ) / )

C...Required gas species
      Character( 16 ), Parameter, Private :: req_so2   = 'SO2'
      Character( 16 ), Parameter, Private :: req_hno3  = 'HNO3'
      Character( 16 ), Parameter, Private :: req_n2o5  = 'N2O5'
      Character( 16 ), Parameter, Private :: req_co2   = 'CO2'
      Character( 16 ), Parameter, Private :: req_nh3   = 'NH3'
      Character( 16 ), Parameter, Private :: req_h2o2  = 'H2O2'
      Character( 16 ), Parameter, Private :: req_o3    = 'O3'
      Character( 16 ), Parameter, Private :: req_foa   = 'FOA'
      Character( 16 ), Parameter, Private :: req_mhp   = 'MHP'
      Character( 16 ), Parameter, Private :: req_paa   = 'PAA'
      Character( 16 ), Parameter, Private :: req_h2so4 = 'H2SO4'
      Character( 16 ), Parameter, Private :: req_hcl   = 'HCL'  
      Character( 16 ), Parameter, Private :: req_gly   = 'GLY'  
      Character( 16 ), Parameter, Private :: req_mgly  = 'MGLY' 
      Character( 16 ), Parameter, Private :: req_ho    = 'HO'   

C...Pointers for the AQCHEM array GAS

      Integer :: lso2     ! Sulfur Dioxide
      Integer :: lhno3    ! Nitric Acid
      Integer :: ln2o5    ! Dinitrogen Pentoxide
      Integer :: lco2     ! Carbon Dioxide
      Integer :: lnh3     ! Ammonia
      Integer :: lh2o2    ! Hydrogen Perioxide
      Integer :: lo3      ! Ozone
      Integer :: lfoa     ! Formic Acid
      Integer :: lmhp     ! Methyl Hydrogen Peroxide
      Integer :: lpaa     ! Peroxyacidic Acid
      Integer :: lh2so4   ! Sulfuric Acid
      Integer :: lhcl     ! Hydrogen Chloride
      Integer :: lgly     ! Glyoxal
      Integer :: lmgly    ! Methylglyoxal
      Integer :: lho      ! OH, hydroxyl radical

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...Aerosol surrogates species 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Type sur_aer_type
         Character( 16 ) :: name( nmodes )     ! surrogate name
         Real            :: backgnd            ! background value
         Real            :: molwt              ! molecular weight
         Real            :: density            ! density (kg/m^3) 
      End Type sur_aer_type

      Type ( sur_aer_type ) :: srgaer( naer )

C... Additional Aerosol surrogate species to add to aerosol list
      Type ( sur_aer_type ), private :: addsrgaer( NADDAER ) = (/
C                         Aitken         Accum             Coarse       Background   MolWt    Density
C                      -------------   -------------   -------------    ---------- --------  ---------
     & sur_aer_type((/'             ','             ','CACO3        '/),     0.0,    100.1,   2000.0),
     & sur_aer_type((/'             ','             ','MGCO3        '/),     0.0,     84.3,   2000.0),
     & sur_aer_type((/'             ','             ','A3FE         '/),   0.010,     55.8,   2000.0),
     & sur_aer_type((/'             ','             ','B2MN         '/),   0.005,     54.9,   2000.0),
     & sur_aer_type((/'             ','             ','K            '/),     0.0,     39.1,   2000.0),
     & sur_aer_type((/'NUM_AITKEN   ','NUM_ACCUM    ','NUM_COARSE   '/),     0.0,      0.0,      1.0),
     & sur_aer_type((/'SRF_AITKEN   ','SRF_ACCUM    ','SRF_COARSE   '/),     0.0,      0.0,      1.0)/)

C...Required aerosol species
      Character(16), Parameter, Private :: req_aso4   = 'SO4_ACCUM'
      Character(16), Parameter, Private :: req_anh4   = 'NH4_ACCUM'
      Character(16), Parameter, Private :: req_ano3   = 'NO3_ACCUM'
      Character(16), Parameter, Private :: req_asoa   = 'SOA_ACCUM'
      Character(16), Parameter, Private :: req_aorgc  = 'AORGC_ACCUM'
!     Character(16), Parameter, Private :: req_apoa   = 'POA_ACCUM'
      Character(16), Parameter, Private :: req_apoa   = 'ORGP_ACCUM'
      Character(16), Parameter, Private :: req_aec    = 'EC_ACCUM'  
      Character(16), Parameter, Private :: req_apri   = 'PRI_ACCUM'
      Character(16), Parameter, Private :: req_apric  = 'PRI_COARSE'
      Character(16), Parameter, Private :: req_ana    = 'NA_ACCUM'  
      Character(16), Parameter, Private :: req_acl    = 'CL_ACCUM'  
      Character(16), Parameter, Private :: req_acaco3 = 'CACO3'
      Character(16), Parameter, Private :: req_amgco3 = 'MGCO3'
      Character(16), Parameter, Private :: req_aa3fe  = 'A3FE' 
      Character(16), Parameter, Private :: req_ab2mn  = 'B2MN' 
      Character(16), Parameter, Private :: req_ak     = 'K'   
      Character(16), Parameter, Private :: req_num    = 'NUM_ACCUM' 
      Character(16), Parameter, Private :: req_srf    = 'SRF_ACCUM' 

C...Pointers for the AQCHEM array AEROSOL

      Integer :: lso4     !  Sulfate
      Integer :: lnh4     !  Ammonium
      Integer :: lno3     !  Nitrate
      Integer :: lsoa     !  SOA
      Integer :: lorgc    !  SOA (anthropogenic organic oligomers)
      Integer :: lpoa     !  Primary organic aerosol
      Integer :: lec      !  Elemental carbon
      Integer :: lpri     !  Primary aerosol (Aitken and Accumulation)
      Integer :: lpricor  !  Primary aerosol (coarse)
      Integer :: lna      !  Sodium
      Integer :: lcl      !  Chloride ion
      Integer :: lcaco3   !  Calcium Carbonate aerosol (place holder)
      Integer :: lmgco3   !  Magnesium Carbonate aerosol (place holder)
      Integer :: la3fe    !  Iron aerosol (place holder)
      Integer :: lb2mn    !  Manganese aerosol (place holder)
      Integer :: lk       !  Potassium aerosol (Cl- tracked separately) (place holder)
      Integer :: lnum     !  Number
      Integer :: lsrf     !  Surface area

C...Pointers for the AQCHEM arrays LIQUID and WETDEP

      Integer, PARAMETER :: lacl        = 01  ! Hydrogen ion
      Integer, PARAMETER :: lnh4accl    = 02  ! Ammonium(attributed to accum. mode)
      Integer, Parameter :: lcal        = 03  ! Calcium
      Integer, Parameter :: lnaaccl     = 04  ! Sodium
      Integer, Parameter :: lohl        = 05  ! Hydroxyl radical
      Integer, Parameter :: lso4accl    = 06  ! Sulfate (attributed to accumulation mode)
      Integer, Parameter :: lhso4accl   = 07  ! bisulfate (attributed to accumulation mode)
      Integer, Parameter :: lso3l       = 08  ! sulfite
      Integer, Parameter :: lhso3l      = 09  ! bisulfite
      Integer, Parameter :: lso2l       = 10  ! sulfur dioxide
      Integer, Parameter :: lco3l       = 11  ! carbonate
      Integer, Parameter :: lhco3l      = 12  ! bicarbonate
      Integer, Parameter :: lco2l       = 13  ! carbon dioxide
      Integer, Parameter :: lno3accl    = 14  ! nitrate(attributed to accumulation mode)
      Integer, Parameter :: lnh3l       = 15  ! ammonia
      Integer, Parameter :: lclaccl     = 16  ! chloride ion (attributed to accumulation mode)
      Integer, Parameter :: lh2o2l      = 17  ! hydrogen peroxide
      Integer, Parameter :: lo3l        = 18  ! ozone
      Integer, Parameter :: lfel        = 19  ! iron
      Integer, Parameter :: lmnl        = 20  ! Manganese
      Integer, Parameter :: lal         = 21  ! generalized anion associated with iron
      Integer, Parameter :: lfoal       = 22  ! Formic acid
      Integer, Parameter :: lhco2l      = 23  ! HCOO- ion
      Integer, Parameter :: lmhpl       = 24  ! Methyl hydrogen peroxide
      Integer, Parameter :: lpaal       = 25  ! Peroxyacidic acid
      Integer, Parameter :: lhcll       = 26  ! Hydrogen chloride
      Integer, Parameter :: lpriml      = 27  ! primary aerosol
      Integer, Parameter :: lmgl        = 28  ! Magnesium
      Integer, Parameter :: lkl         = 29  ! potassium
      Integer, Parameter :: lbl         = 30  ! generalized anion associated with manganese
      Integer, Parameter :: lhno3l      = 31  ! nitric acid
      Integer, Parameter :: lprimcorl   = 32  ! coarse-mode primary aerosol
      Integer, Parameter :: lnumcorl    = 33  ! coarse-mode number
      Integer, Parameter :: lts6corl    = 34  ! sulfate (attributed to coarse mode)
      Integer, Parameter :: lnacorl     = 35  ! sodium (attributed to coarse mode)
      Integer, Parameter :: lclcorl     = 36  ! chloride ion (attributed to coarse mode)
      Integer, Parameter :: lno3corl    = 37  ! nitrate (attributed to coarse mode)
      Integer, Parameter :: lnh4corl    = 38  ! ammonium (attributed to coarse mode)
      Integer, Parameter :: lsoal       = 39  ! secondary organic aerosols
      Integer, Parameter :: lorgcl      = 40  ! in-cloud soa
      integer, parameter :: lpoal       = 41  ! primary organic aerosols
      Integer, Parameter :: lecl        = 42  ! elemental carbon
      Integer, Parameter :: lglyl       = 43  ! glyoxal
      Integer, Parameter :: lmglyl      = 44  ! methylglyoxal

      Integer    :: nsrggas( ngas )           ! # surrogates in CGRID for each gas
      Integer    :: lsrggas( ngas,mxsrg )

      Integer    :: nsrgaer( naer,nmodes )    ! # surrogates in CGRID for each aerosol
      Integer    :: lsrgaer( naer,nmodes,mxsrg )

      Integer, Allocatable :: cgrid2aq_map( : )   ! allocated in AQ_DATA_INIT routine

      Logical, Private, Save :: mapped = .False.

      Contains
 
C------------------------------------------------------------------------
      Subroutine aq_data_init()

C  Defines mapping from CGRID for species concentration and moments.
 
C  Revision History:
C     First version was coded in May 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

#ifndef Incl
      Use cgrid_spcs          ! CGRID mechanism species
#else
      Use cgrid_spcs          ! CGRID species number and offsets
#endif
      Use aero_data

      Implicit None

C Includes:
      Include SUBST_IOPARMS   ! I/O parameters definitions
#ifdef Incl
      Include SUBST_GC_SPC              ! gas chemistry species table
      Include SUBST_AE_SPC              ! aerosol species table
      Include SUBST_NR_SPC              ! non-reactive species table
      Include SUBST_TR_SPC              ! tracer species table
      Include SUBST_GC_G2AQ             ! gas to aqueous maping table
      Include SUBST_AE_A2AQ             ! aerosol to aqueous maping table
      Include SUBST_NR_N2AQ             ! non-reactive to aqueous maping table
      Include SUBST_TR_T2AQ             ! tracer to aqueous maping table
#endif

C Functions:
      Integer, External :: index1, indexn

C Local Variables:
      Character( 100 ) :: xmsg
#ifndef Incl
      Character( 16 ), Allocatable, Save :: aename( : )
      Character( 16 ), Allocatable, Save :: cgrid2aq( : )
#else
      Character( 16 ) :: aename( n_ae_a2aq )
      Character( 16 ) :: cgrid2aq( n_gc_g2aqd + n_ae_a2aq + n_nr_n2aq + n_tr_t2aq )
#endif
      Character( 10 ) :: bkground
      Integer :: igas
      Integer :: iaer
      Integer :: isrg
      Integer :: spc
      Integer :: mode
      Integer :: i
      Integer :: n
      Integer :: var
      Logical :: new
      Integer :: lsrg( mxsrg )
      Integer :: n_cgrid2aq
      Integer :: ios

      INTEGER, EXTERNAL :: INIT3
      INTEGER LOGDEV

C-----------------------------------------------------------------------

      LOGDEV = INIT3()

#ifndef Incl
      Allocate ( aename( n_ae_a2aq ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating aename'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If
#endif

C...Build array of AE names that are mapped to AE_A2AQ
      do spc = 1, n_ae_a2aq
         aename( spc ) = ae_spc( ae_a2aq_map( spc ) )
      end do

C...Set number of species in CGRID and used in Aqeuous module
      n_cgrid2aq = n_gc_g2aq + n_ae_a2aq + n_nr_n2aq + n_tr_t2aq

!     Allocate ( cgrid2aq( n_cgrid2aq ),
!    &           cgrid2aq_map( n_cgrid2aq ), stat = ios )
!     If ( ios .Ne. 0 ) Then
!        xmsg = '*** Error allocating cgrid2aq or cgrid2aq_map'
!        Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
!     End If

#ifndef Incl
      Allocate ( cgrid2aq( n_cgrid2aq ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating cgrid2aq'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If
#endif

      Allocate ( cgrid2aq_map( n_cgrid2aq ), stat = ios )
      If ( ios .Ne. 0 ) Then
         xmsg = '*** Error allocating cgrid2aq or cgrid2aq_map'
         Call m3exit( 'aq_data_init', 0, 0, xmsg, xstat3 )
      End If
#ifdef Incl
      Call cgrid_map ( nspcsd, gc_strt, ae_strt, nr_strt, tr_strt )
#endif

C...Map aero species
      Call map_aero()

      Do spc = 1, n_aerospc
         srgaer( spc ) = sur_aer_type( (/' ',' ',' '/), 0.0, 0.0, 0.0 )

         ! Skip H2O and SOIL species
         If ( spc .Eq. ah2o_idx ) Cycle
         If ( spc .Eq. asoil_idx ) Cycle

         Do mode = 1, n_mode

            n = index1( aerospc( spc )%name( mode ), n_ae_a2aq, aename )
            If ( n .Gt. 0 ) Then
               ! Verify that surrogate is new
               new = .True. 
               Do i = 1, spc
                  If ( ae_a2aq( n ) .Eq. srgaer( i )%name( mode ) ) Then
                     new = .False.; Exit
                  End If
               End Do
               If ( new ) Then
                  srgaer( spc )%name( mode ) = ae_a2aq( n )
                  srgaer( spc )%density = aerospc( spc )%density
                  srgaer( spc )%molwt = aerospc_mw( spc )
               End If
            End If

         End Do   ! mode

      End Do   ! spc

C...Add addition aerosol surrogate species
      Do n = 1, naddaer
         spc = n_aerospc + n
         srgaer( spc ) = addsrgaer( n )
      End Do

!  write out SRGAER to verify
!     Do spc = 1, naer
!        Write(*,'(3A,F6.0,2F8.2)') srgaer( spc )%name, srgaer( spc )%density,
!    &                              srgaer( spc )%molwt, srgaer( spc )%backgnd
!     End Do
      
C...Prepare indices for aqueous species

C...Load the CGRID to aqueous species pointers for the gases
      spc = 0
      Do var = 1, n_gc_g2aq
         spc = spc + 1
         cgrid2aq( spc ) = gc_g2aq( var )
         cgrid2aq_map( spc ) = gc_g2aq_map( var ) + gc_strt - 1
      End Do

C...Load the CGRID to aqueous species pointers for the aerosols
      Do var = 1, n_ae_a2aq
         spc = spc + 1
         cgrid2aq( spc ) = ae_a2aq( var )
         cgrid2aq_map( spc ) = ae_a2aq_map( var ) + ae_strt - 1
      End Do

C...Load the CGRID to aqueous species pointers for the non-reactives
      Do var = 1, n_nr_n2aq
         spc = spc + 1
         cgrid2aq( spc ) = nr_n2aq( var )
         cgrid2aq_map( spc ) = nr_n2aq_map( var ) + nr_strt - 1
      End Do

C...load the CGRID to aqueous species pointers for the tracers
      Do var = 1, n_tr_t2aq
         spc = spc + 1
         cgrid2aq( spc ) = tr_t2aq( var )
         cgrid2aq_map( spc ) = tr_t2aq_map( var ) + tr_strt - 1
      End Do

C...Create the pointers from CGRID to the gaseous species needed by AQCHEM
      nsrggas = 0
      Do igas = 1, ngas
         nsrggas( igas ) = indexn( srggas( igas )%name, n_cgrid2aq, cgrid2aq, lsrg )
         If ( nsrggas( igas ) .Eq. 0 ) Then
            Write ( bkground, '(F8.4)' ) srggas( igas )%backgnd
            xmsg = 'no surrogates for aqueous species '
     &           // trim( srggas( igas )%name )
     &           // '...using background value of '
     &           // trim( bkground )
            Call m3mesg( xmsg )
         Else
            Do isrg = 1, nsrggas( igas )
               lsrggas( igas,isrg ) = lsrg( isrg )
            End Do
         End If
      End Do

C...Create the pointers from CGRID to the aerosol species needed by AQCHEM
      nsrgaer = 0
      Do iaer = 1, naer
         Do n = 1, nmodes
            If ( srgaer( iaer )%name( n ) .Ne. ' ' ) Then
               nsrgaer( iaer,n ) = indexn( srgaer( iaer )%name( n ), n_cgrid2aq,
     &                                     cgrid2aq, lsrg )
               If ( nsrgaer( iaer,n ) .Eq. 0 ) Then
                  Write ( bkground, '(f8.4)' ) srgaer( iaer )%backgnd
                  xmsg = 'no surrogates for aqueous species '
     &                 // trim( srgaer( iaer )%name( n ) )
     &                 // '...using background value of '
     &                 // trim( bkground )
                  Call m3mesg( xmsg )
               Else
                  Do isrg = 1, nsrgaer( iaer,n )
                     lsrgaer( iaer,n,isrg ) = lsrg( isrg )
                  End Do
               End If
            End If
         End Do
      End Do

C...Create the pointers for required gas surrogate species needed by AQCHEM
      lso2   = findSgrGas( req_so2 )    ! Sulfur Dioxide
      lhno3  = findSgrGas( req_hno3 )   ! Nitric Acid
      ln2o5  = findSgrGas( req_n2o5 )   ! Dinitrogen Pentoxide
      lco2   = findSgrGas( req_co2 )    ! Carbon Dioxide
      lnh3   = findSgrGas( req_nh3 )    ! Ammonia
      lh2o2  = findSgrGas( req_h2o2 )   ! Hydrogen Perioxide
      lo3    = findSgrGas( req_o3 )     ! Ozone
      lfoa   = findSgrGas( req_foa )    ! Formic Acid
      lmhp   = findSgrGas( req_mhp )    ! Methyl Hydrogen Peroxide
      lpaa   = findSgrGas( req_paa )    ! Peroxyacidic Acid
      lh2so4 = findSgrGas( req_h2so4 )  ! Sulfuric Acid
      lhcl   = findSgrGas( req_hcl )    ! Hydrogen Chloride
      lgly   = findSgrGas( req_gly )    ! Glyoxal
      lmgly  = findSgrGas( req_mgly )   ! Methylglyoxal
      lho    = findSgrGas( req_ho )     ! OH, hydroxyl radical

C...Create the pointers for required aerosol surrogate species needed by AQCHEM
!     lso4   = findSgrAer( req_aso4 )    !  Sulfate
!     lnh4   = findSgrAer( req_anh4 )    !  Ammonium
!     lno3   = findSgrAer( req_ano3 )    !  Nitrate
!     lsoa   = findSgrAer( req_asoa )    !  SOA
!     lorgc  = findSgrAer( req_aorgc )   !  SOA (anthropogenic organic oligomers)
!     lpoa   = findSgrAer( req_apoa )    !  Primary organic aerosol
!     lec    = findSgrAer( req_aec   )   !  Elemental carbon
!     lpri   = findSgrAer( req_apri  )   !  Primary aerosol (Aitken and Accumulation)
!     lpricor= findSgrAer( req_apric )   !  Primary aerosol (coarse)
!     lna    = findSgrAer( req_ana )     !  Sodium
!     lcl    = findSgrAer( req_acl )     !  Chloride ion
!     lcaco3 = findSgrAer( req_acaco3 )  !  Calcium Carbonate aerosol (place holder)
!     lmgco3 = findSgrAer( req_amgco3 )  !  Magnesium Carbonate aerosol (place holder)
!     la3fe  = findSgrAer( req_aa3fe )   !  Iron aerosol (place holder)
!     lb2mn  = findSgrAer( req_ab2mn )   !  Manganese aerosol (place holder)
!     lk     = findSgrAer( req_ak  )     !  Potassium aerosol (Cl- tracked separately)
!     lnum   = findSgrAer( req_num )     !  Number
!     lsrf   = findSgrAer( req_srf )     !  Surface area
      lso4   = findSgrAer( logdev, req_aso4 )    !  Sulfate
      lnh4   = findSgrAer( logdev, req_anh4 )    !  Ammonium
      lno3   = findSgrAer( logdev, req_ano3 )    !  Nitrate
!      lsoa   = findSgrAer( logdev, req_asoa )    !  SOA
      lorgc  = findSgrAer( logdev, req_aorgc )   !  SOA (anthropogenic organic oligomers)
!      lpoa   = findSgrAer( logdev, req_apoa )    !  Primary organic aerosol
      lec    = findSgrAer( logdev, req_aec   )   !  Elemental carbon
      lpri   = findSgrAer( logdev, req_apri  )   !  Primary aerosol (Aitken and Accumulation)
      lpricor= findSgrAer( logdev, req_apric )   !  Primary aerosol (coarse)
      lna    = findSgrAer( logdev, req_ana )     !  Sodium
      lcl    = findSgrAer( logdev, req_acl )     !  Chloride ion
      lcaco3 = findSgrAer( logdev, req_acaco3 )  !  Calcium Carbonate aerosol (place holder)
      lmgco3 = findSgrAer( logdev, req_amgco3 )  !  Magnesium Carbonate aerosol (place holder)
      la3fe  = findSgrAer( logdev, req_aa3fe )   !  Iron aerosol (place holder)
      lb2mn  = findSgrAer( logdev, req_ab2mn )   !  Manganese aerosol (place holder)
      lk     = findSgrAer( logdev, req_ak  )     !  Potassium aerosol (Cl- tracked separately)
      lnum   = findSgrAer( logdev, req_num )     !  Number
      lsrf   = findSgrAer( logdev, req_srf )     !  Surface area

      mapped = .True.
      Return
      End Subroutine aq_data_init

C-----------------------------------------------------------------------
      Function findSgrGas( vname ) Result ( idx )

C  Finds the index of 'required' gas surrogate species in the srggas list

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Implicit None

C Arguments:
      Character( * ) :: vname
      Integer idx

C Includes:
      Include SUBST_IOPARMS   ! I/O parameters definitions

C Local Variables:
      Character( 80 ) :: xmsg
      Integer spc

      Do spc = 1, ngas      
        If ( srggas( spc )%name .Eq. Trim( vname ) ) Then
           idx = spc
           ! Write( *,'(a,i5)' ) vname, idx
           Return
        End If
      End Do

      xmsg = 'Required Species ' // Trim( vname ) // ' Not found in SgrGas names array'
      Call m3exit( 'AQ_DATA', 0, 0, xmsg, xstat3 )

      Return
      End Function findSgrGas

C-----------------------------------------------------------------------
!     Function findSgrAer( vname ) Result ( idx )
      Function findSgrAer( logdev, vname ) Result ( idx )

C  Finds the index of 'required' Aerosol surrogate species in the srgaer list

C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Implicit None

C Arguments:
      Character( * ) :: vname
      Integer idx
      INTEGER LOGDEV

C Includes:
      Include SUBST_IOPARMS   ! I/O parameters definitions

C Local Variables:
      Character( 80 ) :: xmsg
      Integer spc, n

      Do spc = 1, naer
        Do n = 1, nmodes

          WRITE( LOGDEV,* ) SPC, N, VNAME, SRGAER( SPC )%NAME( N )

          If ( srgaer( spc )%name( n ) .Eq. Trim( vname ) ) Then
             idx = spc
             ! Write( *,'(a,i5)' ) vname, idx
             Return
          End If
        End Do
      End Do

      xmsg = 'Required Species ' // Trim( vname ) // ' Not found in SgrAer names array'
      Call m3exit( 'AQ_DATA', 0, 0, xmsg, xstat3 )

      Return
      End Function findSgrAer

      End Module aq_data
