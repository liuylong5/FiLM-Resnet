C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Module aero_data

C  Defines aerosol species arrays and the parameters required in aerosol
C  processing.
 
C  Contains:
C     Subroutine map_aero
C     Subroutine extract_aero
C     Subroutine update_aero
C     Function findAero
 
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C----------------------------------------------------------------------

      Implicit None

C Number of aerosol species and modes

      Integer, Parameter :: n_aerospc = 30      ! number of aero species
      Integer, Parameter :: n_mode = 3          ! number of modes:
                                                ! 1 = Aitken 
                                                ! 2 = accumulation 
                                                ! 3 = coarse 

C Default minimum concentration
      Real,    Parameter :: conmin = 1.0E-30    ! [ ug/m^3 ]

C Flag to obtain coagulation coefficients
C by analytical approximation (True) or by Gauss-Hermite quadrature (False)
!slzdbg      Logical, Parameter :: fastcoag_flag = .True.
      Logical, Parameter :: fastcoag_flag = .False.

C-------------------------------------------------------------------------------------------------------

      Type spcs_type
         Character( 16 ) :: name( n_mode )     ! names of aerosol species for each mode
         Real            :: min_conc( n_mode ) ! minimum concentration values for each mode
         Real            :: density            ! density [ kg/m^3 ]
         Logical         :: isWet              ! wet flag
         Logical         :: volatile           ! volatile flag
         Integer         :: charge             ! electroneutrality charge
         Real            :: visual_idx         ! visual index factor
         Character( 16 ) :: emis               ! emissions names
      End Type spcs_type

      Type( spcs_type ), Parameter :: aerospc( n_aerospc ) = (/
C                  -----------Name--------------     --------Min_Conc------
C                    Aitken    accum    coarse       Aitken  accum   coarse         density  isWet   volatile chrg visual emis
C                  --------- --------- ---------     -------- -------- ------------ ------- -------- -------- ---- ------ -----
     & spcs_type((/'ASO4I           ','ASO4J           ','ASO4K           '/),
     & (/1.0E-12,1.0E-6, conmin/),      1800.0, .false.,  .true.,   2,  3.0, 'PSO4            '),
     & spcs_type((/'ANH4I           ','ANH4J           ','ANH4K           '/),
     & (/conmin, conmin, conmin/),      1800.0, .false.,  .true.,  -1,  3.0, ' '),
     & spcs_type((/'ANO3I           ','ANO3J           ','ANO3K           '/),
     & (/conmin, conmin, conmin/),      1800.0, .false.,  .true.,   1,  3.0, 'PNO3            '),
     & spcs_type((/'AORGPAI         ','AORGPAJ         ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false.,  .true.,   0,  4.0, 'POC             '),
     & spcs_type((/'A25I            ','A25J            ','                '/),
     &(/conmin, conmin, conmin/),      2200.0, .false.,  .true.,   0,  1.0, 'PMFINE          '),
     & spcs_type((/'AECI            ','AECJ            ','                '/),
     & (/conmin, conmin, conmin/),      2200.0, .false.,  .true.,   0, 10.0, 'PEC             '),
     & spcs_type((/'ANAI            ','ANAJ            ','ANAK            '/),
     & (/conmin, conmin, conmin/),      2200.0, .false.,  .true.,  -1,  0.0, ' '),
     & spcs_type((/'ACLI            ','ACLJ            ','ACLK            '/),
     & (/conmin, conmin, conmin/),      2200.0, .false.,  .true.,   1,  0.0, ' '),
     & spcs_type((/'AH2OI           ','AH2OJ           ','AH2OK           '/),
     & (/conmin, conmin, conmin/),      1000.0,  .true.,  .true.,   0,  0.0, ' '),
     & spcs_type((/'                ','AALKJ           ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AXYL1J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AXYL2J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AXYL3J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','ATOL1J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','ATOL2J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','ATOL3J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','ABNZ1J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','ABNZ2J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','ABNZ3J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','ATRP1J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','ATRP2J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AISO1J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AISO2J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AISO3J          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','ASQTJ           ','                '/),
     & (/conmin, conmin, conmin/),      2000.0,  .true.,  .true.,   0,  4.0, ' '),
     & spcs_type((/'                ','AOLGAJ          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','AOLGBJ          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','AORGCJ          ','                '/),
     & (/conmin, conmin, conmin/),      2000.0, .false., .false.,   0,  4.0, ' '),
     & spcs_type((/'                ','                ','ASOIL           '/),
     &(/conmin, conmin, conmin/),      2600.0, .false.,  .true.,   0,  0.0, ' '),
     & spcs_type((/'                ','                ','ACORS           '/),
     & (/conmin, conmin, 1.889544E-05/),2200.0, .false.,  .true.,   0,  0.0, 'PMC             ')/)
     
C Required species
      Character( 16 ), Private, Parameter :: req_so4   = 'ASO4'  
      Character( 16 ), Private, Parameter :: req_nh4   = 'ANH4'  
      Character( 16 ), Private, Parameter :: req_no3   = 'ANO3'  
      Character( 16 ), Private, Parameter :: req_orgpa = 'AORGPA'  
      Character( 16 ), Private, Parameter :: req_na    = 'ANA'   
      Character( 16 ), Private, Parameter :: req_cl    = 'ACL'   
      Character( 16 ), Private, Parameter :: req_h2o   = 'AH2O'  
      Character( 16 ), Private, Parameter :: req_soil  = 'ASOIL'
      Character( 16 ), Private, Parameter :: req_cors  = 'ACORS'  

C Indexes of required species
      Integer :: aso4_idx
      Integer :: anh4_idx
      Integer :: ano3_idx
      Integer :: aorgpa_idx
      Integer :: ana_idx
      Integer :: acl_idx
      Integer :: ah2o_idx 
      Integer :: asoil_idx
      Integer :: acors_idx

      Real :: aerospc_mw( n_aerospc )          ! aero species M.W. (from AE_SPC.EXT) [ g/mol ]
      Real :: aerospc_conc( n_aerospc,n_mode ) ! aero species concentration [ ug/m^3 ]

C Common factors
      Real( 8 ) :: h2ofac                      ! converts mass concentrations [ug/m3] to 3rd moment concentrations [m3/m3]

C-------------------------------------------------------------------------------------------------------

      Type mode_type
         Character( 16 ) :: num_name     ! name of aerosol number variable
         Character( 16 ) :: srf_name     ! name of aerosol surface area variable
         Real            :: min_numconc  ! minimum number concentration
         Real            :: min_m2conc   ! minimum 2nd moment concentration
         Real            :: min_m3conc   ! minimum 3rd moment concentration
         Real            :: min_num      ! minimum number
         Real            :: min_m2       ! minimum 2nd moment
      End Type mode_type

      Type ( mode_type ), Parameter :: aeromode( n_mode ) = (/
C                   number     surface   minimum  minimum minimum  minimum   minimum
C                    name       name     numconc  m2conc  m3conc    num       m2
C                  ----------  ---------- ---------  ---------  ---------  ---------  ------------
     &   mode_type('NUMATKN', 'SRFATKN', conmin, conmin, conmin, 0.298852, 0.524834E-16),
     &   mode_type('NUMACC ', 'SRFACC ', conmin, conmin, conmin, 356.0190, 0.456022E-11),
     &   mode_type('NUMCOR ', 'SRFCOR ', conmin, conmin, conmin, 1.000000, 0.346714E-11)/)

      Real, Parameter :: min_sigma_g = 1.05
      Real, Parameter :: max_sigma_g = 2.50
      Real, Parameter :: def_sigma_g = 1.70

      Real            :: moment0_conc( n_mode )     ! 0th moment concentration
      Real            :: moment2_conc( n_mode )     ! 2nd moment concentration
      Real            :: moment3_conc( n_mode )     ! 3rd moment concentration

C Mass concentration (calculated by GETPAR)
      Real :: aeromode_mass( n_mode )   ! [ ug/m^3 ]

C Particle density (calculated by GETPAR)
      Real :: aeromode_dens( n_mode )   ! [ kg/m^3 ]

C Geometric mean diameter (calculated by GETPAR)
      Real :: aeromode_diam( n_mode )   ! [ m ]

C Log of geometric standard deviation (calculated by GETPAR )
      Real :: aeromode_sdev( n_mode )

C Mapping for loading from and unloading to CGRID array
      Integer :: aerospc_map( n_aerospc,n_mode )  ! indices of aero species to CGRID
      Integer :: aeronum_map( n_mode )      ! indices of aero number variable to CGRID
      Integer :: aerosrf_map( n_mode )      ! indices of aero surf area variable to CGRID

C Private variables for loading from and unloading to CGRID array
      Logical, Private, Save :: mapped = .False.
      Character( 16 ), Private, Save :: pname = 'Aero_Data'

      Contains

C-----------------------------------------------------------------------
      Subroutine map_aero()

C  Defines aerosol mapping from CGRID for species concentration and moments.
 
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Use cgrid_spcs
      Use met_data

      Implicit None

C Includes:
      Include SUBST_AE_SPC    ! aerosol species table
      Include SUBST_IOPARMS   ! I/O parameters definitions

C External Functions:
      Integer index1

C Local Variables:
      Character( 80 ) :: xmsg
      Integer m, n, spc

      if ( .Not. mapped ) Then
         Call cgrid_map ( nspcsd, gc_strt, ae_strt, nr_strt, tr_strt )
         aerospc_mw = 0.0
         aerospc_map = 0

C build mapping to CGRID for each aero spc
         Do m = 1, n_mode
            Do spc = 1, n_aerospc
               If ( aerospc( spc )%name( m ) .Ne. ' ' ) Then
                  n = index1( aerospc( spc )%name( m ), n_ae_spc, ae_spc )
                  If ( n .Eq. 0 ) Then
                     xmsg = 'Species '// aerospc(spc)%name( m ) // ' in aerospc name is not in AE_SPC'
                     Call m3exit( pname, 0, 0, xmsg, xstat3 )
                  End If
                  aerospc_map( spc,m ) = ae_strt - 1 + n
                  If ( aerospc_mw( spc ) .Lt. 0.5 ) Then  ! mw=0 means a new species
                     aerospc_mw( spc ) = ae_molwt( n )
                  Else If ( aerospc_mw( spc ) .Ne. ae_molwt( n ) ) Then
                     xmsg = 'molecular weight of ' // Trim( aerospc( spc )%name( m ) )
     &                    // ' is different from that of the same species'
     &                    // ' in the same or another mode.'
                   Call m3exit( pname, 0, 0, xmsg, xstat3 )
                   End If
               End If
            End Do
         End Do

C Build mapping to CGRID for aero # and surf area variables
         aeronum_map = 0
         aerosrf_map = 0

         Do m = 1, n_mode
            n = index1( aeromode( m )%num_name, n_ae_spc, ae_spc )
            If ( n .Eq. 0 ) Then
               xmsg = 'Species ' // Trim( aeromode( m )%num_name )
     &              //' in aeronum name is not in AE_SPC'
               Call m3exit( pname, 0, 0, xmsg, xstat3 )
            Else
               aeronum_map( m ) = ae_strt - 1 + n
            End If

            n = index1( aeromode( m )%srf_name, n_ae_spc, ae_spc )
            If ( n .Eq. 0 ) Then
               xmsg = 'species ' // Trim( aeromode( m )%srf_name )
     &              // ' in aerosrf name is not in AE_SPC'
               Call m3exit( pname, 0, 0, xmsg, xstat3 )
            Else
               aerosrf_map( m ) = ae_strt - 1 + n
            End If
         End Do

C Find indices of required species
         aso4_idx   = findAero( req_so4,   .True. )
         anh4_idx   = findAero( req_nh4,   .True. )
         ano3_idx   = findAero( req_no3,   .True. )
         aorgpa_idx = findAero( req_orgpa, .True. )
         ana_idx    = findAero( req_na,    .True. )
         acl_idx    = findAero( req_cl,    .True. )
         ah2o_idx   = findAero( req_h2o,   .True. )
         asoil_idx  = findAero( req_soil,  .True. )
         acors_idx  = findAero( req_cors,  .True. )
 
C Compute common factors
         h2ofac = 1.0E-9 * f6dpi / aerospc( ah2o_idx )%density

         mapped = .True.

      End If   ! Not mapped

      Return
      End Subroutine map_aero

C-----------------------------------------------------------------------
      Subroutine extract_aero( conc )

C  Extracts the required aerosol data from CGRID into conc.
 
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Implicit None

C Arguments:
      Real, Intent( In ) :: conc( : )

C Local Variables:
      Character( 80 ) :: xmsg
      Integer m, n, spc 

      If ( .Not. mapped ) Then
         Call map_aero()
      End If

C copy grid cell concentrations of aero species to aerospc_conc
      aerospc_conc = 0.0
      Do m = 1, n_mode
         Do spc = 1, n_aerospc
            n = aerospc_map( spc,m )
            If ( n .ne. 0 ) Then
               aerospc_conc( spc,m ) = Max( conc( n ), aerospc( spc )%min_conc( m ) ) ! ug/m^3
            End If
         End Do
      End Do

C copy grid cell concentrations of aero # and surf area
C convert and assign to moment0_conc and moment2_conc

      moment0_conc = 0.0
      moment2_conc = 0.0

      Do m = 1, n_mode
         n = aeronum_map( m )
         moment0_conc( m ) = Max( conc( n ), aeromode( m )%min_numconc )
         n = aerosrf_map( m )
         moment2_conc( m ) = Max( conc( n ), aeromode( m )%min_m2conc )
      End Do

      Return
      End Subroutine extract_aero

C-----------------------------------------------------------------------
      Subroutine update_aero( conc )

C  Updates aerosol values in CGRID from the conc array.
 
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Use met_data      !fundamental constants, data type definitions, etc.

      Implicit None

C arguments:
      Real, Intent( Out ) :: conc( : )

C includes:
      Include SUBST_IOPARMS   ! I/O parameters definitions

C local variables:
      Character( 80 ) :: xmsg
      Integer m, n, spc 

      If ( .Not. mapped ) Then
         xmsg = 'CGRID Species has not been mapped'
         Call m3exit( pname, 0, 0, xmsg, xstat3 )
      End If

C Copy aerospc_conc back to grid cell concentrations

      Do m = 1, n_mode
         Do spc = 1, n_aerospc
            n = aerospc_map( spc,m )
            If ( n .Ne. 0 ) Then
               conc( n ) = Max( aerospc_conc( spc,m ), aerospc( spc )%min_conc( m ) )
            ! Write( *,'(a,i5,g16.6)' ) aerospc( spc )%name( m ), n, CONC( n )
            End If
         End Do
      End Do

C Copy aero number and surface area back to grid cell concentrations

      Do m = 1, n_mode
         n = aeronum_map( m )
         conc( n ) = Max( moment0_conc( m ), aeromode( m )%min_numconc )
         ! Write( *,'(a,i5,g16.6)' ) aeromode( m )%num_name, n, CONC( n )
         n = aerosrf_map( m )
         conc( n ) = pi * moment2_conc( m ) 
         ! Write( *,'(a,i5,g16.6)' ) aeromode( m )%srf_name, n, CONC( n )
      End Do

      Return
      End Subroutine update_aero   

C-----------------------------------------------------------------------
      Function findAero( vname, required ) Result ( idx )

C  Finds the index of 'required' aerosol species in the aerospc list
 
C  Revision History:
C     First version was coded in April 2010 by Steve Howard with
C     Prakash Bhave, Jeff Young, and Sergey Napelenok.
C-----------------------------------------------------------------------

      Implicit None

C Arguments:
      Character( * ) :: vname
      Logical required
      Integer idx

C Includes:
      Include SUBST_IOPARMS   ! I/O parameters definitions

C Local Variables:
      Character( 80 ) :: xmsg
      Integer spc, n

      Do n = 1, n_mode
         Do spc = 1, n_aerospc
            If ( Index( aerospc( spc )%name( n ), Trim( vname ) ) .Gt. 0 ) Then
               idx = spc
               ! Write( *,'(a,i5)' ) name, idx
               Return
            End If
         End Do
      End Do

      If ( .Not. required ) Then
         idx = 0
         Return
      End If

      xmsg = 'Required Species ' // Trim( vname ) // ' Not found in aerospc names array'
      Call m3exit( pname, 0, 0, xmsg, xstat3 )

      Return
      End Function findAero
      
      End Module aero_data

