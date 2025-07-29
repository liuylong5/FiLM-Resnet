# mechanism and tracer: cb05cl_ae5_aq, trac0

 MODEL =  bwd_intel
 current_folder := $(shell pwd)
 FC   = mpiifort

 F_FLAGS = -O2 -extend-source -save -I.  \
 -I/dssg/home/acct-esehazenet/share/public_model/CMAQ-5.0-Adjoint/lib/stenex/se

##-Ktrap=divz,inv,ovf
##unf is very common in the code

#-Ktrap=fp .eqv. -Ktrap=divz,inv,ovf
#divz -- divide by zero
#inv -- invalid operands
#ovf -- overflow
#unf -- underflow
#inexact -- inexact result
#denorm -- denormalized operands

 CC   = mpicc
 C_FLAGS    = -O2
 CPP_FLAGS = \
 -DINTERPB=PINTERPB \
 -DIncl \
 -Dparallel \
 -Dcluster \
 -DM3WARN=PM3WARN \
 -DSHUT3=PSHUT3 \
 -DWRITE3=PWRITE3 \
 -DSUBST_MODULES=se_MODULES \
 -DSUBST_BARRIER=se_BARRIER \
 -DSUBST_GLOBAL_MAX=se_GLOBAL_MAX \
 -DSUBST_GLOBAL_MIN=se_GLOBAL_MIN \
 -DSUBST_GLOBAL_MIN_DATA=se_GLOBAL_MIN_DATA \
 -DSUBST_GLOBAL_TO_LOCAL_COORD=se_GLOBAL_TO_LOCAL_COORD \
 -DSUBST_GLOBAL_SUM=se_GLOBAL_SUM \
 -DSUBST_GLOBAL_LOGICAL=se_GLOBAL_LOGICAL \
 -DSUBST_LOOP_INDEX=se_LOOP_INDEX \
 -DSUBST_SUBGRID_INDEX=se_SUBGRID_INDEX \
 -DSUBST_HI_LO_BND_PE=se_HI_LO_BND_PE \
 -DSUBST_SUM_CHK=se_SUM_CHK \
 -DSUBST_INIT_ARRAY=se_INIT_ARRAY \
 -DSUBST_COMM=se_COMM \
 -DSUBST_MY_REGION=se_MY_REGION \
 -DSUBST_SLICE=se_SLICE \
 -DSUBST_GATHER=se_GATHER \
 -DSUBST_DATA_COPY=se_DATA_COPY \
 -DSUBST_IN_SYN=se_IN_SYN \

 LINK_FLAGS = 

 SEF90 = -L/dssg/home/acct-esehazenet/share/public_model/CMAQ-5.0-Adjoint/lib/stenex/se -lse_snl
 PARIO = -L/dssg/home/acct-esehazenet/share/public_model/CMAQ-5.0-Adjoint/lib/pario -lpario
 IOAPI = -L/dssg/home/acct-esehazenet/share/public_model/Libs/ioapi-3.2-20200828/Linux2_x86_64ifort -lioapi
 NETCDF = -L/dssg/home/acct-esehazenet/share/public_model/Libs/netcdf/lib -lnetcdf
 NETCDFF = -L/dssg/home/acct-esehazenet/share/public_model/Libs/netcdf/lib -lnetcdff
 MPI = -L/dssg/opt/icelake/linux-centos8-icelake/gcc-8.5.0/intel-oneapi-mpi-2021.4.0-5xzcg7kxlx4ytkz2plqalrg424so43cy/mpi/2021.4.0/lib -lmpifort
 LIBRARIES = $(SEF90) $(PARIO) $(IOAPI) $(NETCDF) $(NETCDFF) $(MPI)

 BASE_INC =/dssg/home/acct-esehazenet/share/public_code/software_source_files/CMAQ-4.7.1/models/include/release
 IOAPI_INC = /dssg/home/acct-esehazenet/share/public_model/Libs/ioapi-3.2-20200828/ioapi/fixed_src
 MECH_INC = $(BASE_INC)/cb05cl_ae5_aq
 TRAC_INC = $(BASE_INC)/trac0

 INCLUDES = \
 -DSUBST_MPICH=\"/dssg/opt/icelake/linux-centos8-icelake/gcc-8.5.0/intel-oneapi-mpi-2021.4.0-5xzcg7kxlx4ytkz2plqalrg424so43cy/mpi/2021.4.0/include/mpif.h\" \
 -DSUBST_GRID_ID= \
 -DSUBST_PE_COMM=\"$(BASE_INC)/PE_COMM.EXT\" \
 -DSUBST_CONST=\"$(BASE_INC)/CONST.EXT\" \
 -DSUBST_FILES_ID=\"$(BASE_INC)/FILES_CTM.EXT\" \
 -DSUBST_EMPR_VD=\"$(BASE_INC)/EMISPRM.vdif.EXT\" \
 -DSUBST_EMPR_CH=\"$(BASE_INC)/EMISPRM.chem.EXT\" \
 -DSUBST_IOPARMS=\"$(IOAPI_INC)/PARMS3.EXT\" \
 -DSUBST_IOFDESC=\"$(IOAPI_INC)/FDESC3.EXT\" \
 -DSUBST_IODECL=\"$(IOAPI_INC)/IODECL3.EXT\" \
 -DSUBST_RXCMMN=\"$(MECH_INC)/RXCM.EXT\" \
 -DSUBST_RXDATA=\"$(MECH_INC)/RXDT.EXT\" \
 -DSUBST_GC_SPC=\"$(MECH_INC)/GC_SPC.EXT\" \
 -DSUBST_GC_EMIS=\"$(MECH_INC)/GC_EMIS.EXT\" \
 -DSUBST_GC_ICBC=\"$(MECH_INC)/GC_ICBC.EXT\" \
 -DSUBST_GC_DIFF=\"$(MECH_INC)/GC_DIFF.EXT\" \
 -DSUBST_GC_DDEP=\"$(MECH_INC)/GC_DDEP.EXT\" \
 -DSUBST_GC_DEPV=\"$(MECH_INC)/GC_DEPV.EXT\" \
 -DSUBST_GC_ADV=\"$(MECH_INC)/GC_ADV.EXT\" \
 -DSUBST_GC_CONC=\"$(MECH_INC)/GC_CONC.EXT\" \
 -DSUBST_GC_G2AE=\"$(MECH_INC)/GC_G2AE.EXT\" \
 -DSUBST_GC_G2AQ=\"$(MECH_INC)/GC_G2AQ.EXT\" \
 -DSUBST_GC_SCAV=\"$(MECH_INC)/GC_SCAV.EXT\" \
 -DSUBST_GC_WDEP=\"$(MECH_INC)/GC_WDEP.EXT\" \
 -DSUBST_AE_SPC=\"$(MECH_INC)/AE_SPC.EXT\" \
 -DSUBST_AE_EMIS=\"$(MECH_INC)/AE_EMIS.EXT\" \
 -DSUBST_AE_ICBC=\"$(MECH_INC)/AE_ICBC.EXT\" \
 -DSUBST_AE_DIFF=\"$(MECH_INC)/AE_DIFF.EXT\" \
 -DSUBST_AE_DDEP=\"$(MECH_INC)/AE_DDEP.EXT\" \
 -DSUBST_AE_DEPV=\"$(MECH_INC)/AE_DEPV.EXT\" \
 -DSUBST_AE_ADV=\"$(MECH_INC)/AE_ADV.EXT\" \
 -DSUBST_AE_CONC=\"$(MECH_INC)/AE_CONC.EXT\" \
 -DSUBST_AE_A2AQ=\"$(MECH_INC)/AE_A2AQ.EXT\" \
 -DSUBST_AE_SCAV=\"$(MECH_INC)/AE_SCAV.EXT\" \
 -DSUBST_AE_WDEP=\"$(MECH_INC)/AE_WDEP.EXT\" \
 -DSUBST_NR_SPC=\"$(MECH_INC)/NR_SPC.EXT\" \
 -DSUBST_NR_EMIS=\"$(MECH_INC)/NR_EMIS.EXT\" \
 -DSUBST_NR_ICBC=\"$(MECH_INC)/NR_ICBC.EXT\" \
 -DSUBST_NR_DIFF=\"$(MECH_INC)/NR_DIFF.EXT\" \
 -DSUBST_NR_DDEP=\"$(MECH_INC)/NR_DDEP.EXT\" \
 -DSUBST_NR_DEPV=\"$(MECH_INC)/NR_DEPV.EXT\" \
 -DSUBST_NR_ADV=\"$(MECH_INC)/NR_ADV.EXT\" \
 -DSUBST_NR_N2AE=\"$(MECH_INC)/NR_N2AE.EXT\" \
 -DSUBST_NR_CONC=\"$(MECH_INC)/NR_CONC.EXT\" \
 -DSUBST_NR_N2AQ=\"$(MECH_INC)/NR_N2AQ.EXT\" \
 -DSUBST_NR_SCAV=\"$(MECH_INC)/NR_SCAV.EXT\" \
 -DSUBST_NR_WDEP=\"$(MECH_INC)/NR_WDEP.EXT\" \
 -DSUBST_TR_SPC=\"$(TRAC_INC)/TR_SPC.EXT\" \
 -DSUBST_TR_EMIS=\"$(TRAC_INC)/TR_EMIS.EXT\" \
 -DSUBST_TR_ICBC=\"$(TRAC_INC)/TR_ICBC.EXT\" \
 -DSUBST_TR_DIFF=\"$(TRAC_INC)/TR_DIFF.EXT\" \
 -DSUBST_TR_DDEP=\"$(TRAC_INC)/TR_DDEP.EXT\" \
 -DSUBST_TR_DEPV=\"$(TRAC_INC)/TR_DEPV.EXT\" \
 -DSUBST_TR_ADV=\"$(TRAC_INC)/TR_ADV.EXT\" \
 -DSUBST_TR_T2AQ=\"$(TRAC_INC)/TR_T2AQ.EXT\" \
 -DSUBST_TR_SCAV=\"$(TRAC_INC)/TR_SCAV.EXT\" \
 -DSUBST_TR_WDEP=\"$(TRAC_INC)/TR_WDEP.EXT\" \
 -DSUBST_PACTL_ID=\"$(BASE_INC)/pa_noop/PA_CTL.EXT\" \
 -DSUBST_PACMN_ID=\"$(BASE_INC)/pa_noop/PA_CMN.EXT\" \
 -DSUBST_PADAT_ID=\"$(BASE_INC)/pa_noop/PA_DAT.EXT\" \

 VPATH = src

 OBJECTS = \
 adStack.o \
 distr_env.o \
 VGRD_DEFN.o \
 CGRID_SPCS.o \
 HGRD_DEFN.o \
 MET_DATA.o \
 AERO_DATA.o \
 PMEM_DEFN.o \
 GRID_CONF.o \
 SOA_DEFN.o \
 PRECURSOR_DATA.o \
 CSQY_DATA.o \
 PAGRD_DEFN.o \
 precursor_data_b.o \
 kppdata_mod.o \
 DEPVVARS.o \
 DEPV_DEFN.o \
 AERO_DESC.o \
 CHK_FILES_DEFS.o \
 CHK_FILES_IMPL.o \
 AERO_EMIS.o \
 EMIS_DEFN.o \
 PCGRID_DEFN.o \
 AQ_DATA.o \
 aqchem_Precision.o \
 VDIFF_MAP.o \
 PHOT_MOD.o \
 aero_data_b.o \
 soa_defn_b.o \
 ADJOINT_FILES.o \
 CHK_FILES.o \
 KPP_Precision.o \
 KPP_Parameters.o \
 KPP_Util.o \
 KPP_Global.o \
 KPP_JacobianSP.o \
 KPP_LinearAlgebra.o \
 KPP_Jacobian.o \
 KPP_Function.o \
 KPP_HessianSP.o \
 KPP_Hessian.o \
 KPP_Integrator_adj.o \
 KPP_Model.o \
 aqchem_Parameters.o \
 aqchem_Global.o \
 aqchem_Rates.o \
 aqchem_Function.o \
 aqchem_JacobianSP.o \
 aqchem_LinearAlgebra.o \
 aqchem_Jacobian.o \
 aqchem_Integrator.o \
 aqchem_HessianSP.o \
 aqchem_Hessian.o \
 aqchem_Integrator_ADJ.o \
 aqchem_initialize_b.o \
 aqchem_Model.o \
 zfdbc.o \
 ADV_CADJ.o \
 zadvyppm_bwd.o \
 yadv_rj_fwd.o \
 yadv_bwd.o \
 xadv_rj_fwd.o \
 xadv_bwd.o \
 write3_distr.o \
 vppm.o \
 vppm_bwd.o \
 volinorg_adj.o \
 aero_emis_b.o \
 emis_defn_b.o \
 vdiffacm2_bwd.o \
 tri_bwd.o \
 syncstep.o \
 subhfile.o \
 subhdomain.o \
 setup_logdev.o \
 set_aeconcmin.o \
 sciproc_bwd.o \
 scavwdep.o \
 scavwdep_b.o \
 rho_j.o \
 rescld.o \
 rescld_b.o \
 readconf.o \
 rdmet.o \
 rddepv.o \
 rdbcon_rj.o \
 rdbcon_adj.o \
 ppmv_msmr.o \
 phot.o \
 pa_update.o \
 par_term.o \
 par_init.o \
 pa_output.o \
 pa_mkhdr.o \
 pa_irr.o \
 pa_irr_ctl.o \
 pa_init.o \
 opvis.o \
 opphot.o \
 opmet.o \
 opemis.o \
 opdiam.o \
 opdepv_diag.o \
 opddep.o \
 newpart3_adj.o \
 n2o5prob_adj.o \
 matrix_bwd.o \
 m3dry.o \
 lstepf.o \
 kppinit.o \
 kppdriver_adj.o \
 kppcalcks.o \
 isoropiaIIcode_adj.o \
 intracoag_gh_adj.o \
 intercoag_gh_adj.o \
 init_lgrid.o \
 indexn.o \
 hveloc.o \
 hppm.o \
 hppm_bwd.o \
 hlconst.o \
 hetchem.o \
 hetchem_b.o \
 hetchem_adj.o \
 hdiff_bwd.o \
 hcontvel.o \
 hcond3_adj.o \
 hcdiff3d.o \
 hadv_bwd.o \
 grdcheck.o \
 getsc_adj.o \
 getpar.o \
 getpar_b.o \
 get_envlist.o \
 getcoags_adj.o \
 getalpha.o \
 getalpha_b.o \
 gas_depv_map.o \
 flcheck.o \
 findex.o \
 eddyx.o \
 driver_bwd.o \
 diffmsg.o \
 deform.o \
 define_receptor.o \
 couple.o \
 couple_cadj.o \
 convcld_acm.o \
 convcld_acm_b.o \
 conv_cgrid_bwd.o \
 compute_flux_adj.o \
 coags.o \
 coags_b.o \
 cldproc_b.o \
 cksummer.o \
 chkgrid.o \
 cgrid_depv.o \
 calc_h2o_adj.o \
 aq_map.o \
 aq_map_b.o \
 aqchem_newton.o \
 aqchem_newton_b.o \
 aqchem_b.o \
 aero_subs.o \
 aero_subs_b.o \
 aero_driver_b.o \
 aero_depv.o \
 advstep.o \
 advbc_map_adj.o \
 adBuffer.o \
 acmcld.o \
 acmcld_b.o \
 acmcld_adj.o \
 aqchem_kpp.o


.SUFFIXES: .F .f .c .f90 .F90



CXX := /dssg/opt/icelake/linux-centos8-icelake/gcc-8.5.0/gcc-11.2.0-va6hwyzhi4tyigabpu5aqp4f5xz55c6p/bin/g++
LIBTORCH_DIR := /dssg/home/acct-esehazenet/hazenet-liuyulong/temp/libtorch
TORCH_CXX_FLAGS := -D_GLIBCXX_USE_CXX11_ABI=0
TORCH_LIBRARIES := -ltorch -lc10 -lkineto -lstdc++fs \
  /dssg/home/acct-esehazenet/hazenet-liuyulong/temp/libtorch/lib/libc10.so \
  /dssg/home/acct-esehazenet/hazenet-liuyulong/temp/libtorch/lib/libkineto.a

libDnnInfer.so: update_lambda.cpp
	$(CXX) -shared $< -o $@ \
		-O3 -fPIC -std=c++17 -g \
		-I${LIBTORCH_DIR}/include/ -L${LIBTORCH_DIR}/lib/ \
		${TORCH_CXX_FLAGS} \
		${TORCH_LIBRARIES} \
		-Wl,-rpath,/dssg/home/acct-esehazenet/hazenet-liuyulong/temp/libtorch/lib

$(MODEL): $(OBJECTS) libDnnInfer.so
	$(FC) $(LINK_FLAGS) $(OBJECTS) $(LIBRARIES) -qopenmp -L./ -lDnnInfer -Wl,-rpath,$(current_folder)/libDnnInfer -o $@


.F90.o:
	$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.F.o:
	$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f90.o:
	$(FC) $(F_FLAGS) -c $<

.f.o:
	$(FC) $(F_FLAGS) -c $<

.c.o:
	$(CC) $(C_FLAGS) -c $<

clean:
	rm -f *.o *.mod
