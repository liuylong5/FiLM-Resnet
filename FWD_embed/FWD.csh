module purge
module load oneapi
export LD_LIBRARY_PATH=../public_model/Libs/netcdf/lib:${LD_LIBRARY_PATH}
export MODEL_PATH=../FWD/2020_12km_embed/path
export LD_LIBRARY_PATH=../FWD/2020_12km_embed:$LD_LIBRARY_PATH 
EXEC=../FWD/2020_12km_embed/fwd_intel 
export work_dir=../FWD/2020_12km_embed 
export MY_LOG_FILE=../FWD/2020_12km_embed/FWD.log
npcol=4
nprow=4
output_dir=../FWD/2020_12km_embed
icon_dir=
bcon_dir=
export GRID_NAME=China_12km
export GRIDDESC=../FWD/2020_12km_embed/GRIDDESC
OCEANpath=../FWD/2020_12km_embed
OCEANfile=ocean_file.nc
export MY_LOG_FILE=../FWD/2020_12km_embed/FWD.log
METpath=../FWD/2020_12km_embed
EMISpath=../FWD/2020_12km_embed
JVALpath=../FWD/2020_12km_embed
let NPROCS=npcol*nprow 
export NPCOL_NPROW="$npcol $nprow" 
START_DATE=$(date -d "2020-01-01" +'%Y-%m-%d')
END_DATE=$(date -d "2020-01-07" +'%Y-%m-%d')
export AVG_CONC_SPCS="ALL"
export ACONC_BLEV_ELEV=" 1 36"
export CTM_EMLAYS=0
export CREATE_CHK=Y

export CTM_MAXSYNC=240
export CTM_MINSYNC=240
export CTM_AERDIAG=Y
export CTM_PHOTDIAG=Y
export KZMIN=Y
export IOAPI_LOG_WRITE=F
export FL_ERR_STOP=F
export CTM_APPL="fwd_intel"
export FLOOR_FILE=$work_dir/floor
export CONC_SPCS="CO NO2 NO NO3 HNO3 HO2 N2O5 NH3 O3 OH SO2 AECI AECJ ASO4I ASO4J ASO4K ANH4I ANH4J ANH4K ANO3I ANO3J ANO3K"
export CONC_BLEV_ELEV=" 1 1"
export CTM_WVEL=T
export IOAPI_CHECK_HEADERS=F

cd $work_dir
START_J=$(date -ud "${START_DATE}" +'%Y%j')
END_J=$(date -ud "${END_DATE}" +'%Y%j')
STTIME=000000        # beginning GMT time (HHMMSS)
NSTEPS=240000        # time duration (HHMMSS) for this run
TSTEP=010000        # output time step interval (HHMMSS)

TODAYG=${START_DATE}
TODAYJ=$(date -ud "${TODAYG}" +'%Y%j')
YESTERDAYG=$(date -ud "${TODAYG}-1days" +'%Y-%m-%d')
while (( "${TODAYJ}" <= "${END_J}" )); do
echo "Forward model of $TODAYG was started!"
EMISfile=emission_${TODAYG}.nc

# if [ "$TODAYJ" = "$START_J" ]; then
# if [ ${#icon_dir} -gt 5 ]; then
# GC_ICpath=${icon_dir}
# GC_ICfile=CGRID_${YESTERDAYG}.nc
# else
# GC_ICpath=${work_dir}
# GC_ICfile=ICON_V5f_${TODAYG}.nc
# fi
# else
# GC_ICpath=${output_dir}
# GC_ICfile=CGRID_${YESTERDAYG}.nc
# fi

if [ "$TODAYJ" = "$START_J" ]; then
    GC_ICpath==${icon_dir}
    GC_ICfile=CGRID_2019-12-31.nc  
else
    GC_ICpath=${output_dir}
    GC_ICfile=CGRID_${YESTERDAYG}.nc
fi


if [ ${#bcon_dir} -gt 5 ]; then
GC_BCpath=${bcon_dir}
GC_BCfile=ACONC_${TODAYG}.nc
else
GC_BCpath=${work_dir}
GC_BCfile=BCON_V5f_${TODAYG}.nc
fi

GC2file=GRIDCRO2D_${TODAYG}.nc
GD2file=GRIDDOT2D_${TODAYG}.nc
MC2file=METCRO2D_${TODAYG}.nc
MD3file=METDOT3D_${TODAYG}.nc
MC3file=METCRO3D_${TODAYG}.nc
MB3file=METBDY3D_${TODAYG}.nc
JVALfile=JTABLE_${TODAYG}.dat

TR_DVpath=$METpath
TR_DVfile=$MC2file
AE_ICpath=$GC_ICpath
NR_ICpath=$GC_ICpath
TR_ICpath=$GC_ICpath
AE_ICfile=$GC_ICfile
NR_ICfile=$GC_ICfile
TR_ICfile=$GC_ICfile
AE_BCpath=$GC_BCpath
NR_BCpath=$GC_BCpath
TR_BCpath=$GC_BCpath
AE_BCfile=$GC_BCfile
NR_BCfile=$GC_BCfile
TR_BCfile=$GC_BCfile


if [ $?EMISpath2 ]; then
export EMIS_SUP=$EMISpath2/$EMISfile2
if [ ! -e $EMIS_SUP ]; then
echo " $EMIS_SUP not used "
fi
fi

if [ $?TR_DVpath ]; then
export DEPV_TRAC_1=$TR_DVpath/$TR_DVfile
if [ ! -e $DEPV_TRAC_1 ]; then
echo " $DEPV_TRAC_1 not found "
exit 1
fi
fi

if [ $?TR_EMpath ]; then
export EMIS_TRAC_1=$TR_EMpath/$TR_EMfile
if [ ! -e $EMIS_TRAC_1 ]; then
echo " $EMIS_TRAC_1 not found "
exit 1
fi
fi

# for OCEAN file ...

if [ $?OCEANpath ]; then
export OCEAN_1=$OCEANpath/$OCEANfile
if [ ! -e $OCEAN_1 ]; then
echo  " $OCEAN_1 not found "
exit 1
fi
fi

export EMIS_1=$EMISpath/$EMISfile
export INIT_GASC_1=$GC_ICpath/$GC_ICfile
export BNDY_GASC_1=$GC_BCpath/$GC_BCfile
export INIT_AERO_1=$AE_ICpath/$AE_ICfile
export BNDY_AERO_1=$AE_BCpath/$AE_BCfile
export INIT_NONR_1=$NR_ICpath/$NR_ICfile
export BNDY_NONR_1=$NR_BCpath/$NR_BCfile
export INIT_TRAC_1=$TR_ICpath/$TR_ICfile
export BNDY_TRAC_1=$TR_BCpath/$TR_BCfile

export GRID_DOT_2D=$METpath/$GD2file
export GRID_CRO_2D=$METpath/$GC2file
export MET_CRO_2D=$METpath/$MC2file
export MET_DOT_3D=$METpath/$MD3file
export MET_CRO_3D=$METpath/$MC3file
export MET_BDY_3D=$METpath/$MB3file
export LAYER_FILE=$METpath/$MC3file
export XJ_DATA=$JVALpath/$JVALfile

flist="$EMIS_1
$INIT_GASC_1
$BNDY_GASC_1
$INIT_AERO_1
$BNDY_AERO_1
$INIT_NONR_1
$BNDY_NONR_1
$INIT_TRAC_1
$BNDY_TRAC_1
$GRID_DOT_2D
$GRID_CRO_2D
$MET_CRO_2D
$MET_DOT_3D
$MET_CRO_3D
$MET_BDY_3D
$LAYER_FILE
$XJ_DATA"

for file in $flist; do
if [ ! -e $file ]; then
echo " $file not found "
exit 1
fi
done

# action if the output files already exist ...
#> output files and directories
CONCfile="CONC_"${TODAYG}.nc               # CTM_CONC_1
ACONCfile="ACONC_"${TODAYG}.nc              # CTM_ACONC_1
CGRIDfile="CGRID_"${TODAYG}.nc              # CTM_CGRID_1
SGRIDfile="SGRID_"${TODAYG}.nc              # sens file
DD1file="DRYDEP_"${TODAYG}.nc             # CTM_DRY_DEP_1
WD1file="WETDEP1_"${TODAYG}.nc            # CTM_WET_DEP_1
WD2file="WETDEP2_"${TODAYG}.nc            # CTM_WET_DEP_2
SS1file="SSEMIS1_"${TODAYG}.nc            # CTM_SSEMIS_1
AV1file="AEROVIS_"${TODAYG}.nc            # CTM_VIS_1
AD1file="AERODIAM_"${TODAYG}.nc           # CTM_DIAM_1
PA1file="PA_1_"${TODAYG}.nc               # CTM_IPR_1
PA2file="PA_2_"${TODAYG}.nc               # CTM_IPR_2
PA3file="PA_3_"${TODAYG}.nc               # CTM_IPR_3
IRR1file="IRR_1_"${TODAYG}.nc              # CTM_IRR_1
IRR2file="IRR_2_"${TODAYG}.nc              # CTM_IRR_2
IRR3file="IRR_3_"${TODAYG}.nc              # CTM_IRR_3
RJ1file="RJ_1_"${TODAYG}.nc               # CTM_RJ_1
RJ2file="RJ_2_"${TODAYG}.nc               # CTM_RJ_2

export CTM_CONC_1="$output_dir/$CONCfile -v"
export A_CONC_1="$output_dir/$ACONCfile -v"
export S_CGRID="$output_dir/$CGRIDfile -v"
export S_SGRID="$output_dir/$SGRIDfile -v"
export CTM_DRY_DEP_1="$output_dir/$DD1file -v"
export CTM_WET_DEP_1="$output_dir/$WD1file -v"
export CTM_WET_DEP_2="$output_dir/$WD2file -v"
export CTM_SSEMIS_1="$output_dir/$SS1file -v"
export CTM_VIS_1="$output_dir/$AV1file -v"
export CTM_DIAM_1="$output_dir/$AD1file -v"
export CTM_IPR_1="$output_dir/$PA1file -v"
export CTM_IPR_2="$output_dir/$PA2file -v"
export CTM_IPR_3="$output_dir/$PA3file -v"
export CTM_IRR_1="$output_dir/$IRR1file -v"
export CTM_IRR_2="$output_dir/$IRR2file -v"
export CTM_IRR_3="$output_dir/$IRR3file -v"
export CTM_RJ_1="$output_dir/$RJ1file -v"
export CTM_RJ_2="$output_dir/$RJ2file -v"
export DBG_IGRID="$output_dir/DBG_iCONVCLD_${TODAYG}.nc -v"
export DBG_RGRID="$output_dir/DBG_rCONVCLD_${TODAYG}.nc -v"
export ADJ_CHEM_CHK="$output_dir/CHEM_CHK_${TODAYG}.nc -v"
export ADJ_CHEM_END="$output_dir/CHEM_END_${TODAYG}.nc -v"
export ADJ_CHEM_START="$output_dir/CHEM_STA_${TODAYG}.nc -v"
export ADJ_VDIFF_CHK="$output_dir/VDIFF_CHK_${TODAYG}.nc -v"
export ADJ_HA_RHOJ_CHK="$output_dir/HA_RHOJ_CHK_${TODAYG}.nc -v"
export ADJ_VA_RHOJ_CHK="$output_dir/VA_RHOJ_CHK_${TODAYG}.nc -v"
export ADJ_AERO_CHK="$output_dir/AERO_CHK_${TODAYG}.nc -v"
export ADJ_HADV_CHK="$output_dir/HADV_CHK_${TODAYG}.nc -v"
export ADJ_VADV_CHK="$output_dir/VADV_CHK_${TODAYG}.nc -v"
export ADJ_CLD_CHK="$output_dir/CLD_CHK_${TODAYG}.nc -v"
export ADJ_EMIS_CHK="$output_dir/EMIS_CHK_${TODAYG}.nc -v"
export ADJ_EMISU_CHK="$output_dir/EMISU_CHK_${TODAYG}.nc -v"
export CTM_XFIRST_OUT="$output_dir/XFIRST_${TODAYG}"

flist="$CTM_CONC_1
$A_CONC_1
$CTM_DRY_DEP_1
$CTM_WET_DEP_1
$CTM_WET_DEP_2
$CTM_SSEMIS_1
$CTM_VIS_1
$CTM_DIAM_1
$CTM_IPR_1
$CTM_IPR_2
$CTM_IPR_3
$CTM_IRR_1
$CTM_IRR_2
$CTM_IRR_3
$CTM_RJ_1
$CTM_RJ_2
$DBG_IGRID
$DBG_RGRID
$ADJ_CHEM_CHK
$ADJ_CHEM_END
$ADJ_CHEM_START
$ADJ_VDIFF_CHK
$ADJ_HA_RHOJ_CHK
$ADJ_VA_RHOJ_CHK
$ADJ_AERO_CHK
$ADJ_HADV_CHK
$ADJ_VADV_CHK
$ADJ_CLD_CHK
$ADJ_EMIS_CHK
$ADJ_EMISU_CHK
$CTM_XFIRST_OUT"

for file in $flist; do
if [ $file != '-v' ]; then
if [ -e $file ]; then
/bin/rm $file
fi
fi
done

export CTM_STDATE=$TODAYJ
export CTM_STTIME=$STTIME
export CTM_RUNLEN=$NSTEPS
export CTM_TSTEP=$TSTEP
export CTM_PROGNAME='fwd_intel'

rm CTM* >&/dev/null
{ time -p mpirun -np $NPROCS $EXEC; } > ${MY_LOG_FILE}
if [ $? -ne 0 ]; then
echo "Model run failed for $TODAYG"
exit 1  # Exit the script if the model run fails
fi
mv FLOOR_* ${FLOOR_FILE} >&/dev/null
echo "Forward model of ${TODAYG} was finished"

#> Increment both Gregorian and Julian Days
TODAYG=$(date -ud "${TODAYG}+1days" +'%Y-%m-%d') #> Add a day for tomorrow
TODAYJ=$(date -ud "${TODAYG}" +'%Y%j')
YESTERDAYG=$(date -ud "${TODAYG}-1days" +'%Y-%m-%d')
done
echo "Forward model for all days were finished!!"
