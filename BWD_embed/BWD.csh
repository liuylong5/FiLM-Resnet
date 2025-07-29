module purge
module load oneapi 
export LD_LIBRARY_PATH=../public_model/Libs/netcdf/lib:${LD_LIBRARY_PATH} 
export MODEL_PATH=../BWD/bwd_12km_2020_embed/path
export LD_LIBRARY_PATH=../BWD/bwd_12km_2020_embed:$LD_LIBRARY_PATH
EXEC=../BWD/bwd_12km_2020_embed/bwd_intel 
work_dir=../BWD/bwd_12km_2020_embed 
export MY_LOG_FILE=../BWD/bwd_12km_2020_embed/BWD.log
FWDDIR=../FWD/2020_12km_embed
npcol=4
nprow=4
OUTDIR=../BWD/bwd_12km_2020_embed 
START_DATE=$(date -d "2020-01-01" +'%Y-%m-%d') 
END_DATE=$(date -d "2020-01-01" +'%Y-%m-%d') 

let NPROCS=npcol*nprow 
START_J=$(date -ud "${START_DATE}" +'%Y%j') 
END_J=$(date -ud "${END_DATE}" +'%Y%j') 
export NPCOL_NPROW="$npcol $nprow" 
export ADJ_LGRID_FREQ="SYNC_STEP" 
export GRID_NAME=China_12km
export GRIDDESC=../FWD/2020_12km/GRIDDESC 
export CTM_MAXSYNC=240
export CTM_MINSYNC=240
export KZMIN=Y  
export IOAPI_LOG_WRITE=F 
export FL_ERR_STOP=F  
export CTM_APPL='bwd_intel' 
export FLOOR_FILE=$work_dir/floor 

TODAYG=$(date -d "$START_DATE" +%Y-%m-%d)
TODAYJ=$(date -ud "${TODAYG}" +'%Y%j')
while (( "${TODAYJ}" >= "${END_J}" )); do
echo "Backward model of $TODAYG was started!"
NEXTDAYG=$(date -ud "${TODAYG}+1days" +'%Y-%m-%d')
if [ "$TODAYJ" != "$START_J" ]; then
export INIT_LGRID_1=$OUTDIR/lgrid_$NEXTDAYG.nc
if [ ! -f "$INIT_LGRID_1" ]; then
echo " Error: Initial LGRID file, "$INIT_LGRID_1" is invalid"
exit 1
fi
fi
export EMIS_1="$FWDDIR/emission_${TODAYG}.nc -v"
export BNDY_GASC_1="" #"$FWDDIR/BCON_V5f_${TODAYG} -v"
export BNDY_AERO_1="" #"$FWDDIR/BCON_V5f_${TODAYG} -v"
export BNDY_NONR_1="" #"$FWDDIR/BCON_V5f_${TODAYG} -v"
export BNDY_TRAC_1="" #"$FWDDIR/BCON_V5f_${TODAYG} -v"
export GRID_DOT_2D="$FWDDIR/GRIDDOT2D_${TODAYG}.nc -v"
export GRID_CRO_2D="$FWDDIR/GRIDCRO2D_${TODAYG}.nc -v"
export MET_CRO_2D="$FWDDIR/METCRO2D_${TODAYG}.nc -v"
export MET_DOT_3D="$FWDDIR/METDOT3D_${TODAYG}.nc -v"
export MET_CRO_3D="$FWDDIR/METCRO3D_${TODAYG}.nc -v"
export MET_BDY_3D="$FWDDIR/METBDY3D_${TODAYG}.nc -v"
export LAYER_FILE="$FWDDIR/METCRO3D_${TODAYG}.nc -v"
export CTM_CONC_1="$FWDDIR/CONC_${TODAYG}.nc -v"
export A_CONC_1="$FWDDIR/ACONC_${TODAYG}.nc -v"
export S_CGRID="$FWDDIR/CGRID_${TODAYG}.nc -v"
export CTM_DRY_DEP_1="$FWDDIR/DRYDEP_${TODAYG}.nc -v"
export CTM_WET_DEP_1="$FWDDIR/WETDEP1_${TODAYG}.nc -v"
export CTM_WET_DEP_2="$FWDDIR/WETDEP2_${TODAYG}.nc -v"
export CTM_SSEMIS_1="$FWDDIR/SSEMIS1_${TODAYG}.nc -v"
export CTM_VIS_1="$FWDDIR/AEROVIS_${TODAYG}.nc -v"
export CTM_DIAM_1="$FWDDIR/AERODIAM_${TODAYG}.nc -v"
export CTM_IPR_1="$FWDDIR/PA_1_${TODAYG}.nc -v"
export CTM_IPR_2="$FWDDIR/PA_2_${TODAYG}.nc -v"
export CTM_IPR_3="$FWDDIR/PA_3_${TODAYG}.nc -v"
export CTM_IRR_1="$FWDDIR/IRR_1_${TODAYG}.nc -v"
export CTM_IRR_2="$FWDDIR/IRR_2_${TODAYG}.nc -v"
export CTM_IRR_3="$FWDDIR/IRR_3_${TODAYG}.nc -v"
export CTM_RJ_1="$FWDDIR/RJ_1_${TODAYG}.nc -v"
export CTM_RJ_2="$FWDDIR/RJ_2_${TODAYG}.nc -v"
export DBG_IGRID="$FWDDIR/DBG_iCONVCLD_${TODAYG} -v"
export DBG_RGRID="$FWDDIR/DBG_rCONVCLD_${TODAYG} -v"
export CTM_CONC_FWD="$FWDDIR/CONC_${TODAYG}.nc -v"
export ADJ_CHEM_CHK="$FWDDIR/CHEM_STA_${TODAYG}.nc -v"
export ADJ_VDIFF_CHK="$FWDDIR/VDIFF_CHK_${TODAYG}.nc -v"
export ADJ_HA_RHOJ_CHK="$FWDDIR/HA_RHOJ_CHK_${TODAYG}.nc -v"
export ADJ_VA_RHOJ_CHK="$FWDDIR/VA_RHOJ_CHK_${TODAYG}.nc -v"
export ADJ_AERO_CHK="$FWDDIR/AERO_CHK_${TODAYG}.nc -v"
export ADJ_HADV_CHK="$FWDDIR/HADV_CHK_${TODAYG}.nc -v"
export ADJ_VADV_CHK="$FWDDIR/VADV_CHK_${TODAYG}.nc -v"
export ADJ_CLD_CHK="$FWDDIR/CLD_CHK_${TODAYG}.nc -v"
export ADJ_EMIS_CHK="$FWDDIR/EMIS_CHK_${TODAYG}.nc -v"
export ADJ_EMISU_CHK="$FWDDIR/EMISU_CHK_${TODAYG}.nc -v"
export CTM_XFIRST_IN="$FWDDIR/XFIRST_${TODAYG}"
export ADJ_LGRID="$OUTDIR/lgrid_${TODAYG}.nc"  #adjoint sens. output
export ADJ_LGRID_EM="$OUTDIR/lgrid_em_${TODAYG}.nc"  #adjoint sens. output
export CTM_STDATE=$TODAYJ
export CTM_TSTEP=010000
export CTM_STTIME=000000
export CTM_RUNLEN=240000
export ADJ_FRC=Y 
export ADJ_FORCE="ADJ_FORCE_${TODAYG}.nc -v" 

cd $work_dir
/bin/rm $ADJ_LGRID >&/dev/null
rm CTM* >&/dev/null
{ time -p mpirun -np $NPROCS $EXEC; } > ${MY_LOG_FILE}
if [ $? -ne 0 ]; then
echo "Model run failed for $TODAYG"
exit 1  # Exit the script if the model run fails
fi
mv FLOOR_* ${FLOOR_FILE} >&/dev/null
echo "Backward model of ${TODAYG} was finished!"

#> Increment both Gregorian and Julian Days
TODAYG=$(date -ud "${TODAYG}-1days" +'%Y-%m-%d') #> Add a day for tomorrow
TODAYJ=$(date -ud "${TODAYG}" +'%Y%j')
done
echo "Backward model for all days were finished!!"
