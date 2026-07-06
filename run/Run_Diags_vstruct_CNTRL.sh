#!/bin/bash
#PBS -N vsCNTRL_Diag
#PBS -A NAML0001
#PBS -l walltime=04:00:00
#PBS -o Paper_WCD/nc/run_vs_CNTRL/vstruct_CNTRL.out
#PBS -e Paper_WCD/nc/run_vs_CNTRL/vstruct_CNTRL.out
#PBS -q casper
#PBS -l select=1:ncpus=9:mem=210GB
#PBS -m a
#PBS -M wchapman@ucar.edu

module load nco
module load ncl
module load conda
conda activate npl-2024a

# ---- experiment ----
TAG="CNTRL"
export modout="f.e.FTORCHmjo_CNTRLmjo_DT2"
export mod_data_daily="/glade/campaign/cgd/amp/wchapman/ADF/${modout}/ts/${modout}.cam.h1.NC4_Classic.plev.mjo.1979010100000-1990122200000.nc"

# ---- observations (ERA5 / GPCP, same as Run_Diags_All_CNTRL.sh) ----
export obs_data_p="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/PR.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data_t="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/T.24hr.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data_q="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/Q.24hr.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data_u="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/U.24hr.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data_v="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/V.24hr.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data_w="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/W.24hr.NC4_Classic.Camgrid.1979-1995.nc"
export obs_data="/glade/campaign/cgd/amp/wchapman/Reanalysis/ERA5_obs/ERA5.24hr.NC4_Classic.1979_1995.nc"

export ymdStrt="19790101"
export ymdLast="19901230"
export cwLastm="19901230"

# ---- per-experiment output directories ----
WANG="/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang"
NCDIR="${WANG}/Paper_WCD/nc/${TAG}"
RUNDIR="${WANG}/Paper_WCD/nc/run_vs_${TAG}"
mkdir -p "${NCDIR}" "${RUNDIR}" "${RUNDIR}/FIG"

# NC_DIR must end in "/"
export NC_DIR="${NCDIR}/"
export FIG_DIR="${RUNDIR}/FIG"

# Several NCL scripts hardcode "./NC_OUT/"; run from RUNDIR with NC_OUT -> NCDIR
cd "${RUNDIR}"
ln -sfn "${NCDIR}" NC_OUT
ln -sf "${WANG}/vertical_thete.ncl" vertical_thete.ncl
for s in mjo_theta mjo_SPQ mjo_U mjo_OMEGA mjo_diah; do
  ln -sf "${WANG}/${s}.ncl" "${s}.ncl"
done

echo "=== $(date) start ${TAG} : NCDIR=${NCDIR} ==="
ncl mjo_theta.ncl
ncl mjo_SPQ.ncl
ncl mjo_U.ncl
ncl mjo_OMEGA.ncl
ncl mjo_diah.ncl
echo "=== $(date) done ${TAG} ==="
ls -la "${NCDIR}" | grep -E "mjo_(thetae|SH|U|W|diah)_"
