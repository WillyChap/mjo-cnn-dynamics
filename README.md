# Dynamical Signatures of a State-Dependent CNN Correction to MJO Propagation

Analysis and figure code for the Weather and Climate Dynamics (WCD) manuscript
*"Dynamical Signatures of a State-Dependent CNN Correction to Madden–Julian Oscillation
Propagation"* (Chapman and Berner). The study evaluates a CNN-based, state-dependent
model-error correction to CAM6 (coupled via FTorch) through process-oriented MJO diagnostics,
comparing a control (CNTRL), a climatological correction (MEAN), and the combined correction
(CNN) against ERA5. It is a follow-on to Chapman and Berner (2025, GRL, 10.1029/2024GL114106).

## Repository layout

```
figure_scripts/    Python figure builders + the metric/bootstrap analysis scripts
ncl/               NCL diagnostics that produce the regressed NetCDF products
run/               PBS drivers that run the NCL diagnostics per experiment
lib/               Python helpers (lanczos_filter.py, DIAH.py)
third_party/       Vendored Wheeler–Kiladis package (wk_spectra; see its LICENSE)
```

The manuscript source and the internal reproduction guide are kept outside this repository
while the paper is under preparation. This repository holds the analysis and figure code.

## Reproduce the figures

Figures are built by the scripts in `figure_scripts/` using the NCAR `npl-2024a` conda
environment:

```bash
module load conda
E=/glade/u/apps/opt/conda/envs/npl-2024a
conda run --no-capture-output -p $E python figure_scripts/make_panels.py blmc925 precip diah700 conindex u850
conda run --no-capture-output -p $E python figure_scripts/fig_spectra.py
conda run --no-capture-output -p $E python figure_scripts/fig_hovmoller.py
conda run --no-capture-output -p $E python figure_scripts/fig_heating_vstruct.py
conda run --no-capture-output -p $E python figure_scripts/fig_ept_vstruct.py
conda run --no-capture-output -p $E python figure_scripts/fig_condheat_append.py
conda run --no-capture-output -p $E python figure_scripts/fig_increment.py   # needs the INC_SUBSET file
```

Quantitative results: `compute_structural_metrics.py` (Table 3, pattern correlation and
normalized centered RMSE vs ERA5), `bootstrap_ew_season.py` (Table 2, the eastward:westward
spectral-power ratio and its 90% winter-block-bootstrap CI), `compute_kr_boxes.py` (the
fragility check on the Kelvin/Rossby amplitude ratio), and `compute_increment_metrics.py`
(Sect. 3.5, what the applied CNN wind increment does to the MJO; with only eleven winters it
reports an exact per-winter sign test and a leave-one-winter-out jackknife rather than relying
on a block bootstrap).

`fig_increment.py` is the one diagnostic of the correction itself rather than of the model's
response: it regresses the applied CNN wind tendency (`cb24cnn_U`, `cb24cnn_V`, archived only in
the fullCNN experiment) onto the MJO precipitation index. It reads the CNN experiment's
pressure-level timeseries directly, not the NCL products, and needs an `ncks` subset first:

```bash
module load nco
ncks -O -4 -L1 -v cb24cnn_U,cb24cnn_V,U,V,PRECT \
     -d lat,-25.0,25.0 -d lon,40.0,220.0 -d level,9,12 \
     /glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_fullCNN_DT2/ts/f.e.FTORCHmjo_fullCNN_DT2.cam.h1.NC4_Classic.plev.mjo.1979010100000-1990123100000.nc \
     inc_subset.nc
export INC_SUBSET=$PWD/inc_subset.nc
```

## Data

The large NetCDF products are not included; they are regenerated from source data on the
NCAR GLADE filesystem by the PBS drivers in `run/`, which export the data paths each NCL
diagnostic expects.
The model output, the trained CNN, and the FTorch coupling follow Chapman and Berner (2025);
ERA5 is from the ECMWF/Copernicus Climate Data Store; the prescribed SST/sea-ice boundary
condition is the standard CESM2/CAM6 merged dataset (Hurrell et al., 2008).

**Note on paths.** These scripts were run on NCAR HPC (Casper/Derecho) and contain absolute
GLADE paths. They document the exact analysis as performed; reproducing on another system
requires updating those paths to the local data.

## Third-party code

`third_party/wk_spectra/` is a vendored copy of the Wheeler–Kiladis spectral-analysis package
(from the `CMJO_Diagnostics_Tool`), used read-only to compute the wavenumber–frequency spectra.
See `third_party/wk_spectra/LICENSE` for its terms.

## Citation

Chapman, W. E. and Berner, J.: *Dynamical Signatures of a State-Dependent CNN Correction to
Madden–Julian Oscillation Propagation*, Weather and Climate Dynamics (in preparation).
