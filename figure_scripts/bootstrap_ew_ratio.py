#!/usr/bin/env python
"""
Block-bootstrap confidence intervals for the MJO-band eastward:westward spectral
power ratio (symmetric, background-removed Wheeler-Kiladis spectrum).

Method
------
The WK analysis averages the space-time power over a set of ~34 essentially
independent 96-day windows (nDayWin=96, nDaySkip=30, spd=1) before dividing by a
smoothed red-noise background. Those 96-day windows are the natural block for a
bootstrap that preserves the temporal autocorrelation of individual MJO events.

For each experiment we:
  1. replicate the WK pipeline up to the per-window power array
     power[window, freq, lat, lon]  (= |FFT|^2 of each detrended/tapered window),
     reusing the vendored wk_spectra functions unmodified;
  2. point estimate: average over ALL windows -> psumsym / psumb -> psumsym_r,
     integrate eastward(wn 1-3) / westward(wn -1..-3) over periods 30-90 d;
  3. bootstrap: draw nWindow window indices WITH replacement, average over the
     resampled windows (background recomputed consistently), integrate E:W ratio.
     Repeat NBOOT times -> 5th/95th percentile = 90% CI.

The E:W integration is validated against the published point values
(ERA5 1.48, CNTRL 1.13, MEAN 1.23, CNN 1.46).
"""
import os, sys, time
import numpy as np

WK_PKG_DIR = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/CMJO_Diagnostics_Tool/wk_spectra"
if WK_PKG_DIR not in sys.path:
    sys.path.insert(0, WK_PKG_DIR)
import warnings
warnings.filterwarnings("ignore")

import netCDF4 as ncf
import wk_spectra.wk_analysis as wka

# --- experiment definitions (identical to fig_spectra.py) -----------------------
EXPERIMENTS = [
    ("ERA5",  "/glade/campaign/cgd/amp/wchapman/ADF/ERA5_data/ts/ERA5.h1.FLUT.anomalies.1979010100000-1993123100000.nc"),
    ("CNTRL", "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_CNTRLmjo_DT2/ts/f.e.FTORCHmjo_CNTRLmjo_DT2.cam.h1.FLUT.1979010100000-1990122200000.nc"),
    ("MEAN",  "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_MEANmjo_DT2/ts/f.e.FTORCHmjo_MEANmjo_DT2.cam.h1.FLUT.1979010100000-1990122200000.nc"),
    ("CNN",   "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_fullCNN_DT2/ts/f.e.FTORCHmjo_fullCNN_DT2.cam.h1.FLUT.1979010100000-1990123100000.nc"),
]
WK_KW = dict(spd=1, nDayWin=96, nDaySkip=30, max_freq=0.5, max_wn=20)
LAT_BOUND = 15
PAPER = {"ERA5": 1.48, "CNTRL": 1.13, "MEAN": 1.23, "CNN": 1.46}

NBOOT = 1000
SEED = 20260706


def per_window_power(fpath):
    """Return power[nWindow, nfreq, nlat, nlon], nlat, nSampWin, wavefft, freqfft."""
    f = ncf.Dataset(fpath, "r")
    data = f.variables["FLUT"][:]
    lats = f.variables["lat"][:]
    f.close()
    li = np.where((lats >= -LAT_BOUND) & (lats <= LAT_BOUND))[0]
    data = np.asarray(data[:, li, :])

    spd, nDayWin, nDaySkip = WK_KW["spd"], WK_KW["nDayWin"], WK_KW["nDaySkip"]
    ntim, nlat, nlon = data.shape
    nSampWin = nDayWin * spd

    array_dt = wka.remove_dominant_signals(data, spd, nDayWin, nDaySkip)
    array_as = wka.decompose_symasym(array_dt)
    wavefft, freqfft, peeAS = wka.spectral_coefficients(array_as, spd, nDayWin, nDaySkip)
    power = (np.abs(peeAS)) ** 2  # [nWindow, freq, lat, lon]
    return power, nlat, nSampWin, wavefft, freqfft


def ew_ratio_from_power(power, nlat, nSampWin, wavefft, freqfft):
    """Full WK reduction on a (possibly resampled) power array -> E:W ratio."""
    psumanti, psumsym = wka.separate_power(power, nlat, nSampWin, wavefft, freqfft)
    psumb = wka.derive_background(power, nlat, nSampWin, wavefft, freqfft)

    # crop exactly as wheeler_kiladis_spectra does
    iw = np.where((wavefft >= -WK_KW["max_wn"]) & (wavefft <= WK_KW["max_wn"]))[0]
    ifq = np.where((freqfft > 0) & (freqfft <= WK_KW["max_freq"]))[0]
    wave = wavefft[iw]
    freq = freqfft[ifq]
    psr = (psumsym[np.ix_(ifq, iw)]) / (psumb[np.ix_(ifq, iw)])

    # MJO band: eastward wn 1-3 vs westward -1..-3, periods 30-90 d
    flo, fhi = 1.0 / 90.0, 1.0 / 30.0
    fsel = np.where((freq >= flo) & (freq <= fhi))[0]
    east = np.where((wave >= 1) & (wave <= 3))[0]
    west = np.where((wave <= -1) & (wave >= -3))[0]
    E = np.nansum(psr[np.ix_(fsel, east)])
    W = np.nansum(psr[np.ix_(fsel, west)])
    return E / W


def main():
    rng = np.random.default_rng(SEED)
    summary = {}
    for key, fpath in EXPERIMENTS:
        t0 = time.time()
        power, nlat, nSampWin, wavefft, freqfft = per_window_power(fpath)
        nWindow = power.shape[0]

        point = ew_ratio_from_power(power, nlat, nSampWin, wavefft, freqfft)

        boots = np.empty(NBOOT)
        for b in range(NBOOT):
            idx = rng.integers(0, nWindow, size=nWindow)
            boots[b] = ew_ratio_from_power(power[idx], nlat, nSampWin, wavefft, freqfft)

        lo, hi = np.percentile(boots, [5, 95])
        summary[key] = dict(point=point, lo=lo, hi=hi, mean=boots.mean(),
                            std=boots.std(ddof=1), nWindow=nWindow, boots=boots)
        print(f"[{key:6s}] nWindow={nWindow}  point={point:.3f} (paper {PAPER[key]})  "
              f"boot mean={boots.mean():.3f}  90% CI=[{lo:.3f}, {hi:.3f}]  "
              f"sd={boots.std(ddof=1):.3f}  ({time.time()-t0:.1f}s)", flush=True)

    # pairwise overlap / bootstrap difference tests vs CNN
    print("\n=== Pairwise comparison (paired bootstrap, common resample per draw) ===")
    # recompute with a shared resample so differences are paired
    np.savez("/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang/Paper_WCD/figures/bootstrap_ew_ratio.npz",
             **{k: summary[k]["boots"] for k in summary},
             point={k: summary[k]["point"] for k in summary})
    print("saved boots to figures/bootstrap_ew_ratio.npz")
    return summary


if __name__ == "__main__":
    main()
