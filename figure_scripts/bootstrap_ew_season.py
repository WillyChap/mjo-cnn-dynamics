#!/usr/bin/env python
"""
Block-bootstrap confidence intervals for the MJO-band eastward:westward spectral
power ratio, resampling INDEPENDENT seasonal (boreal-winter) blocks.

Motivation
----------
bootstrap_ew_ratio.py resampled the individual ~96-day WK analysis windows with
replacement.  A reviewer correctly noted that those windows are NOT mutually
independent: consecutive windows sit inside the same boreal cold season / year and
share the same interannual state (ENSO, background MJO activity), so resampling
them individually treats correlated windows as independent and makes the interval
too narrow.

Fix
---
The natural independent block is a whole boreal-winter-centred year.  The WK
pipeline builds non-overlapping 96-day windows spaced (nSampWin+nSampSkip)=126
days apart (spd=1, nDayWin=96, nDaySkip=30).  Each window is assigned to a
winter-year block from its centre date:

    block = year      if centre-month >= 7   (Jul(Y)..Jun(Y+1) -> centred on DJF)
    block = year - 1  otherwise

so every block spans one boreal cold season plus the surrounding shoulder seasons
and contains ~3 consecutive windows.  We resample the blocks (not the windows)
with replacement, keeping each block's windows together, and rebuild the WK
reduction (symmetric power / smoothed red-noise background) on the resampled
window set exactly as the point estimate does.

For each experiment:
  1. replicate the WK pipeline up to power[window, freq, lat, lon]
     (reusing the vendored wk_spectra functions unmodified);
  2. assign every window to a winter-year block from its centre date;
  3. point estimate: average over ALL windows -> psumsym / psumb -> ratio,
     integrate eastward(wn 1-3)/westward(wn -1..-3) over periods 30-90 d;
  4. bootstrap: draw nBlocks blocks WITH replacement, gather the windows of the
     drawn blocks, rebuild psumsym/psumb on that window set, integrate E:W ratio;
     repeat NBOOT times -> 5th/95th percentile = 90% CI.

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
    """Return power[nWindow,...], nlat, nSampWin, wavefft, freqfft, block[nWindow]."""
    f = ncf.Dataset(fpath, "r")
    data = f.variables["FLUT"][:]
    lats = f.variables["lat"][:]
    t = f.variables["time"]
    dts = ncf.num2date(t[:], t.units, getattr(t, "calendar", "standard"))
    f.close()
    li = np.where((lats >= -LAT_BOUND) & (lats <= LAT_BOUND))[0]
    data = np.asarray(data[:, li, :])

    spd, nDayWin, nDaySkip = WK_KW["spd"], WK_KW["nDayWin"], WK_KW["nDaySkip"]
    ntim, nlat, nlon = data.shape
    nSampWin = nDayWin * spd
    nSampSkip = nDaySkip * spd
    step = nSampWin + nSampSkip
    nWindow = (ntim - nSampWin) // step + 1

    # winter-year block label from each window's centre date
    blocks = np.empty(nWindow, dtype=int)
    for nw in range(nWindow):
        ci = nw * step + nSampWin // 2
        d = dts[ci]
        blocks[nw] = d.year if d.month >= 7 else d.year - 1

    array_dt = wka.remove_dominant_signals(data, spd, nDayWin, nDaySkip)
    array_as = wka.decompose_symasym(array_dt)
    wavefft, freqfft, peeAS = wka.spectral_coefficients(array_as, spd, nDayWin, nDaySkip)
    power = (np.abs(peeAS)) ** 2  # [nWindow, freq, lat, lon]
    return power, nlat, nSampWin, wavefft, freqfft, blocks


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
        power, nlat, nSampWin, wavefft, freqfft, blocks = per_window_power(fpath)
        nWindow = power.shape[0]

        # group window indices by block
        ublocks = np.unique(blocks)
        nBlocks = len(ublocks)
        win_by_block = {b: np.where(blocks == b)[0] for b in ublocks}

        point = ew_ratio_from_power(power, nlat, nSampWin, wavefft, freqfft)

        boots = np.empty(NBOOT)
        for b in range(NBOOT):
            drawn = ublocks[rng.integers(0, nBlocks, size=nBlocks)]
            idx = np.concatenate([win_by_block[bb] for bb in drawn])
            boots[b] = ew_ratio_from_power(power[idx], nlat, nSampWin, wavefft, freqfft)

        lo, hi = np.percentile(boots, [5, 95])
        summary[key] = dict(point=point, lo=lo, hi=hi, mean=boots.mean(),
                            std=boots.std(ddof=1), nWindow=nWindow,
                            nBlocks=nBlocks, boots=boots)
        print(f"[{key:6s}] nBlocks={nBlocks} nWindow={nWindow}  point={point:.3f} "
              f"(paper {PAPER[key]})  boot mean={boots.mean():.3f}  "
              f"90% CI=[{lo:.3f}, {hi:.3f}]  sd={boots.std(ddof=1):.3f}  "
              f"({time.time()-t0:.1f}s)", flush=True)

    # CI-overlap check vs CNN
    print("\n=== 90% CI overlap with CNN [{:.3f}, {:.3f}] ===".format(
        summary["CNN"]["lo"], summary["CNN"]["hi"]))
    cnn = summary["CNN"]
    for k in ("CNTRL", "MEAN", "ERA5"):
        s = summary[k]
        overlap = (s["lo"] <= cnn["hi"]) and (cnn["lo"] <= s["hi"])
        print(f"  CNN vs {k:5s} [{s['lo']:.3f}, {s['hi']:.3f}]: "
              f"{'OVERLAP' if overlap else 'DISJOINT'}")

    outpath = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang/Paper_WCD/figures/bootstrap_ew_season.npz"
    np.savez(outpath,
             **{k + "_boots": summary[k]["boots"] for k in summary},
             point=np.array([summary[k]["point"] for k in summary]),
             lo=np.array([summary[k]["lo"] for k in summary]),
             hi=np.array([summary[k]["hi"] for k in summary]),
             nBlocks=np.array([summary[k]["nBlocks"] for k in summary]),
             nWindow=np.array([summary[k]["nWindow"] for k in summary]),
             keys=np.array(list(summary.keys())))
    print(f"\nsaved boots to {outpath}")
    return summary


if __name__ == "__main__":
    main()
