#!/usr/bin/env python
"""
Wheeler-Kiladis wavenumber-frequency spectra COMPARISON figure (MJO paper).

Computes the symmetric (and antisymmetric) background-removed WK spectra for four
experiments and assembles a single 2x2 comparison figure with a shared colorbar,
matsuno dispersion curves, and a marked MJO band:

    a) ERA5   (observations)
    b) CNTRL  (control)
    c) MITA1.0        (climatological / mean correction)
    d) MITA1.0CNN0.3  (CNN correction)

Environment:
    module load conda && conda activate npl-2024a

Run:
    python fig_spectra.py

The vendored wk_spectra package is used *read only* to compute the spectra; the
figure is hand-assembled here so all four panels share one figure and colorbar.
"""

import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle

# --- make the vendored package importable (do NOT modify it) --------------------
WK_PKG_DIR = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/CMJO_Diagnostics_Tool/wk_spectra"
if WK_PKG_DIR not in sys.path:
    sys.path.insert(0, WK_PKG_DIR)

import warnings
warnings.filterwarnings("ignore")

import wk_spectra.wk_analysis as wka
import wk_spectra.matsuno_plot as mp
from wk_spectra import nclcmaps

# --------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------
OUTDIR = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang/Paper_WCD"
FIGDIR = os.path.join(OUTDIR, "figures")
CACHE = os.path.join(FIGDIR, "fig_spectra_cache.npz")
os.makedirs(FIGDIR, exist_ok=True)

# panel key -> (label, FLUT file, latname)
# NOTE: the ERA5 scratch path in the driver notebook was purged; the surviving
# copy on campaign is the (annual-cycle-removed) anomaly file. The WK analysis
# removes the annual cycle and long-term trend internally and detrends/tapers
# each window, so using ERA5 FLUT anomalies is equivalent for the spectrum.
EXPERIMENTS = [
    ("ERA5",
     "ERA5",
     "/glade/campaign/cgd/amp/wchapman/ADF/ERA5_data/ts/ERA5.h1.FLUT.anomalies.1979010100000-1993123100000.nc"),
    ("CNTRL",
     "CNTRL",
     "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_CNTRLmjo_DT2/ts/f.e.FTORCHmjo_CNTRLmjo_DT2.cam.h1.FLUT.1979010100000-1990122200000.nc"),
    ("MEAN",
     "MEAN",
     "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_MEANmjo_DT2/ts/f.e.FTORCHmjo_MEANmjo_DT2.cam.h1.FLUT.1979010100000-1990122200000.nc"),
    ("CNN",
     "CNN",
     "/glade/campaign/cgd/amp/wchapman/ADF/f.e.FTORCHmjo_fullCNN_DT2/ts/f.e.FTORCHmjo_fullCNN_DT2.cam.h1.FLUT.1979010100000-1990123100000.nc"),
]

# WK parameters (identical to the per-experiment driver notebooks)
WK_KW = dict(spd=1, nDayWin=96, nDaySkip=30, max_freq=0.5, max_wn=20)
LAT_BOUND = 15

# Plot ranges
MAX_WN_PLOT = 15
MAX_FREQ_PLOT = 0.5
HE = [12, 25, 50]            # equivalent depths for dispersion curves
CPD_LINES = [3, 6, 30]       # period reference lines (days)

# Contour levels (shared across panels) -- from the driver notebooks
CLEVELS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.15, 1.2, 1.25,
           1.3, 1.35, 1.4, 1.45, 1.5, 1.6]


# --------------------------------------------------------------------------------
# Compute spectra (with caching so re-plotting is cheap)
# --------------------------------------------------------------------------------
def compute_all():
    results = {}
    for key, label, fpath in EXPERIMENTS:
        print(f"[compute] {key}: {fpath}", flush=True)
        b = wka.wk_analysis()
        b.import_netcdf(file=fpath, varname="FLUT", latname="lat", latBound=LAT_BOUND)
        b.wheeler_kiladis_spectra(**WK_KW)
        results[key] = dict(
            psumsym_r=b.wk_spectra["psumsym_r"],
            psumanti_r=b.wk_spectra["psumanti_r"],
            wavefft=b.wk_spectra["wavefft"],
            freqfft=b.wk_spectra["freqfft"],
        )
        del b
    return results


def load_or_compute():
    if os.path.exists(CACHE):
        print(f"[cache] loading {CACHE}", flush=True)
        z = np.load(CACHE, allow_pickle=True)
        return z["results"].item()
    results = compute_all()
    np.savez(CACHE, results=np.array(results, dtype=object))
    print(f"[cache] wrote {CACHE}", flush=True)
    return results


# --------------------------------------------------------------------------------
# Plotting helpers
# --------------------------------------------------------------------------------
def crop(spec):
    """Crop a spectrum dict to the plotting window and return X, Y, Z."""
    wave = spec["wavefft"]
    freq = spec["freqfft"]
    iw = np.where((wave >= -MAX_WN_PLOT - 0.1) & (wave <= MAX_WN_PLOT + 0.1))[0]
    ifq = np.where((freq > 0) & (freq <= MAX_FREQ_PLOT + 0.1))[0]
    w = wave[iw]
    f = freq[ifq]
    return np.meshgrid(w, f), iw, ifq


def draw_dispersion(ax, is_sym=True):
    modes = mp.matsuno_modes_wk(he=HE, n=[1], max_wn=MAX_WN_PLOT)
    for h in modes:
        df = modes[h]
        if is_sym:
            cols = [f"Kelvin(he={h}m)", f"ER(n=1,he={h}m)",
                    f"EIG(n=1,he={h}m)", f"WIG(n=1,he={h}m)"]
        else:
            cols = [f"MRG(he={h}m)", f"EIG(n=0,he={h}m)"]
        for c in cols:
            ax.plot(df.index.values, df[c].values, color="k",
                    linestyle="--", linewidth=0.7)


def draw_mjo_band(ax):
    """Mark the MJO band: zonal wavenumber 1-3, period 30-90 days (eastward)."""
    f_lo, f_hi = 1.0 / 90.0, 1.0 / 30.0   # 0.0111 - 0.0333 cpd
    rect = Rectangle((1, f_lo), 2, f_hi - f_lo, fill=False,
                     edgecolor="magenta", linewidth=2.2, zorder=6)
    ax.add_patch(rect)
    ax.text(3.15, f_hi + 0.004, "MJO", color="magenta", fontsize=15,
            fontweight="bold", zorder=6)


def style_axis(ax, XY, Z, cmap, norm, panel, title):
    (X, Y) = XY
    cset = ax.contourf(X, Y, Z, levels=CLEVELS, cmap=cmap, norm=norm, extend="max")
    lines = [l for l in CLEVELS if l >= 1.1]
    ax.contour(X, Y, Z, levels=lines, colors="k", linewidths=0.4)

    ax.axvline(x=0, color="k", linestyle="--", linewidth=0.6)

    # period reference lines
    for d in CPD_LINES:
        y = 1.0 / d
        if y <= MAX_FREQ_PLOT:
            ax.axhline(y=y, color="0.35", linestyle=":", linewidth=0.7)
            ax.text(-MAX_WN_PLOT + 0.3, y + 0.006, f"{d} d", size=12,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    ax.set_xlim(-MAX_WN_PLOT, MAX_WN_PLOT)
    ax.set_ylim(0.02, MAX_FREQ_PLOT)
    ax.tick_params(labelsize=14)
    ax.set_title(f"{panel}) {title}", fontsize=16, fontweight="bold", loc="left")
    return cset


def assemble(results, component="sym", outname=None):
    """component in {'sym','anti'}."""
    is_sym = component == "sym"
    key_r = "psumsym_r" if is_sym else "psumanti_r"

    cmap = nclcmaps.cmap("amwg_blueyellowred")
    norm = mpl.colors.BoundaryNorm(CLEVELS, cmap.N)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), sharex=True, sharey=True)
    panels = ["a", "b", "c", "d"]
    cset = None

    for ax, (key, label, _), pan in zip(axes.ravel(), EXPERIMENTS, panels):
        spec = results[key]
        XY, iw, ifq = crop(spec)
        Z = spec[key_r][np.ix_(ifq, iw)]
        cset = style_axis(ax, XY, Z, cmap, norm, pan, label)
        draw_dispersion(ax, is_sym=is_sym)
        draw_mjo_band(ax)

    # eastward / westward annotation on bottom row
    for ax in axes[1, :]:
        ax.set_xlabel("Zonal wavenumber", fontsize=16, fontweight="bold")
        ax.text(MAX_WN_PLOT * 0.45, -0.055, "EASTWARD", fontsize=12,
                fontweight="bold", ha="center")
        ax.text(-MAX_WN_PLOT * 0.55, -0.055, "WESTWARD", fontsize=12,
                fontweight="bold", ha="center")
    for ax in axes[:, 0]:
        ax.set_ylabel("Frequency (CPD)", fontsize=16, fontweight="bold")

    comp_name = "Symmetric" if is_sym else "Antisymmetric"
    # no figure title: the manuscript caption serves this role

    fig.tight_layout(rect=[0, 0, 0.90, 0.96])
    cax = fig.add_axes([0.92, 0.12, 0.02, 0.74])
    cb = fig.colorbar(cset, cax=cax, ticks=CLEVELS, boundaries=CLEVELS,
                      spacing="uniform", extend="max")
    cb.set_label("Power / background", fontsize=15)
    cb.ax.tick_params(labelsize=13)

    if outname is None:
        outname = os.path.join(
            FIGDIR, f"fig_spectra_{'symmetric' if is_sym else 'antisymmetric'}.png")
    fig.savefig(outname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {outname}", flush=True)
    return outname


def main():
    results = load_or_compute()
    sym = assemble(results, "sym")
    anti = assemble(results, "anti")
    for p in (sym, anti):
        sz = os.path.getsize(p)
        print(f"[done] {p}  ({sz/1e6:.2f} MB)", flush=True)


if __name__ == "__main__":
    main()
