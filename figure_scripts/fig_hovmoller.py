#!/usr/bin/env python
"""
MJO propagation figure (lag-longitude / Hovmoller correlation).

CLIVAR MJO diagnostic: lag-longitude correlation of 20-100 day band-pass
filtered precipitation averaged 10S-10N, correlated against the band-pass
filtered precipitation index over the equatorial Indian Ocean
(10S-5N, 75-100E), with 850 hPa zonal-wind correlation overlaid as contours.
Eastward MJO propagation across the Maritime Continent appears as a positive
slope with increasing lag.

Four panels:
    a) ERA5   (OBS)
    b) CNTRL  (f.e.FTORCHmjo_CNTRLmjo_DT2)
    c) MEAN (f.e.FTORCHmjo_MEANmjo_DT2, climatological/MEAN correction)
    d) CNN (f.e.FTORCHmjo_fullCNN_DT2, CNN correction)

Adapted from CMJO_Diagnostics_Tool/09_mjoxcor_lag_season.py and
Panel_XLAG_FIGURE.ipynb (the extra CNNmjo experiment is dropped here).

Period: 1979-1990.  Anomalies computed on the fly (daily climatology removed)
from the raw daily timeseries, then band-pass filtered 20-100 days.
"""

import os
import copy
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import ticker
from cartopy.mpl.ticker import LongitudeFormatter

import warnings
warnings.filterwarnings("ignore")


# ----------------------------------------------------------------------------
# Filter / correlation helpers (copied verbatim from
# CMJO_Diagnostics_Tool/py so this script is self-contained)
# ----------------------------------------------------------------------------
def bandpass_filter_via_fft_1d(time_series, low_period, high_period, sampling_interval=1):
    """Band-pass a 1-D series via FFT. Keeps periods in [high_period, low_period]."""
    fft_data = np.fft.fft(time_series)
    n = len(fft_data)
    frequencies = np.fft.fftfreq(n, d=sampling_interval)
    f_low = 1 / high_period
    f_high = 1 / low_period
    mask = (np.abs(frequencies) >= f_low) & (np.abs(frequencies) <= f_high)
    return np.fft.ifft(fft_data * mask).real


def bandpass_filter_via_fft_2d(data, low_period, high_period, sampling_interval=1):
    """Band-pass each column (time, space) via FFT."""
    filtered = np.empty_like(data, dtype=float)
    n = data.shape[0]
    frequencies = np.fft.fftfreq(n, d=sampling_interval)
    f_low = 1 / high_period
    f_high = 1 / low_period
    mask = (np.abs(frequencies) >= f_low) & (np.abs(frequencies) <= f_high)
    for i in range(data.shape[1]):
        fft_data = np.fft.fft(data[:, i])
        filtered[:, i] = np.fft.ifft(fft_data * mask).real
    return filtered


def lagged_correlation(x, y, maxlag):
    """Lagged correlations between two pandas.Series (as in py)."""
    correlations = {}
    for lag in range(-maxlag, maxlag + 1):
        if lag >= 0:
            correlations[lag] = x.corr(y.shift(lag))
        else:
            correlations[lag] = x.shift(-lag).corr(y)
    return correlations

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
OUTDIR = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang/Paper_WCD"
FIG_OUT = os.path.join(OUTDIR, "figures", "fig_hovmoller.png")

# Region / filter settings (same as the CMJO diagnostic template)
LAT_IO = slice(-10, 5)      # Indian-Ocean index latitude band
LON_IO = slice(75, 100)     # Indian-Ocean index longitude band
LAT_NS = slice(-10, 10)     # equatorial band for the Hovmoller average
F_LOW = 100                 # period (days) - low-frequency cutoff
F_HIGH = 20                 # period (days) - high-frequency cutoff
MAXLAG = 30                 # +/- lags (days)
T0, T1 = "1979-01-01", "1990-12-31"   # analysis period

# Experiments: (label, panel_letter, precip_file, precip_var, u_file, u_var)
ADF_C = "/glade/campaign/cgd/amp/wchapman/ADF"
EXPS = [
    ("a) ERA5", "a",
     f"{ADF_C}/ERA5_data/ts/ERA5.h1.PR.1979010100000-1993123100000.nc", "PR",
     f"{ADF_C}/ERA5_data/ts/ERA5.h1.U850.anomalies.1979010100000-1993123100000.nc", "U850"),
    ("b) CNTRL", "b",
     f"{ADF_C}/f.e.FTORCHmjo_CNTRLmjo_DT2/ts/f.e.FTORCHmjo_CNTRLmjo_DT2.cam.h1.PRECT.1979010100000-1990122200000.nc", "PRECT",
     f"{ADF_C}/f.e.FTORCHmjo_CNTRLmjo_DT2/ts/f.e.FTORCHmjo_CNTRLmjo_DT2.cam.h1.U850.1979010100000-1990122200000.nc", "U850"),
    ("c) MEAN", "c",
     f"{ADF_C}/f.e.FTORCHmjo_MEANmjo_DT2/ts/f.e.FTORCHmjo_MEANmjo_DT2.cam.h1.PRECT.1979010100000-1990122200000.nc", "PRECT",
     f"{ADF_C}/f.e.FTORCHmjo_MEANmjo_DT2/ts/f.e.FTORCHmjo_MEANmjo_DT2.cam.h1.U850.1979010100000-1990122200000.nc", "U850"),
    ("d) CNN", "d",
     f"{ADF_C}/f.e.FTORCHmjo_fullCNN_DT2/ts/f.e.FTORCHmjo_fullCNN_DT2.cam.h1.PRECT.1979010100000-1990123100000.nc", "PRECT",
     f"{ADF_C}/f.e.FTORCHmjo_fullCNN_DT2/ts/f.e.FTORCHmjo_fullCNN_DT2.cam.h1.U850.1979010100000-1990123100000.nc", "U850"),
]


def remove_daily_climatology(da):
    """Remove the long-term daily (seasonal-cycle) climatology."""
    clim = da.groupby("time.dayofyear").mean("time")
    return da.groupby("time.dayofyear") - clim


def compute_lagged_corr_lon(precip_file, precip_var, u_file, u_var):
    """Compute the lag-longitude correlation of band-pass filtered precip and
    U850 against the Indian-Ocean precip index (boreal-winter seasons averaged).
    Returns (lons, lags, xcorr_P[lag,lon], xcorr_U[lag,lon])."""
    print(f"  precip: {os.path.basename(precip_file)} ({precip_var})")
    print(f"  u850  : {os.path.basename(u_file)} ({u_var})")

    dsp = xr.open_dataset(precip_file)[[precip_var]].sel(time=slice(T0, T1))
    dsu = xr.open_dataset(u_file)[[u_var]].sel(time=slice(T0, T1))

    # normalize to a plain daily datetime index (CFTime / cftime tolerant)
    dsp = dsp.sortby("time")
    dsu = dsu.sortby("time")

    weights = np.cos(np.deg2rad(dsp.lat))

    # remove daily climatology (seasonal cycle) -> anomalies
    p_anom = remove_daily_climatology(dsp[precip_var])
    u_anom = remove_daily_climatology(dsu[u_var])

    # ----- Indian Ocean precip index -----
    PIO = (p_anom * weights).sel(lat=LAT_IO, lon=LON_IO).mean(["lat", "lon"])
    PIO_arr = np.array(PIO)
    PIO_arr[np.isnan(PIO_arr)] = np.nanmean(PIO_arr)
    PIO_f = bandpass_filter_via_fft_1d(PIO_arr, F_HIGH, F_LOW, sampling_interval=1)
    PIO_f = xr.DataArray(PIO_f, coords={"time": PIO.time}, dims="time")

    # ----- equatorial-band longitude average, band-pass filtered -----
    P_tl = (p_anom * weights).sel(lat=LAT_NS).mean("lat")   # (time, lon)
    U_tl = (u_anom * weights).sel(lat=LAT_NS).mean("lat")   # (time, lon)
    P_tl_f = P_tl.copy(data=bandpass_filter_via_fft_2d(
        np.array(P_tl), F_HIGH, F_LOW, sampling_interval=1))
    U_tl_f = U_tl.copy(data=bandpass_filter_via_fft_2d(
        np.array(U_tl), F_HIGH, F_LOW, sampling_interval=1))

    lons = dsp["lon"].values
    nlon = lons.size

    years = np.unique(P_tl_f["time.year"].values)
    seasons = years[:-2]   # drop last 2 (incomplete Nov-May window), as template
    nlag = 2 * MAXLAG + 1

    xc_P = np.zeros([len(seasons), nlon, nlag])
    xc_U = np.zeros([len(seasons), nlon, nlag])

    lags = None
    for ss, sd in enumerate(seasons):
        tsl = slice(f"{sd}-11-01", f"{sd + 1}-05-31")
        x = pd.Series(np.array(PIO_f.sel(time=tsl)))
        Pseas = np.array(P_tl_f.sel(time=tsl))   # (t, lon)
        Useas = np.array(U_tl_f.sel(time=tsl))
        for vv in range(nlon):
            yP = pd.Series(Pseas[:, vv])
            cP = lagged_correlation(yP, x, MAXLAG)
            xc_P[ss, vv, :] = np.array(list(cP.values()))
            yU = pd.Series(Useas[:, vv])
            cU = lagged_correlation(yU, x, MAXLAG)
            xc_U[ss, vv, :] = np.array(list(cU.values()))
        if lags is None:
            lags = np.array(list(cP.keys()))

    meanP = np.nanmean(xc_P, axis=0).T   # (lag, lon)
    meanU = np.nanmean(xc_U, axis=0).T
    return lons, lags, meanP, meanU


# ----------------------------------------------------------------------------
# Colormap (matches the CMJO diagnostic template)
# ----------------------------------------------------------------------------
def build_cmap():
    cmap = plt.cm.RdYlBu_r
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = cmap(1)
    cmaplist[cmap.N - 1] = cmap(0.99)
    for ii in range(120, 136):
        cmaplist[ii] = [1, 1, 1, 1]
    cmap = cmap.from_list("My cmap", cmaplist, cmap.N)
    cmap.set_under([0.3, 0.00, 0.1, 1.0])
    cmap.set_over("k")
    return cmap


def main():
    os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True)

    results = []
    for label, letter, pf, pv, uf, uv in EXPS:
        print(f"=== {label} ===")
        lons, lags, mP, mU = compute_lagged_corr_lon(pf, pv, uf, uv)
        results.append((label, letter, lons, lags, mP, mU))

    cmap = build_cmap()
    clevels = np.arange(-1.0, 1.05, 0.05)
    contlevels = np.arange(-1, 1.1, 0.1)
    contlevels_emph = [-0.3, 0.3]
    norm = mpl.colors.BoundaryNorm(clevels, cmap.N)

    # 5 m/s eastward phase-speed reference line
    # 5 m/s = 432 km/day; at the equator ~111.32 km/deg -> ~3.88 deg/day
    deg_per_day = 5.0 * 86400.0 / 1000.0 / 111.32
    ref_lag = np.array([0, 20])
    ref_lon0 = 80.0                       # anchor at IO longitude, lag 0
    ref_lon = ref_lon0 + deg_per_day * ref_lag

    fig, axes = plt.subplots(1, 4, sharey=True, figsize=(20, 7))
    bbox_props = dict(fc="white", ec="k", lw=2)
    lon_formatter = LongitudeFormatter(number_format=".0f")

    for ax, (label, letter, lons, lags, mP, mU) in zip(axes, results):
        ax.contourf(lons, lags, mP, levels=clevels, cmap=cmap, norm=norm, extend="both")
        ax.contour(lons, lags, mU, levels=contlevels_emph, colors="k", alpha=0.8, linewidths=2)
        cl = ax.contour(lons, lags, mU, levels=contlevels, colors="k", alpha=0.3)
        ax.clabel(cl, colors="k", fontsize=12, inline=True, fmt="%1.1f")
        ax.grid(True, alpha=0.2)
        ax.set_xlabel("Longitude", fontsize=22)
        ax.tick_params(labelsize=16)
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.set_xlim([40, 180])
        ax.set_ylim([-25, 25])
        # Maritime-Continent markers
        ax.axvline(x=90, color="black", linestyle="--", alpha=0.7)
        ax.axvline(x=150, color="black", linestyle="--", alpha=0.7)
        # 5 m/s eastward phase-speed reference line
        ax.plot(ref_lon, ref_lag, color="magenta", linestyle="-", lw=2.5, alpha=0.9)
        ax.text(0.98, 0.98, label, transform=ax.transAxes, ha="right", va="top",
                fontsize=22, bbox=bbox_props)

    axes[0].set_ylabel("Lag (days)", fontsize=22)
    axes[-1].text(0.97, 0.30, "5 m s$^{-1}$", transform=axes[-1].transAxes,
                  ha="right", va="center", fontsize=13, color="magenta")

    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.18, top=0.92, wspace=0.09)
    cbar_ax = fig.add_axes([0.30, 0.02, 0.40, 0.035])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, extend="both",
                                   spacing="proportional", ticks=clevels,
                                   boundaries=clevels, orientation="horizontal")
    cbar_ax.set_xlabel("Correlation (precip shaded, U850 contours)", size=18)
    cb.ax.tick_params(labelsize=13)
    cb.locator = ticker.MaxNLocator(nbins=10)
    cb.update_ticks()

    plt.savefig(FIG_OUT, bbox_inches="tight", dpi=300)
    print(f"saved figure to {FIG_OUT}")


if __name__ == "__main__":
    main()
