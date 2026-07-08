#!/usr/bin/env python
"""
check_thetae_bolton.py

Test the manuscript's claim (Appendix A4) that the first-order equivalent potential
temperature used for the moist-instability index does not materially change the regressed
patterns relative to the exact pseudo-adiabatic formulation of Bolton (1980).

The paper's index is I = theta_e(850 hPa) - theta_e(400 hPa), regressed on the standardized
Indian Ocean precipitation index. We compute I with two formulations of theta_e from the same
T, q, p, regress both, and compare the resulting maps.

  First order (Eq. A4 of the manuscript):
      theta_e = (T + Lv/cp * r) * (p0/p)**(R/cp),        r = q/(1-q)

  Bolton (1980), Eqs. (15) and (38) -- exact pseudo-adiabatic:
      e   = p*q / (0.622 + 0.378*q)                       [vapour pressure, hPa]
      T_L = 2840 / (3.5*ln T - ln e - 4.805) + 55         [LCL temperature, K]
      theta_e = T * (1000/p)**(0.2854*(1 - 0.28*r))
                  * exp( (3376/T_L - 2.54) * r * (1 + 0.81*r) )

Filtering, seasons, index, and regression follow the paper exactly (see fig_increment.py).

Input : an ncks subset of the CNN plev.mjo timeseries holding T, Q at 400 and 850 hPa + PRECT:
    ncks -O -4 -L1 -v T,Q,PRECT -d lat,-25.,25. -d lon,40.,220. -d level,6,10,4 \
         <exp>.cam.h1.NC4_Classic.plev.mjo.*.nc  tq_subset.nc
Output: Paper_WCD/figures/thetae_bolton_check.txt
"""

import os
import sys
import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig_increment import (lanczos_bandpass_weights, bandpass, season_indices, stack_seasons,
                           regress, area_mean, NWGT, FCA, FCB, LAT_IO)

ROOT = os.path.dirname(HERE)
FIGDIR = os.path.join(ROOT, "figures")

DEFAULT = os.environ.get("TQ_SUBSET",
    "/glade/derecho/scratch/wchapman/tmp/claude-32755/"
    "-glade-work-wchapman-DA-ML-CESML-AI-Paper-Figures-Wang/"
    "593824dd-faa9-4cc0-8c37-90e1784b2e57/scratchpad/tq_subset.nc")

LV, CP, RD, P0 = 2.51e6, 1004.0, 287.0, 1000.0   # manuscript's A4 constants
IDX_BOX = (80.0, 100.0)
DOM_LAT, DOM_LON = (-20, 20), (40, 180)


def theta_e_first_order(T, q, p):
    r = q / (1.0 - q)
    return (T + (LV / CP) * r) * (P0 / p) ** (RD / CP)


def theta_e_bolton(T, q, p):
    """Bolton (1980) Eqs. 15 and 38. p in hPa, q in kg/kg."""
    r = q / (1.0 - q)
    e = p * q / (0.622 + 0.378 * q)                      # hPa
    e = np.maximum(e, 1e-10)
    T_L = 2840.0 / (3.5 * np.log(T) - np.log(e) - 4.805) + 55.0
    return (T * (1000.0 / p) ** (0.2854 * (1.0 - 0.28 * r))
            * np.exp((3376.0 / T_L - 2.54) * r * (1.0 + 0.81 * r)))


def metrics(a, b, lat, lon):
    """Centered pattern correlation and normalized centered RMSE of b against a."""
    jm = (lat >= DOM_LAT[0]) & (lat <= DOM_LAT[1])
    im = (lon >= DOM_LON[0]) & (lon <= DOM_LON[1])
    w = np.cos(np.deg2rad(lat))[jm][:, None]
    A, B = a[np.ix_(jm, im)], b[np.ix_(jm, im)]
    W = np.broadcast_to(w, A.shape)
    am = np.sum(A * W) / np.sum(W)
    bm = np.sum(B * W) / np.sum(W)
    Ac, Bc = A - am, B - bm
    sa = np.sqrt(np.sum(W * Ac ** 2) / np.sum(W))
    sb = np.sqrt(np.sum(W * Bc ** 2) / np.sum(W))
    r = np.sum(W * Ac * Bc) / np.sum(W) / (sa * sb)
    ncrmse = np.sqrt(np.sum(W * (Ac - Bc) ** 2) / np.sum(W)) / sa
    return r, ncrmse, sa, sb


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    ds = xr.open_dataset(path, decode_times=True)
    lat, lon = ds["lat"].values, ds["lon"].values
    lev = ds["level"].values
    i850 = int(np.argmin(np.abs(lev - 850.0)))
    i400 = int(np.argmin(np.abs(lev - 400.0)))
    print(f"[thetae] levels {lev}; using 850 at idx {i850}, 400 at idx {i400}")

    T = ds["T"].values
    q = ds["Q"].values
    p = lev[None, :, None, None] * np.ones_like(T)

    out = {}
    for name, fn in (("first_order", theta_e_first_order), ("bolton", theta_e_bolton)):
        th = fn(T, q, p)
        out[name] = th[:, i850] - th[:, i400]          # I(t, lat, lon)
        print(f"[thetae] {name:11s}: mean I = {np.nanmean(out[name]):8.3f} K, "
              f"sd = {np.nanstd(out[name]):6.3f} K")

    w = lanczos_bandpass_weights(NWGT, FCA, FCB)
    starts = season_indices(ds["time"].values)
    pr = bandpass(ds["PRECT"].values * 86400.0 * 1000.0, w, 0)
    pr_da = xr.DataArray(pr, coords={"time": ds["time"], "lat": ds["lat"], "lon": ds["lon"]},
                         dims=("time", "lat", "lon"))
    idx = stack_seasons(area_mean(pr_da, LAT_IO[0], LAT_IO[1], *IDX_BOX), starts, 0)
    idx = idx / np.nanstd(idx, ddof=1)

    reg = {}
    for name in out:
        f = stack_seasons(bandpass(out[name], w, 0), starts, 0)
        reg[name] = regress(idx, f)
        print(f"[thetae] regressed {name:11s}: max |I'| = {np.nanmax(np.abs(reg[name])):.4f} K per sd")

    r, ncrmse, s_fo, s_bo = metrics(reg["first_order"], reg["bolton"], lat, lon)

    L = []
    P = L.append
    P("Sensitivity of the moist-instability index to the theta_e formulation")
    P("=" * 76)
    P(f"source            : {os.path.basename(path)} (CNN experiment)")
    P(f"winters           : {len(starts)}, 180-day seasons from 1 Nov")
    P(f"domain for metrics: {DOM_LON[0]}-{DOM_LON[1]}E, {DOM_LAT[0]}-{DOM_LAT[1]}N")
    P("")
    P("I = theta_e(850 hPa) - theta_e(400 hPa), regressed on the standardized Indian Ocean")
    P("precipitation index (K per index standard deviation).")
    P("")
    P("  first-order (manuscript Eq. A4) vs exact pseudo-adiabatic (Bolton 1980, Eqs. 15, 38)")
    P("")
    P(f"  centered pattern correlation of the two regressed maps : {r:.4f}")
    P(f"  normalized centered RMSE (Bolton vs first-order)       : {ncrmse:.4f}")
    P(f"  spatial sd, first-order                                : {s_fo:.4f} K per sd")
    P(f"  spatial sd, Bolton                                     : {s_bo:.4f} K per sd")
    P(f"  amplitude ratio (Bolton / first-order)                 : {s_bo/s_fo:.4f}")
    P("")
    if r > 0.99 and abs(s_bo / s_fo - 1.0) < 0.15:
        P("VERDICT: the two formulations give the same regressed pattern to within a few percent;")
        P("the manuscript's claim that the first-order form does not materially change the")
        P("regressed patterns is supported.")
    else:
        P("VERDICT: the formulations differ materially. The manuscript's claim is NOT supported")
        P("and must be softened or the exact formulation adopted.")
    txt = "\n".join(L)
    print("\n" + txt)
    os.makedirs(FIGDIR, exist_ok=True)
    with open(os.path.join(FIGDIR, "thetae_bolton_check.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
