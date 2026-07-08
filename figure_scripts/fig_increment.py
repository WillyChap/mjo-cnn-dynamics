#!/usr/bin/env python
"""
fig_increment.py

Regress the CNN wind-increment tendency ITSELF onto the MJO precipitation index.

Every other figure in this paper diagnoses the model's *response* to the correction.
This one diagnoses the *correction*: does the learned increment that CAM6 actually applies
(cb24cnn_U, cb24cnn_V; "U/V DAMLining Tendency", m s^-1 s^-1) have a coherent structure
phase-locked to MJO convection, or does it merely add unstructured intraseasonal variance?

Following the Panel_Figures/Model_Increments.ipynb design, the increment is vertically
integrated over the lower troposphere (700-1000 hPa, mass-weighted by interface thickness)
before regression. The layer integral is essential: single-level increments are dominated
by grid-scale structure, and their horizontal divergence is not planetary-scale coherent.
Metrics are insensitive to dropping 1000 hPa (conv_E = +0.52 for 700-1000, +0.57 for
700-925, +0.68 for 850 alone; 11/11 winters positive in each case).

Conventions match the NCL diagnostics (mjo_cb24cnnU/V.ncl, mjo_blmc925.ncl) so the result is
directly comparable to Figs. 3 and 9:
  * 20-70 day Lanczos band-pass, 141 weights, sigma = 1  (Duchon 1979)
  * filter the FULL daily record, THEN extract 180-day seasons beginning 1 November
  * standardize the index over the concatenated seasons; regress each field on it at zero lag
  * increments converted to m s^-1 day^-1 (x 86400)

Two reference boxes, both 10S-10N:
  IO = 80-100E   (the canonical index used by every other figure in the paper)
  MC = 120-140E  (moves the reference east: does the increment structure track convection?)

Inference note: there are only 11 independent winters. A percentile bootstrap over 11
season-blocks runs narrow, and the 141-weight filter spreads information +/-70 days so
sub-season blocks are not independent. Significance is therefore reported by
compute_increment_metrics.py as an exact per-winter sign test plus a leave-one-winter-out
jackknife, not as a block bootstrap.

Input: a lat/lon/level subset of the CNN experiment's plev.mjo timeseries:
    ncks -O -4 -L1 -v cb24cnn_U,cb24cnn_V,U,V,PRECT -d lat,-25.,25. -d lon,40.,220. -d level,9,12 \
         <exp>.cam.h1.NC4_Classic.plev.mjo.1979010100000-1990123100000.nc  inc_subset.nc
Pass its path as argv[1] or via the INC_SUBSET environment variable.

Output: Paper_WCD/figures/fig_increment.png  (+ fig_increment_cache.npz)
"""

import os
import sys
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.ndimage import convolve1d
import warnings
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGDIR = os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

DEFAULT_SUBSET = os.environ.get(
    "INC_SUBSET",
    "/glade/derecho/scratch/wchapman/tmp/claude-32755/"
    "-glade-work-wchapman-DA-ML-CESML-AI-Paper-Figures-Wang/"
    "593824dd-faa9-4cc0-8c37-90e1784b2e57/scratchpad/inc_subset.nc")

# --- analysis constants (mirror the NCL scripts) ------------------------------
NWGT   = 141
FCA    = 1.0 / 70.0
FCB    = 1.0 / 20.0
NDAY   = 180
MON0, DAY0 = 11, 1
LAT_IO = (-10.0, 10.0)
BOXES  = {"IO": (80.0, 100.0), "MC": (120.0, 140.0)}
LEVS   = np.array([700.0, 850.0, 925.0, 1000.0])
DP     = np.array([75.0, 112.5, 75.0, 37.5])      # interface thicknesses -> mass weights
EXTENT = [40, 180, -20, 20]
AEARTH = 6.371e6

# Colour choice. Row 2 is the 850-ish zonal wind, the same field as Fig. 9, so it reuses
# Fig. 9's RdBu_r: the same variable should read the same way across the paper. Row 1 is a
# new quantity (a tendency, not a state) and gets PRGn, a ColorBrewer diverging,
# colorblind-safe map distinct from the BrBG teal-brown of Figs. 3-4.
CMAP_CONV = "PRGn"    # convergence tendency: green = the correction generates convergence
CMAP_U    = "RdBu_r"  # MJO zonal wind anomaly (matches Fig. 9)

# The raw increment carries substantial grid-scale variance (only ~45% of the equatorial
# u-increment variance is at zonal wavenumber <= 6), so it must be smoothed before it is
# differentiated. The scalar box metrics in compute_increment_metrics.py use the Gauss
# theorem instead and need no smoothing; they confirm the sign of the smoothed field.
SMOOTH_PRE  = 12      # smth9 passes applied to the increment before taking its divergence
SMOOTH_POST = 3       # smth9 passes applied to the resulting divergence, for display


def lanczos_bandpass_weights(nwt, fca, fcb):
    """Duchon (1979) Lanczos band-pass weights; matches NCL filwgts_lanczos(ihp=2, sigma=1)."""
    if nwt % 2 == 0:
        raise ValueError("nwt must be odd")
    nn = (nwt - 1) // 2
    w = np.zeros(nwt)
    w[nn] = 2.0 * (fcb - fca)
    k = np.arange(1, nn + 1)
    sigma = np.sin(np.pi * k / nn) / (np.pi * k / nn)
    lobe = (np.sin(2 * np.pi * fcb * k) / (np.pi * k)
            - np.sin(2 * np.pi * fca * k) / (np.pi * k)) * sigma
    w[nn + 1:] = lobe
    w[:nn] = lobe[::-1]
    return w


def bandpass(arr, w, axis=0):
    """Apply the weights along `axis`; the nn edge points on each end become NaN,
    exactly as NCL's wgt_runave_leftdim leaves them missing."""
    nn = (len(w) - 1) // 2
    out = convolve1d(arr, w, axis=axis, mode="constant", cval=0.0)
    sl = [slice(None)] * arr.ndim
    sl[axis] = slice(0, nn);      out[tuple(sl)] = np.nan
    sl[axis] = slice(-nn, None);  out[tuple(sl)] = np.nan
    return out


def season_indices(time):
    """Start indices of each 180-day season beginning 1 Nov, fully inside the record."""
    months = np.array([t.month for t in time])
    days = np.array([t.day for t in time])
    starts = np.where((months == MON0) & (days == DAY0))[0]
    return np.array([s for s in starts if s + NDAY <= len(time)])


def stack_seasons(arr, starts, axis=0):
    return np.concatenate([np.take(arr, np.arange(s, s + NDAY), axis=axis) for s in starts], axis=axis)


def regress(index, field):
    """Zero-lag least-squares slope of field (t,...) on index (t,). NaN-safe."""
    x = index - np.nanmean(index)
    xb = x.reshape((-1,) + (1,) * (field.ndim - 1))
    return np.nansum(xb * (field - np.nanmean(field, axis=0)), axis=0) / np.nansum(x * x)


def layer_mean(ds, var):
    """Mass-weighted 700-1000 hPa mean (m s^-1 s^-1 for increments, m s^-1 for winds)."""
    wts = DP / DP.sum()
    return np.einsum("tzyx,z->tyx", ds[var].sel(level=LEVS).values, wts)


def sph_divergence(u, v, lat, lon):
    phi = np.deg2rad(lat)[:, None]
    dlam = np.deg2rad(np.gradient(lon))[None, :]
    dphi = np.deg2rad(np.gradient(lat))[:, None]
    cos = np.cos(phi)
    return (np.gradient(u, axis=1) / dlam + np.gradient(v * cos, axis=0) / dphi) / (AEARTH * cos)


def smth9(a, p=0.50, q=0.25):
    """NCL smth9: 9-point smoother, centre weight p, neighbour weight q."""
    c = np.pad(a, 1, mode="edge")
    cross = c[:-2, 1:-1] + c[2:, 1:-1] + c[1:-1, :-2] + c[1:-1, 2:] - 4 * a
    diag = c[:-2, :-2] + c[:-2, 2:] + c[2:, :-2] + c[2:, 2:] - 4 * a
    return a + (p / 4.0) * cross + (q / 4.0) * diag


def smoothn(a, n):
    for _ in range(n):
        a = smth9(a)
    return a


def area_mean(da, lat0, lat1, lon0, lon1):
    sub = da.sel(lat=slice(lat0, lat1), lon=slice(lon0, lon1))
    w = np.cos(np.deg2rad(sub.lat))
    return sub.weighted(w).mean(dim=("lat", "lon")).values


def compute(path):
    print(f"[fig_increment] reading {path}")
    ds = xr.open_dataset(path, decode_times=True)
    lat, lon = ds["lat"].values, ds["lon"].values
    w = lanczos_bandpass_weights(NWGT, FCA, FCB)
    starts = season_indices(ds["time"].values)
    print(f"[fig_increment] {len(starts)} winters, {NDAY}-day seasons from "
          f"{ds['time'].values[starts[0]]}")

    iu = bandpass(layer_mean(ds, "cb24cnn_U") * 86400.0, w, 0)   # m s^-1 day^-1
    iv = bandpass(layer_mean(ds, "cb24cnn_V") * 86400.0, w, 0)
    bu = bandpass(layer_mean(ds, "U"), w, 0)                     # m s^-1
    pr = bandpass(ds["PRECT"].values * 86400.0 * 1000.0, w, 0)   # mm day^-1

    pr_da = xr.DataArray(pr, coords={"time": ds["time"], "lat": ds["lat"], "lon": ds["lon"]},
                         dims=("time", "lat", "lon"))
    ius, ivs, bus, prs = (stack_seasons(a, starts, 0) for a in (iu, iv, bu, pr))

    out = {}
    for name, (l0, l1) in BOXES.items():
        idx = stack_seasons(area_mean(pr_da, LAT_IO[0], LAT_IO[1], l0, l1), starts, 0)
        idx = idx / np.nanstd(idx, ddof=1)
        ru, rv, rbu, rp = (regress(idx, a) for a in (ius, ivs, bus, prs))
        su, sv = smoothn(ru, SMOOTH_PRE), smoothn(rv, SMOOTH_PRE)
        conv = -sph_divergence(su, sv, lat, lon) * 1e6           # 1e-6 s^-1 day^-1
        # plot the smoothed increment vectors too, so the arrows match the field they imply
        out[name] = dict(u=su, v=sv, conv=smoothn(conv, SMOOTH_POST),
                         bu=smth9(rbu), precip=smth9(rp))
        print(f"[fig_increment] {name}: |dV|max={np.nanmax(np.hypot(ru, rv)):.3f} m/s/day, "
              f"conv range [{np.nanmin(out[name]['conv']):+.2f}, {np.nanmax(out[name]['conv']):+.2f}], "
              f"U' max={np.nanmax(np.abs(rbu)):.2f} m/s, precip max={np.nanmax(rp):.2f} mm/day")

    np.savez(os.path.join(FIGDIR, "fig_increment_cache.npz"), lat=lat, lon=lon,
             **{f"{n}_{k}": out[n][k] for n in BOXES for k in ("u", "v", "conv", "bu", "precip")})
    return out, lat, lon


def plot(res, lat, lon):
    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(13.4, 5.8))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.07, hspace=0.30)

    rows = [("conv", CMAP_CONV, np.arange(-0.9, 0.95, 0.15),
             "$-\\nabla\\!\\cdot\\,\\overline{\\delta\\mathbf{V}}$\n[$10^{-6}$ s$^{-1}$ day$^{-1}$]"),
            ("bu", CMAP_U, np.arange(-1.2, 1.25, 0.15),
             "$\\overline{U'}$\n[m s$^{-1}$]")]
    # Column identity (which index box) is given in the caption, not on the panels.
    cols = [("IO", ""), ("MC", "")]
    tags = "abcd"
    mappables = {}
    skip = 5

    for r, (key, cmap, clevs, cbl) in enumerate(rows):
        for c, (box, boxlab) in enumerate(cols):
            ax = fig.add_subplot(gs[r, c], projection=proj)
            d = res[box]
            ax.coastlines("50m", linewidth=0.5)
            mappables[key] = ax.contourf(lon, lat, d[key], levels=clevs, cmap=cmap,
                                         extend="both", transform=proj)
            pr = d["precip"]
            ax.contour(lon, lat, pr, levels=[-1.5, -0.75], colors="0.35",
                       linewidths=0.8, linestyles="dashed", transform=proj)
            ax.contour(lon, lat, pr, levels=[0.75, 1.5], colors="k",
                       linewidths=1.3, transform=proj)
            q = ax.quiver(lon[::skip], lat[::skip], d["u"][::skip, ::skip], d["v"][::skip, ::skip],
                          pivot="middle", color="0.12", scale=8.5, width=0.0038,
                          headwidth=4.0, headlength=4.5, alpha=0.8, transform=proj)
            if r == 0 and c == 1:
                ax.quiverkey(q, 0.78, 1.11, 0.3, r"0.3 m s$^{-1}$ day$^{-1}$",
                             labelpos="E", coordinates="axes", fontproperties={"size": 11})
            l0, l1 = BOXES[box]
            ax.add_patch(Rectangle((l0, LAT_IO[0]), l1 - l0, LAT_IO[1] - LAT_IO[0],
                                   lw=1.6, edgecolor="crimson", facecolor="none",
                                   transform=proj, zorder=6))
            ax.set_xticks(np.arange(60, 181, 40), crs=proj)
            ax.set_yticks(np.arange(-20, 21, 10), crs=proj)
            ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".0f"))
            ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".0f"))
            ax.tick_params(labelsize=13)
            if c == 1:
                ax.set_yticklabels([])
            ax.set_extent(EXTENT, crs=proj)
            lbl = f"{tags[r * 2 + c]}){(' ' + boxlab) if boxlab else ''}"
            ax.text(0.985, 0.90, lbl, transform=ax.transAxes,
                    ha="right", va="top", fontsize=13,
                    bbox=dict(fc="white", ec="k", lw=0.8))

    for r, (key, cmap, clevs, cbl) in enumerate(rows):
        cax = fig.add_axes([0.925, 0.545 - 0.435 * r, 0.013, 0.315])
        cb = fig.colorbar(mappables[key], cax=cax, orientation="vertical")
        cb.set_label(cbl, fontsize=10.5, labelpad=8)
        cb.ax.tick_params(labelsize=9.5)

    out = os.path.join(FIGDIR, "fig_increment.png")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[fig_increment] wrote {out}")


if __name__ == "__main__":
    res, lat, lon = compute(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SUBSET)
    plot(res, lat, lon)
