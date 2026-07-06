#!/usr/bin/env python
"""
compute_kr_boxes.py

Robust box-averaged Kelvin-like/Rossby-like amplitude ratio (r_KR) from the
regressed 850 hPa zonal wind, replacing the fragile single-gridpoint extrema
metric max(-U850)/max(U850).

Self-contained: opens the same nc products that make_panels.py "u850" uses,
reproduces reg_coef_n (per-gridpoint linregress slope of field on the
standardized Indian-Ocean precip index), then computes cos-lat-weighted box
averages and a season-block bootstrap CI. Does NOT import make_panels.

File/var conventions copied from make_panels.py:
  model: <modout>.mjo_U850_model_all.850mb.80.100.nc  field 'u_timem'  index 'piom'
  obs:   <modout>.mjo_U850_obs_all.850mb.nc            field 'u_time'   index 'pio'
"""
import os
import numpy as np
import xarray as xr
from scipy.stats import linregress
import warnings
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NCDIR = os.path.join(ROOT, "nc")

DATASETS = [
    dict(key="ERA5",  dir="CNN",   modout="f.e.FTORCHmjo_fullCNN_DT2", obs=True),
    dict(key="CNTRL", dir="CNTRL", modout="f.e.FTORCHmjo_CNTRLmjo_DT2", obs=False),
    dict(key="MEAN",  dir="MEAN",  modout="f.e.FTORCHmjo_MEANmjo_DT2",  obs=False),
    dict(key="CNN",   dir="CNN",   modout="f.e.FTORCHmjo_fullCNN_DT2",  obs=False),
]

SUFFIX_OBS   = "mjo_U850_obs_all.850mb.nc"
SUFFIX_MODEL = "mjo_U850_model_all.850mb.80.100.nc"

# analysis domain for the OLD extrema metric
DOM_LAT = (-10.0, 10.0)
DOM_LON = (40.0, 180.0)

# box definitions: (name, (lonW,lonE))  all on 10S-10N equatorial band
PRIMARY = dict(west=(60.0, 90.0),  east=(110.0, 150.0))
ALT1    = dict(west=(50.0, 90.0),  east=(120.0, 160.0))
ALT2    = dict(west=(60.0, 100.0), east=(120.0, 160.0))
BOXSETS = [("primary", PRIMARY), ("alt1", ALT1), ("alt2", ALT2)]

NBOOT = 1000
WINLEN = 180  # days per winter block


def load(ds_def):
    """Return (field[t,lat,lon], idx[t], lat, lon) subset to equatorial domain."""
    suffix = SUFFIX_OBS if ds_def["obs"] else SUFFIX_MODEL
    path = os.path.join(NCDIR, ds_def["dir"], f"{ds_def['modout']}.{suffix}")
    ds = xr.open_dataset(path)
    idx_name = "piom" if "piom" in ds else "pio"
    fld_name = "u_timem" if "u_timem" in ds else "u_time"
    idx = ds[idx_name].values.squeeze().astype(float)
    lat = ds["lat"].values
    lon = ds["lon"].values
    jlat = np.where((lat >= DOM_LAT[0]) & (lat <= DOM_LAT[1]))[0]
    ilon = np.where((lon >= DOM_LON[0]) & (lon <= DOM_LON[1]))[0]
    field = ds[fld_name].values[:, jlat, :][:, :, ilon].astype(float)
    ds.close()
    return field, idx, lat[jlat], lon[ilon], path


def reg_slope_vec(x, y):
    """Vectorized per-gridpoint OLS slope of y[t,lat,lon] on x[t].
    Mathematically identical to scipy linregress(x, y).slope."""
    x = np.asarray(x, float)
    xd = x - x.mean()
    denom = (xd * xd).sum()
    yd = y - y.mean(axis=0, keepdims=True)
    return (xd[:, None, None] * yd).sum(axis=0) / denom


def box_mean(slopes, lat, lon, lonrange):
    """cos-lat-weighted area mean of slope map over full lat band, lon in range."""
    ii = np.where((lon >= lonrange[0]) & (lon <= lonrange[1]))[0]
    w = np.cos(np.deg2rad(lat))
    sub = slopes[:, ii]                 # (lat, nlon_box)
    wcol = w[:, None] * np.ones((1, sub.shape[1]))
    return np.nansum(sub * wcol) / np.nansum(wcol)


def kr_box(slopes, lat, lon, boxset):
    uw = box_mean(slopes, lat, lon, boxset["west"])
    ue = box_mean(slopes, lat, lon, boxset["east"])
    return (-ue) / uw, uw, ue


def kr_extrema(slopes):
    return np.max(-slopes) / np.max(slopes)


def main():
    print("=" * 78)
    print("Box-averaged Kelvin/Rossby amplitude ratio  r_KR_box = -mean(U_east)/mean(U_west)")
    print("cos-lat weighted, 10S-10N equatorial band")
    print("=" * 78)

    results = {}
    for d in DATASETS:
        field, idx, lat, lon, path = load(d)
        nwin = field.shape[0] // WINLEN
        clean = (field.shape[0] % WINLEN == 0)

        # full-precision slope map via scipy linregress per gridpoint (reg_coef_n)
        ny, nx = field.shape[1], field.shape[2]
        slopes_lr = np.zeros((ny, nx))
        for j in range(ny):
            for i in range(nx):
                slopes_lr[j, i] = linregress(idx, field[:, j, i]).slope
        # cross-check the vectorized version we use in the bootstrap
        slopes_vec = reg_slope_vec(idx, field)
        maxdiff = np.max(np.abs(slopes_lr - slopes_vec))

        extr = kr_extrema(slopes_lr)
        boxvals = {}
        for bname, bset in BOXSETS:
            r, uw, ue = kr_box(slopes_lr, lat, lon, bset)
            boxvals[bname] = (r, uw, ue)

        # season-block bootstrap on the PRIMARY boxes
        cis = {}
        if clean:
            f4 = field.reshape(nwin, WINLEN, ny, nx)
            x4 = idx.reshape(nwin, WINLEN)
            rng = np.random.default_rng(12345)
            draws = {n: np.empty(NBOOT) for n, _ in BOXSETS}
            for b in range(NBOOT):
                sel = rng.integers(0, nwin, size=nwin)
                fb = f4[sel].reshape(nwin * WINLEN, ny, nx)
                xb = x4[sel].reshape(nwin * WINLEN)
                sl = reg_slope_vec(xb, fb)          # slope map once per draw
                for bname, bset in BOXSETS:
                    draws[bname][b], _, _ = kr_box(sl, lat, lon, bset)
            for bname, _ in BOXSETS:
                cis[bname] = (np.percentile(draws[bname], 5),
                              np.percentile(draws[bname], 95))
            boot_mode = f"season-block ({nwin} winters resampled)"
        else:
            boot_mode = "SKIPPED (time not divisible by 180)"

        results[d["key"]] = dict(extr=extr, box=boxvals, ci=cis,
                                 nwin=nwin, clean=clean, boot=boot_mode,
                                 maxdiff=maxdiff)

    # ---- report -----------------------------------------------------------
    print(f"\nvectorized-vs-linregress slope max abs diff (should be ~0):")
    for k in results:
        print(f"   {k:6s}: {results[k]['maxdiff']:.3e}   bootstrap: {results[k]['boot']}")

    print("\n--- OLD extrema metric  max(-U850)/max(U850)  (cross-check) ---")
    print("   expected: ERA5 1.06  CNTRL 0.84  MEAN 0.70  CNN 0.92")
    for k in ["ERA5", "CNTRL", "MEAN", "CNN"]:
        print(f"   {k:6s}: {results[k]['extr']:.3f}")

    print("\n--- PRIMARY boxes  West 60-90E / East 110-150E ---")
    print(f"   {'ds':6s} {'r_KR_box':>9s} {'90% CI':>18s}   {'U_west':>8s} {'U_east':>8s}")
    for k in ["ERA5", "CNTRL", "MEAN", "CNN"]:
        r, uw, ue = results[k]["box"]["primary"]
        ci = results[k]["ci"].get("primary", (np.nan, np.nan))
        print(f"   {k:6s} {r:9.3f}  [{ci[0]:6.3f}, {ci[1]:6.3f}]   {uw:8.3f} {ue:8.3f}")

    print("\n--- SENSITIVITY across box choices ---")
    hdr = "   {:6s} " + " ".join(f"{n:>9s}" for n, _ in BOXSETS) + "     range"
    print(hdr.format("ds"))
    for k in ["ERA5", "CNTRL", "MEAN", "CNN"]:
        vals = [results[k]["box"][n][0] for n, _ in BOXSETS]
        line = f"   {k:6s} " + " ".join(f"{v:9.3f}" for v in vals)
        line += f"     {min(vals):.3f}-{max(vals):.3f}"
        print(line)

    print("\n--- box definitions ---")
    for n, b in BOXSETS:
        print(f"   {n:8s} West {b['west']} East {b['east']}")

    print("\n--- CI for all box sets ---")
    for k in ["ERA5", "CNTRL", "MEAN", "CNN"]:
        parts = []
        for n, _ in BOXSETS:
            ci = results[k]["ci"].get(n, (np.nan, np.nan))
            parts.append(f"{n}:[{ci[0]:.3f},{ci[1]:.3f}]")
        print(f"   {k:6s} " + "  ".join(parts))


if __name__ == "__main__":
    main()
