#!/usr/bin/env python
"""
compute_structural_metrics.py

Objective structural metrics that replace subjective "closer to ERA5" visual
claims with numbers, for the MJO regression-map panels built by make_panels.py.

For each MAP diagnostic (blmc925, u850, diah700, conindex, precip)
the band-pass-filtered field is regressed onto the equatorial Indian-Ocean
precipitation index EXACTLY as in make_panels.py (load_slopes / reg_coef_n are
imported and reused, so the regression is identical). For each model experiment
(CNTRL, MEAN, CNN) the regressed map is compared with the ERA5 reference map
over the Indo-Pacific MJO domain (make_panels EXTENT = 40-180E, 20S-20N) using
cos(lat) area weighting:

  1. Uncentered pattern correlation with the ERA5 map.
  2. Centered (anomaly) pattern correlation with the ERA5 map.
  3. Centered RMS error normalized by the ERA5 map's spatial standard deviation.
  4. Amplitude-weighted "centroid" longitude of the POSITIVE anomaly along the
     equatorial band (10S-10N) -- an east/west extension measure. Reported for
     ERA5 and each experiment.

Additionally, if the Hovmoller (lag-longitude) inputs used by fig_hovmoller.py
are loadable, an objective eastward-propagation metric is computed: the
longitude of the precip correlation ridge as a function of lag, fit to an
eastward phase speed (m/s) and the longitude the ridge reaches by lag +20 days.

Run:  conda run -p /glade/u/apps/opt/conda/envs/npl-2024a python \
        figure_scripts/compute_structural_metrics.py
Output: prints a table; also writes figure_scripts/structural_metrics.txt
"""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# Reuse make_panels' exact loading + regression logic (do NOT re-derive it).
from make_panels import (load_slopes, DIAGNOSTICS, EXPERIMENTS, NCDIR, EXTENT)

import warnings
warnings.filterwarnings("ignore")

# Map diagnostics to score (subset of make_panels DIAGNOSTICS that are true maps)
MAP_DIAGS = ["blmc925", "u850", "diah700", "conindex", "precip"]

# Model experiments (exclude obs entry) and the obs entry, from make_panels.
OBS_EXP = next(e for e in EXPERIMENTS if e["obs"])
MODEL_EXPS = [e for e in EXPERIMENTS if not e["obs"]]

LONMIN, LONMAX, LATMIN, LATMAX = EXTENT      # 40, 180, -20, 20
EQ_LATMIN, EQ_LATMAX = -10.0, 10.0           # equatorial band for centroid


def _domain_mask(lons, lats, latmin, latmax):
    lo = (lons >= LONMIN) & (lons <= LONMAX)
    la = (lats >= latmin) & (lats <= latmax)
    return la, lo


def _weights2d(lats_sub):
    w = np.cos(np.deg2rad(lats_sub))
    w = np.clip(w, 0.0, None)
    return w[:, None]


def slopes_for(exp, cfg):
    """Return (slopes, lons, lats) for one experiment/diagnostic using the
    make_panels loader, or None if the file is missing."""
    suffix = cfg["suffix_obs"] if exp["obs"] else cfg["suffix_model"]
    path = os.path.join(NCDIR, exp["dir"], f"{exp['modout']}.{suffix}")
    if not os.path.exists(path):
        return None
    idx_path = None
    isuf = cfg.get("index_suffix_obs") if exp["obs"] else cfg.get("index_suffix")
    if isuf:
        idx_path = os.path.join(NCDIR, exp["dir"], f"{exp['modout']}.{isuf}")
    slopes, lons, lats, _ = load_slopes(
        path, cfg["field_candidates"], level=cfg.get("level"), index_path=idx_path)
    return slopes, np.asarray(lons), np.asarray(lats)


def pattern_metrics(model, obs, lons, lats):
    """cos(lat)-weighted uncentered corr, centered corr, and centered RMSE
    normalized by the obs (ERA5) weighted spatial std, over the MJO domain."""
    la, lo = _domain_mask(lons, lats, LATMIN, LATMAX)
    m = model[np.ix_(la, lo)].astype(float)
    o = obs[np.ix_(la, lo)].astype(float)
    w = np.broadcast_to(_weights2d(lats[la]), m.shape)

    good = np.isfinite(m) & np.isfinite(o) & (w > 0)
    m, o, w = m[good], o[good], w[good]
    W = w.sum()

    # uncentered pattern correlation
    unc = (w * m * o).sum() / np.sqrt((w * m * m).sum() * (w * o * o).sum())

    # weighted means -> centered fields
    mbar = (w * m).sum() / W
    obar = (w * o).sum() / W
    ma, oa = m - mbar, o - obar
    num = (w * ma * oa).sum()
    den = np.sqrt((w * ma * ma).sum() * (w * oa * oa).sum())
    cen = num / den

    # centered RMSE (Taylor E') normalized by obs weighted spatial std
    erms = np.sqrt((w * (ma - oa) ** 2).sum() / W)
    ostd = np.sqrt((w * oa * oa).sum() / W)
    nrmse = erms / ostd

    # amplitude (weighted spatial std) ratio model/obs -- for the overshoot note
    mstd = np.sqrt((w * ma * ma).sum() / W)
    amp_ratio = mstd / ostd
    return unc, cen, nrmse, amp_ratio


def positive_centroid_lon(field, lons, lats):
    """Amplitude-weighted mean longitude of positive values over 10S-10N,
    40-180E (cos-lat weighted)."""
    la, lo = _domain_mask(lons, lats, EQ_LATMIN, EQ_LATMAX)
    sub = field[np.ix_(la, lo)].astype(float)
    lon_sub = lons[lo]
    lat_sub = lats[la]
    w = np.broadcast_to(_weights2d(lat_sub), sub.shape)
    lon2d = np.broadcast_to(lon_sub[None, :], sub.shape)
    pos = np.isfinite(sub) & (sub > 0)
    wp = (w * sub)[pos]                 # amplitude * area weight
    if wp.sum() <= 0:
        return np.nan
    return (lon2d[pos] * wp).sum() / wp.sum()


# ----------------------------------------------------------------------------
# Build the table
# ----------------------------------------------------------------------------
def build_map_table():
    lines = []
    exp_labels = {e["key"]: e["label"].split(")")[-1].strip() for e in EXPERIMENTS}

    header = (f"{'diag':<9} {'exp':<6} {'unc_corr':>9} {'cen_corr':>9} "
              f"{'norm_RMSE':>9} {'amp_ratio':>9} {'centroid_lon':>12}")
    results = {}   # diag -> dict

    for diag in MAP_DIAGS:
        if diag not in DIAGNOSTICS:
            lines.append(f"# {diag}: not defined in make_panels DIAGNOSTICS -- skipped")
            continue
        cfg = DIAGNOSTICS[diag]
        obs = slopes_for(OBS_EXP, cfg)
        if obs is None:
            lines.append(f"# {diag}: ERA5/obs map MISSING -- skipped")
            continue
        obs_sl, olons, olats = obs
        obs_centroid = positive_centroid_lon(obs_sl, olons, olats)

        block = [header]
        # ERA5 reference row (self-correlation is 1 by construction; show centroid)
        block.append(f"{diag:<9} {'ERA5':<6} {1.0:9.3f} {1.0:9.3f} "
                     f"{0.0:9.3f} {1.0:9.3f} {obs_centroid:12.1f}")
        rows = {"ERA5": dict(centroid=obs_centroid)}
        for exp in MODEL_EXPS:
            got = slopes_for(exp, cfg)
            if got is None:
                block.append(f"{diag:<9} {exp_labels[exp['key']]:<6} "
                             f"{'MISSING':>9}")
                continue
            msl, mlons, mlats = got
            if msl.shape != obs_sl.shape:
                block.append(f"{diag:<9} {exp_labels[exp['key']]:<6} "
                             f"grid-mismatch model{msl.shape} obs{obs_sl.shape}")
                continue
            unc, cen, nrmse, amp = pattern_metrics(msl, obs_sl, mlons, mlats)
            cen_lon = positive_centroid_lon(msl, mlons, mlats)
            block.append(f"{diag:<9} {exp_labels[exp['key']]:<6} "
                         f"{unc:9.3f} {cen:9.3f} {nrmse:9.3f} {amp:9.3f} "
                         f"{cen_lon:12.1f}")
            rows[exp_labels[exp["key"]]] = dict(
                unc=unc, cen=cen, nrmse=nrmse, amp=amp, centroid=cen_lon)
        results[diag] = rows
        lines.append("\n".join(block))
        lines.append("")
    return "\n".join(lines), results


# ----------------------------------------------------------------------------
# Hovmoller eastward-propagation metric (optional; reuses fig_hovmoller inputs)
# ----------------------------------------------------------------------------
def hovmoller_metric():
    out = ["", "=" * 78,
           "EASTWARD-PROPAGATION METRIC (lag-longitude precip correlation ridge)",
           "=" * 78]
    try:
        import fig_hovmoller as fh
    except Exception as e:
        out.append(f"could not import fig_hovmoller: {e}")
        return "\n".join(out)

    # eastward-propagation window (lags 0..+20 days), ridge tracked over 40-180E
    lag_lo, lag_hi = 0, 20
    ridge_lon_range = (60.0, 180.0)   # avoid the IO index anchor edge at 40E

    out.append("def: at each lag, longitude of max precip correlation within "
               f"{ridge_lon_range[0]:.0f}-{ridge_lon_range[1]:.0f}E; slope of "
               f"lon vs lag over lag {lag_lo}..{lag_hi}d -> phase speed; "
               "lon reached at lag +20d = eastward extent.")
    out.append(f"{'exp':<10} {'speed_deg/day':>13} {'speed_m/s':>10} "
               f"{'lon@lag0':>9} {'lon@lag20':>10}")

    for label, letter, pf, pv, uf, uv in fh.EXPS:
        try:
            lons, lags, mP, mU = fh.compute_lagged_corr_lon(pf, pv, uf, uv)
        except Exception as e:
            out.append(f"{label:<10} FAILED: {e}")
            continue
        lons = np.asarray(lons)
        lags = np.asarray(lags)
        lo_mask = (lons >= ridge_lon_range[0]) & (lons <= ridge_lon_range[1])
        lag_sel = (lags >= lag_lo) & (lags <= lag_hi)
        sel_lags = lags[lag_sel]
        ridge_lons = []
        for il, lg in enumerate(lags):
            if not lag_sel[il]:
                continue
            row = mP[il, :].copy()
            row[~lo_mask] = -np.inf
            ridge_lons.append(lons[np.argmax(row)])
        ridge_lons = np.array(ridge_lons, dtype=float)
        # linear fit lon = a + b*lag
        b, a = np.polyfit(sel_lags, ridge_lons, 1)
        speed_degday = b
        speed_ms = b * 111.32 * 1000.0 / 86400.0
        lon0 = a + b * 0
        lon20 = a + b * 20
        clean = label.split(")")[-1].strip()
        out.append(f"{clean:<10} {speed_degday:13.2f} {speed_ms:10.2f} "
                   f"{lon0:9.1f} {lon20:10.1f}")
    return "\n".join(out)


def main():
    table, _ = build_map_table()
    full = ["MJO REGRESSION-MAP STRUCTURAL METRICS vs ERA5 "
            "(Indo-Pacific domain 40-180E, 20S-20N, cos-lat weighted)",
            "unc_corr=uncentered pattern corr; cen_corr=centered pattern corr; "
            "norm_RMSE=centered RMSE / ERA5 spatial std;",
            "amp_ratio=model spatial std / ERA5 spatial std (>1 = overshoot); "
            "centroid_lon=amp-weighted lon of positive anomaly 10S-10N.",
            "", table]
    out_text = "\n".join(full)

    # Hovmoller part (optional)
    hov = ""
    do_hov = ("--no-hov" not in sys.argv)
    if do_hov:
        try:
            hov = hovmoller_metric()
        except Exception as e:
            hov = f"\nEASTWARD-PROPAGATION METRIC: skipped ({e})"
    else:
        hov = "\nEASTWARD-PROPAGATION METRIC: skipped (--no-hov)"
    out_text = out_text + "\n" + hov

    print(out_text)
    dest = os.path.join(HERE, "structural_metrics.txt")
    with open(dest, "w") as fh:
        fh.write(out_text + "\n")
    print(f"\n[compute_structural_metrics] wrote {dest}")


if __name__ == "__main__":
    main()
