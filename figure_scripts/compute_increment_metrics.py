#!/usr/bin/env python
"""
compute_increment_metrics.py

Quantify what the CNN wind increment does to the MJO, with inference appropriate to a
record of only 11 independent winters.

WHY NOT A BLOCK BOOTSTRAP.  The regression composites use a 141-weight, 20-70 day Lanczos
band-pass, which spreads information +/-70 days. Sub-season blocks are therefore not
independent, and the only defensible block is a whole winter. With 11 blocks a percentile
bootstrap has an effective sample size of 11 and runs narrow. We instead use inference that
is exact at n = 11:
  * an exact two-sided binomial SIGN TEST over the 11 winters (each winter regressed on its
    own index, so the winters are independent realizations of the same statistic), and
  * a leave-one-winter-out JACKKNIFE to show that no single winter carries the result.
The block bootstrap is still reported for comparison, flagged as anticonservative.

METRICS (all from the 700-1000 hPa mass-weighted layer mean, regressed on the standardized
Indian Ocean precipitation index, per index standard deviation):

  conv_E   box-mean convergence of the increment, 100-130E / 10S-10N  [1e-6 s^-1 day^-1]
  conv_W   the same, 55-80E / 10S-10N (west of convection)
  work     <U' dU' + V' dV'> over 40-180E, 20S-20N  [m2 s-2 day-1 per (index sd)^2]
           positive => the increment does positive work on the MJO circulation
  prop     <dU'/dx dU' + dV'/dx dV'>; negative => the tendency displaces the wind
           pattern eastward. Reported as suggestive only (see the sign test).

Box means use the Gauss divergence theorem on the box boundary rather than an average of
pointwise derivatives, so interior grid-scale noise in the increment cancels exactly.

Usage:  python compute_increment_metrics.py [inc_subset.nc]
Writes: Paper_WCD/figures/increment_metrics.txt
"""

import os
import sys
import numpy as np
import xarray as xr
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig_increment import (lanczos_bandpass_weights, bandpass, season_indices, layer_mean,
                           area_mean, NWGT, FCA, FCB, NDAY, LAT_IO, AEARTH, DEFAULT_SUBSET)

ROOT = os.path.dirname(HERE)
FIGDIR = os.path.join(ROOT, "figures")

IDX_BOX = (80.0, 100.0)
BOX_E = (-10, 10, 100, 130)
BOX_W = (-10, 10, 55, 80)
DOM = dict(lat=(-20, 20), lon=(40, 180))
NBOOT = 4000


def sign_test(vals, expect_positive):
    """Exact two-sided binomial test that the per-winter values share the expected sign."""
    n = len(vals)
    k = sum(v > 0 for v in vals) if expect_positive else sum(v < 0 for v in vals)
    p_one = sum(comb(n, j) for j in range(k, n + 1)) / 2 ** n
    return k, min(1.0, 2 * p_one)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SUBSET
    ds = xr.open_dataset(path, decode_times=True)
    lat, lon = ds["lat"].values, ds["lon"].values
    w = lanczos_bandpass_weights(NWGT, FCA, FCB)
    starts = season_indices(ds["time"].values)
    NS = len(starts)
    years = [ds["time"].values[s].year for s in starts]

    iu = bandpass(layer_mean(ds, "cb24cnn_U") * 86400.0, w, 0)
    iv = bandpass(layer_mean(ds, "cb24cnn_V") * 86400.0, w, 0)
    bu = bandpass(layer_mean(ds, "U"), w, 0)
    bv = bandpass(layer_mean(ds, "V"), w, 0)
    pr = bandpass(ds["PRECT"].values * 86400.0 * 1000.0, w, 0)
    pr_da = xr.DataArray(pr, coords={"time": ds["time"], "lat": ds["lat"], "lon": ds["lon"]},
                         dims=("time", "lat", "lon"))
    idx_full = area_mean(pr_da, LAT_IO[0], LAT_IO[1], *IDX_BOX)

    # per-season partial sums: slope over any subset of winters in O(NS)
    sl = [np.arange(s, s + NDAY) for s in starts]
    X = np.stack([idx_full[k] for k in sl])
    FLD = {n: np.stack([a[k] for k in sl]) for n, a in
           (("iu", iu), ("iv", iv), ("bu", bu), ("bv", bv))}
    C_s, D_s = X.sum(1), (X ** 2).sum(1)
    A_s = {n: np.einsum("st,styx->syx", X, f) for n, f in FLD.items()}
    B_s = {n: f.sum(1) for n, f in FLD.items()}

    def slopes(pick):
        """Regression slope on the index restandardized within `pick` (a set of winters)."""
        N = NDAY * len(pick)
        C, D = C_s[pick].sum(), D_s[pick].sum()
        den = D - C * C / N
        return {n: (A_s[n][pick].sum(0) - (C / N) * B_s[n][pick].sum(0)) / np.sqrt(den * (N - 1))
                for n in FLD}

    def box_conv(u, v, lat0, lat1, lon0, lon1):
        """Area-mean CONVERGENCE from the boundary flux (Gauss). 1e-6 s^-1 day^-1."""
        j = np.where((lat >= lat0) & (lat <= lat1))[0]
        i = np.where((lon >= lon0) & (lon <= lon1))[0]
        phi, lam = np.deg2rad(lat[j]), np.deg2rad(lon[i])
        U, V = u[np.ix_(j, i)], v[np.ix_(j, i)]
        flux = ((np.trapz(U[:, -1], phi) - np.trapz(U[:, 0], phi))
                + (np.trapz(V[-1, :] * np.cos(phi[-1]), lam)
                   - np.trapz(V[0, :] * np.cos(phi[0]), lam)))
        area = (lam[-1] - lam[0]) * (np.sin(phi[-1]) - np.sin(phi[0]))
        return -flux / (area * AEARTH) * 1e6

    cw = np.cos(np.deg2rad(lat))[:, None]
    dm = (lat >= DOM["lat"][0]) & (lat <= DOM["lat"][1])
    lm_ = (lon >= DOM["lon"][0]) & (lon <= DOM["lon"][1])

    def awmean(a):
        aa = a[np.ix_(dm, lm_)]
        ww = np.broadcast_to(cw[dm, :], aa.shape)
        return np.nansum(aa * ww) / np.nansum(ww)

    dlam = np.deg2rad(np.gradient(lon))[None, :]
    coslat = np.cos(np.deg2rad(lat))[:, None]
    ddx = lambda a: np.gradient(a, axis=1) / dlam / (AEARTH * coslat)

    def metrics(pick):
        s = slopes(pick)
        return dict(conv_E=box_conv(s["iu"], s["iv"], *BOX_E),
                    conv_W=box_conv(s["iu"], s["iv"], *BOX_W),
                    work=awmean(s["bu"] * s["iu"] + s["bv"] * s["iv"]),
                    prop=awmean(ddx(s["bu"]) * s["iu"] + ddx(s["bv"]) * s["iv"]))

    full = metrics(np.arange(NS))
    per = [metrics(np.array([k])) for k in range(NS)]

    # KE amplification timescale
    f = slopes(np.arange(NS))
    ke = awmean(f["bu"] ** 2 + f["bv"] ** 2)
    tau = ke / full["work"] / 2.0 if full["work"] > 0 else np.nan

    rng = np.random.default_rng(0)
    boot = {k: np.empty(NBOOT) for k in full}
    for b in range(NBOOT):
        m = metrics(rng.integers(0, NS, NS))
        for k in full:
            boot[k][b] = m[k]

    L = []
    P = L.append
    P("Increment metrics: what the CNN wind correction does to the MJO")
    P("=" * 78)
    P(f"source            : {os.path.basename(path)}")
    P(f"winters           : {NS} ({years[0]}--{years[-1]}), 180-day seasons from 1 Nov")
    P(f"layer             : 700-1000 hPa mass-weighted mean")
    P(f"index             : 20-70d band-pass precip, {IDX_BOX[0]:.0f}-{IDX_BOX[1]:.0f}E, "
      f"{LAT_IO[0]:.0f}-{LAT_IO[1]:.0f}N, standardized")
    P("")
    P("Per-winter values (each winter regressed on its own standardized index)")
    P(f"{'winter':>8} {'conv_E':>9} {'conv_W':>9} {'work':>11} {'prop (1e-8)':>13}")
    for y, m in zip(years, per):
        P(f"{y:>8} {m['conv_E']:+9.3f} {m['conv_W']:+9.3f} {m['work']:+11.4f} {m['prop']*1e8:+13.3f}")
    P("")
    P("Full-record value, exact two-sided binomial sign test over the 11 winters,")
    P("and leave-one-winter-out jackknife range")
    P("-" * 78)
    spec = [("conv_E", True,  "convergence east of convection (100-130E)  > 0"),
            ("conv_W", False, "divergence  west of convection (55-80E)    < 0"),
            ("work",   True,  "increment work on the MJO circulation      > 0"),
            ("prop",   False, "eastward displacement tendency             < 0")]
    for key, pos, label in spec:
        vals = [m[key] for m in per]
        k, p = sign_test(vals, pos)
        jk = [metrics(np.delete(np.arange(NS), i))[key] for i in range(NS)]
        lo, hi = np.percentile(boot[key], [5, 95])
        flag = "" if p < 0.05 else "   [NOT significant at 5%]"
        P(f"{label}")
        P(f"   full record = {full[key]:+.4g}")
        P(f"   sign test   = {k}/{NS} winters, p = {p:.4f}{flag}")
        P(f"   jackknife   = [{min(jk):+.4g}, {max(jk):+.4g}]  sign stable: "
          f"{all(np.sign(v) == np.sign(full[key]) for v in jk)}")
        P(f"   block boot  = 90% CI [{lo:+.4g}, {hi:+.4g}]  (n=11 blocks; anticonservative)")
        P("")
    P(f"MJO circulation variance <|V'|^2> = {ke:.4f} m2 s-2 per (index sd)^2")
    P(f"implied kinetic-energy e-folding amplification time = {tau:.1f} days")
    P("")
    P("Units: conv_E/conv_W in 1e-6 s^-1 day^-1; work in m2 s-2 day-1 per (index sd)^2;")
    P("prop in m s-2 day-1 per (index sd)^2. Box means use the Gauss theorem on the box")
    P("boundary, so interior grid-scale noise in the increment cancels.")

    txt = "\n".join(L)
    print(txt)
    out = os.path.join(FIGDIR, "increment_metrics.txt")
    with open(out, "w") as fh:
        fh.write(txt + "\n")
    print(f"\n[compute_increment_metrics] wrote {out}")


if __name__ == "__main__":
    main()
