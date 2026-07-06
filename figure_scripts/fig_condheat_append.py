#!/usr/bin/env python
"""
fig_condheat_append.py  (appendix figure)

Equatorial longitude-height cross-section of the MJO-regressed *condensational*
temperature tendency DTCOND (CAM "T tendency - moist processes" = total heating
from ZM deep convection + shallow convection + large-scale/stratiform
condensation, i.e. the model's parameterized moist heating, NOT deep convection
alone), for CONTROL, MEAN, and the CNN-corrected simulation, plus the two
physically distinct differences that carry the paper's central logic:
    * CNN - MEAN  = the STATE-DEPENDENT effect (the primary contrast: what the
      flow-dependent CNN correction does beyond its time-mean),
    * MEAN - CNTRL = the CLIMATOLOGICAL effect (what the time-mean correction
      alone does relative to the uncorrected control).

DTCOND is read from the ADF model files and put through the *identical*
mjo_ape.ncl pipeline used for the Q1 heating figure, then regressed onto the same
standardized precip index (piom) stored in mjo_ape_model_all.nc:
    * meridional average 5S-5N
    * 20-70 day Lanczos band-pass (nWgt=141), matching filwgts_lanczos
    * winter seasons: 180-day blocks starting 1 Nov, concatenated
    * regression onto standardized piom, scaled to a +3 sigma MJO event
    * K/s -> K/day

CNTRL and CNN carry DTCOND on the 13 standard pressure levels in their
`plev.mjo` ADF files. The MEAN experiment's `plev.mjo` file does NOT contain
DTCOND (its vstruct run never processed it), so MEAN DTCOND is reconstructed from
the raw model-level ADF file `*.cam.h1.DTCOND.*.nc` by hybrid->pressure
interpolation (geocat interp_hybrid_to_pressure, log method) onto the same 13
pressure levels using the file's bundled hyam/hybm/P0/PS -- the same operation
ADF/vinth2p performs, so this is a regeneration of the identical field, not
synthetic data.

CNN canonical = fullCNN DT2 (period-matched to CNTRL DT2). DTCOND is present in
the DT2 ADF files, so no DT8 fallback is needed.

Run: conda run -p /glade/u/apps/opt/conda/envs/npl-2024a python fig_condheat_append.py
Out: Paper_WCD/figures/fig_condheat_append.png
"""
import os
import sys
import glob
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import warnings
warnings.filterwarnings("ignore")

WANG = "/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang"
sys.path.insert(0, WANG)
from lanczos_filter import lanczos_lowpass_weights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NCDIR = os.path.join(ROOT, "nc")
FIGDIR = os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

ADF_DIR = "/glade/campaign/cgd/amp/wchapman/ADF/{tag}/ts/"


def find_adf(tag):
    pat = os.path.join(ADF_DIR.format(tag=tag),
                       f"{tag}.cam.h1.NC4_Classic.plev.mjo.*.nc")
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(pat)
    return hits[0]


def find_raw_dtcond(tag):
    """Raw (model hybrid-level) DTCOND time series, with hyam/hybm/P0/PS."""
    pat = os.path.join(ADF_DIR.format(tag=tag), f"{tag}.cam.h1.DTCOND.*.nc")
    hits = sorted(glob.glob(pat))
    if not hits:
        raise FileNotFoundError(pat)
    return hits[0]


EXPS = {
    "CNTRL": ("CNTRL", "f.e.FTORCHmjo_CNTRLmjo_DT2"),
    "MEAN":  ("MEAN",  "f.e.FTORCHmjo_MEANmjo_DT2"),
    "CNN":   ("CNN",   "f.e.FTORCHmjo_fullCNN_DT2"),
}
# Experiments whose plev.mjo file lacks DTCOND -> interpolate from the raw file.
HYBRID_EXPS = {"MEAN"}
LATS, LATN = -5.0, 5.0
NDAY, MMDD_START = 180, 1101   # winter season: 180 days from 1 Nov
NWGT = 141
LON0, LON1 = 40, 180
MJO_CENTER = 90.0
# 13 standard pressure levels used in the CNTRL/CNN plev.mjo files (hPa).
PLEV_HPA = np.array([80., 100., 150., 200., 250., 300., 400.,
                     500., 600., 700., 850., 925., 1000.])


def bandpass_weights():
    # 20-70 day band pass = lowpass(1/20) - lowpass(1/70); matches NCL filwgts_lanczos
    lp_hi = lanczos_lowpass_weights(NWGT, 1.0 / 20.0).data
    lp_lo = lanczos_lowpass_weights(NWGT, 1.0 / 70.0).data
    return lp_hi - lp_lo


def load_dtcond_plev(tag):
    """DTCOND already on pressure levels in the plev.mjo ADF file (CNTRL/CNN)."""
    ds = xr.open_dataset(find_adf(tag))
    d = ds["DTCOND"].sel(lat=slice(LATS, LATN)).mean("lat")  # (time, level, lon)
    a = d.values * 86400.0                                   # K/day
    lev = ds["level"].values
    lon = ds["lon"].values
    tt = xr.decode_cf(ds[["time"]]).time
    return a, lon, lev, tt.dt.month.values, tt.dt.day.values


def load_dtcond_hybrid(tag, chunk=500):
    """DTCOND from the raw model-level ADF file, interpolated to PLEV_HPA.

    The MEAN plev.mjo file has no DTCOND, so reconstruct the pressure-level
    field from the native hybrid-level file the same way ADF/vinth2p does
    (log-linear hybrid->pressure using the bundled hyam/hybm/P0/PS).
    """
    import geocat.comp as gc
    ds = xr.open_dataset(find_raw_dtcond(tag)).sel(lat=slice(LATS, LATN))
    plev_pa = PLEV_HPA * 100.0
    p0 = float(ds["P0"])
    nt = ds.sizes["time"]
    outs = []
    for s in range(0, nt, chunk):
        sub = ds.isel(time=slice(s, s + chunk))
        op = gc.interp_hybrid_to_pressure(
            sub["DTCOND"], sub["PS"], sub["hyam"], sub["hybm"],
            p0=p0, new_levels=plev_pa, method="log")      # (time, plev, lat, lon)
        outs.append(op.mean("lat").values)                # (time, plev, lon)
    a = np.concatenate(outs, axis=0) * 86400.0            # K/day
    lon = ds["lon"].values
    tt = xr.decode_cf(ds[["time"]]).time
    return a, lon, PLEV_HPA.copy(), tt.dt.month.values, tt.dt.day.values


def process_dtcond(folder, tag):
    if folder in HYBRID_EXPS:
        a, lon, lev, mm, dd = load_dtcond_hybrid(tag)
    else:
        a, lon, lev, mm, dd = load_dtcond_plev(tag)

    # band-pass filter along time
    w = bandpass_weights()
    filt = np.apply_along_axis(lambda s: np.convolve(s, w, mode="same"), 0, a)

    # winter blocks: 180-day segments starting on 1 Nov (mm=11, dd=01)
    starts = np.where((mm == 11) & (dd == 1))[0]
    starts = [s for s in starts if s + NDAY <= a.shape[0]]
    blocks = [filt[s:s + NDAY] for s in starts]
    win = np.concatenate(blocks, axis=0)                    # (nsea*180, lev, lon)

    # stored precip index (already NCL-filtered, winter-concatenated, standardized)
    ap = xr.open_dataset(os.path.join(NCDIR, folder, f"{tag}.mjo_ape_model_all.nc"))
    piom = ap["piom"].values.astype(float)
    piom = np.where(piom > 1e30, np.nan, piom)
    n = min(len(piom), win.shape[0])
    piom, win = piom[:n], win[:n]

    # standardize index, regress field onto it, scale to +3 sigma event
    good = np.isfinite(piom)
    x = piom[good]
    x = (x - x.mean()) / x.std()
    Y = win[good]                                           # (n, lev, lon)
    Ya = Y - np.nanmean(Y, axis=0)
    slope = (x[:, None, None] * Ya).sum(0) / (x * x).sum()  # (lev, lon)
    reg = slope * 3.0
    # mild longitudinal 1-2-1 smoothing for presentation
    k = np.array([1.0, 2.0, 1.0]); k /= k.sum()
    reg = np.apply_along_axis(lambda s: np.convolve(s, k, mode="same"), 1, reg)
    return reg, lon, lev, len(starts)


def main():
    fields = {}
    lon = lev = None
    for name, (folder, tag) in EXPS.items():
        sl, lon, lev, nsea = process_dtcond(folder, tag)
        fields[name] = sl
        print(f"{name}: nseasons={nsea}  peak={np.nanmax(sl):.2f}  "
              f"min={np.nanmin(sl):.2f} K/day")

    # State-dependent (primary) and climatological differences.
    diff_cnn_mean = fields["CNN"] - fields["MEAN"]      # state-dependent effect
    diff_mean_cntrl = fields["MEAN"] - fields["CNTRL"]  # climatological effect
    # (CNN - CNTRL remains available if needed for supplementary comparison.)
    diff_cnn_cntrl = fields["CNN"] - fields["CNTRL"]    # noqa: F841

    m = (lon >= LON0) & (lon <= LON1)
    dcm = diff_cnn_mean[:, m]
    print("CNN-MEAN diff: peak={:.3f} at lev={:.0f}hPa lon={:.0f}E ; "
          "min={:.3f} at lev={:.0f}hPa lon={:.0f}E".format(
              np.nanmax(dcm),
              lev[np.unravel_index(np.nanargmax(dcm), dcm.shape)[0]],
              lon[m][np.unravel_index(np.nanargmax(dcm), dcm.shape)[1]],
              np.nanmin(dcm),
              lev[np.unravel_index(np.nanargmin(dcm), dcm.shape)[0]],
              lon[m][np.unravel_index(np.nanargmin(dcm), dcm.shape)[1]]))

    # Top row: absolute regressed DTCOND. Bottom row: the two differences.
    panels = [
        (0, 0, "a) CNTRL", fields["CNTRL"], "heat"),
        (0, 1, "b) MEAN", fields["MEAN"], "heat"),
        (0, 2, "c) CNN", fields["CNN"], "heat"),
        (1, 0, "d) CNN $-$ MEAN", diff_cnn_mean, "diff"),
        (1, 1, "e) MEAN $-$ CNTRL", diff_mean_cntrl, "diff"),
    ]
    clevs_h = np.array([-0.9, -0.5, -0.1, 0.3, 0.7, 1.1, 1.5, 1.9, 2.3, 2.7, 3.1])
    norm_h = TwoSlopeNorm(vmin=-0.9, vcenter=0.1, vmax=3.1)
    clevs_d = np.linspace(-0.8, 0.8, 17)

    yticks = [1000, 850, 600, 400, 250, 150]
    xticks = np.arange(40, 181, 40)

    fig = plt.figure(figsize=(14, 8.2))
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.22, hspace=0.30,
                           left=0.06, right=0.985, top=0.95, bottom=0.10)
    cf_h = cf_d = None
    for r, c, label, arr, kind in panels:
        ax = fig.add_subplot(gs[r, c])
        if kind == "heat":
            cf_h = ax.contourf(lon[m], lev, arr[:, m], levels=clevs_h,
                               cmap="RdBu_r", norm=norm_h, extend="both")
        else:
            cf_d = ax.contourf(lon[m], lev, arr[:, m], levels=clevs_d,
                               cmap="RdBu_r", extend="both")
        ax.axvline(MJO_CENTER, ls=":", color="k", lw=1.3)
        ax.set_ylim(1000, 100)
        ax.set_xlim(LON0, LON1)
        ax.set_yticks(yticks); ax.set_yticklabels([str(v) for v in yticks])
        ax.set_xticks(xticks); ax.set_xticklabels([f"{int(v)}°E" for v in xticks])
        ax.tick_params(labelsize=14)
        ax.grid(True, alpha=0.25)
        if c == 0:
            ax.set_ylabel("hPa", fontsize=16)
        ax.text(0.96, 0.94, label, transform=ax.transAxes, ha="right", va="top",
                fontsize=16, bbox=dict(fc="white", ec="k", lw=1.0))

    # Two vertical colorbars occupying the empty bottom-right cell (gs[1, 2]).
    pos = gs[1, 2].get_position(fig)
    cax1 = fig.add_axes([pos.x0 + 0.010, pos.y0 + 0.04,
                         0.016, pos.height - 0.08])
    cb1 = fig.colorbar(cf_h, cax=cax1, orientation="vertical", ticks=clevs_h)
    cb1.set_label(r"DTCOND regressed  [K day$^{-1}$]", fontsize=15)
    cb1.ax.tick_params(labelsize=13)
    cax2 = fig.add_axes([pos.x0 + pos.width * 0.55, pos.y0 + 0.04,
                         0.016, pos.height - 0.08])
    cb2 = fig.colorbar(cf_d, cax=cax2, orientation="vertical", ticks=clevs_d[::2])
    cb2.set_label(r"difference  [K day$^{-1}$]", fontsize=15)
    cb2.ax.tick_params(labelsize=13)

    # no figure title: the manuscript caption serves this role
    out = os.path.join(FIGDIR, "fig_condheat_append.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
