#!/usr/bin/env python
"""
fig_heating_vstruct.py

Equatorial longitude-height cross-section of the MJO-regressed diabatic heating
structure, as a clean 2x2 comparison:
    a) ERA5   b) CNTRL   c) MEAN (MEAN)   d) CNN (CNN)

Shading : diabatic heating rate Q1 [K/day]        (mjo_diah_{model,obs}.nc)
Vectors : zonal-vertical circulation (U, -omega)  (mjo_U_ / mjo_W_{model,obs}.nc)

Clean 2x2 fork of Panel_Figures/Diabatic_Heating.ipynb (pyramid layout dropped,
"partial CNN" MITA1.0MJO experiment dropped). Each product carries the 5S-5N-
averaged, 20-70 day band-passed time-lon field plus the standardized Indian-Ocean
precip index (piom model / pio obs); fields are regressed onto that index here
and scaled to a +3 sigma MJO event.  Q1 x3; U x3; omega x3*100 (and negated for
the vertical component so up = ascent).

Run: module load conda && conda run --no-capture-output \
     -p /glade/u/apps/opt/conda/envs/npl-2024a python fig_heating_vstruct.py
Out: Paper_WCD/figures/fig_heating_vstruct.png
"""
import os
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from scipy.stats import linregress
import warnings
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NCDIR = os.path.join(ROOT, "nc")
FIGDIR = os.path.join(ROOT, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# ---- exact helpers copied from Panel_Figures/Diabatic_Heating.ipynb --------
def reg_coef_n(x, y, dims_x, dims_y):
    if x.shape[dims_x] != y.shape[dims_y]:
        raise ValueError("The specified dimensions of x and y must be of the same length.")
    if x.ndim > 1:
        x = np.moveaxis(x, dims_x, 0)
    if y.ndim > 1:
        y = np.moveaxis(y, dims_y, 0)
    if x.ndim > 1 and y.ndim > 1:
        other_dims = x.shape[1:]
    elif x.ndim > 1:
        other_dims = x.shape[1:]
    elif y.ndim > 1:
        other_dims = y.shape[1:]
    else:
        other_dims = ()
    slopes = np.zeros(other_dims)
    intercepts = np.zeros(other_dims)
    it = np.nditer(slopes, flags=['multi_index'])
    while not it.finished:
        idx = it.multi_index
        if x.ndim > 1 and y.ndim > 1:
            idx_x = (slice(None),) + idx
            slope, intercept, r, p, se = linregress(x[idx_x], y[idx_x])
        elif x.ndim > 1:
            idx_x = (slice(None),) + idx
            slope, intercept, r, p, se = linregress(x[idx_x], y)
        elif y.ndim > 1:
            idx_y = (slice(None),) + idx
            slope, intercept, r, p, se = linregress(x, y[idx_y])
        else:
            slope, intercept, r, p, se = linregress(x, y)
        slopes[idx] = slope
        intercepts[idx] = intercept
        it.iternext()
    return slopes, intercepts


def mb_to_km(levs):
    return ((1 - ((levs / 1013.25) ** 0.190284)) * 145366.45 * .3048) / 1000


# ---- shading colormap: exact ListedColormap from cell 7 --------------------
clevs = np.array([-0.9, -0.5, -0.1, 0.1, 0.5, 0.9, 1.3, 1.7, 2.1, 2.5])
newcolors = np.array([[0.31666667, 0.51372549, 0.73431373, 1.],
                      [0.5627451, 0.76470588, 0.86666667, 1.],
                      [1., 1., 1., 1.],
                      [0.99754902, 0.92401961, 0.63382353, 1.],
                      [0.99509804, 0.82941176, 0.51862745, 1.],
                      [0.99264706, 0.70686275, 0.40343137, 1.],
                      [0.9745098, 0.55490196, 0.32156863, 1.],
                      [0.91895425, 0.34771242, 0.22614379, 1.],
                      [0.85577342, 0.21481481, 0.16514161, 1.],
                      [0.75599129, 0.10457516, 0.15119826, 1.]])
cmapzzz = ListedColormap(newcolors)
cmapzzz.set_over((0.64705882, 0., 0.14901961, 1.))
cmapzzz.set_under((0.19215686, 0.21176471, 0.58431373, 1.))
bbox_props = dict(fc="white", ec="k", lw=1)

EXPS = {
    "CNTRL": "f.e.FTORCHmjo_CNTRLmjo_DT2",
    "MEAN":  "f.e.FTORCHmjo_MEANmjo_DT2",
    "CNN":   "f.e.FTORCHmjo_fullCNN_DT2",
}
PANELS = [
    ("a) ERA5",         "CNTRL", True),
    ("b) CNTRL",        "CNTRL", False),
    ("c) MEAN",      "MEAN",  False),
    ("d) CNN", "CNN",   False),
]
skip1 = slice(None, None, 7)


def fpath(tag, product):
    return os.path.join(NCDIR, tag, f"{EXPS[tag]}.{product}.nc")


def _reg(tag, product, ivar, fvar):
    ds = xr.open_dataset(fpath(tag, product))
    slopes, _ = reg_coef_n(ds[ivar].values, ds[fvar].values, 0, 0)
    return slopes, ds["lon"].values, ds["level"].values


def load_panel(tag, is_obs):
    """Return (lons, levels, Q1 [K/day], U [m/s], omega [hPa/day-ish scaled])."""
    if is_obs:
        ivar, fvar = "pio", "THETAE"
        q1, lons, levels = _reg(tag, "mjo_diah_obs", "pio", "q1")
        w, _, _ = _reg(tag, "mjo_W_obs", ivar, fvar)
        u, _, _ = _reg(tag, "mjo_U_obs", ivar, fvar)
    else:
        ivar, fvar = "piom", "THETAEm"
        q1, lons, levels = _reg(tag, "mjo_diah_model", "piom", "q1m")
        w, _, _ = _reg(tag, "mjo_W_model", ivar, fvar)
        u, _, _ = _reg(tag, "mjo_U_model", ivar, fvar)
    xxs = u * 3
    yys = w * 3 * 100
    return lons, levels, q1 * 3, xxs, yys


def main():
    fig = plt.figure(figsize=(11, 7))
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.14, hspace=0.20,
                           left=0.075, right=0.985, top=0.965, bottom=0.135)
    axes = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
    f = None
    for ax, (label, tag, is_obs) in zip(axes, PANELS):
        lons, levels, q1, xxs, yys = load_panel(tag, is_obs)
        yy = mb_to_km(levels)
        f = ax.contourf(lons, yy, q1, levels=clevs, cmap=cmapzzz, extend='both')
        ax.quiver(lons[skip1], yy, xxs[:, skip1], -1 * yys[:, skip1],
                  pivot="middle", color='gray', scale=50, alpha=0.45,
                  width=0.0035, headwidth=5)
        ax.set_xlim([10, 180])
        ax.set_xticks([40, 80, 120, 160])
        ax.set_xticklabels([f"{v}°E" for v in [40, 80, 120, 160]])
        ax.grid(True, alpha=0.3)
        ax.set_yticks(yy[::2])
        ax.set_yticklabels(list(levels.astype(int))[::2])
        ax.set_ylim([0, yy[-1]])
        ax.tick_params(labelsize=14)
        ax.plot([90, 90], [0, 20], color='k', linestyle=':')
        ax.text(0.98, 0.97, label, transform=ax.transAxes, ha='right', va='top',
                fontsize=16, bbox=bbox_props)
    for ax in (axes[1], axes[3]):
        plt.setp(ax.get_yticklabels(), visible=False)
    for ax in (axes[0], axes[1]):
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in (axes[0], axes[2]):
        ax.set_ylabel("hPa", fontsize=16)

    cbar_ax = fig.add_axes([0.30, 0.055, 0.42, 0.02])
    cbar = fig.colorbar(f, cax=cbar_ax, orientation='horizontal', ticks=clevs)
    cbar.ax.tick_params(labelsize=13)
    cbar.set_label(label='[K/day]', fontsize=15)
    out = os.path.join(FIGDIR, "fig_heating_vstruct.png")
    fig.savefig(out, bbox_inches='tight', dpi=400)
    print("wrote", out)


if __name__ == "__main__":
    main()
