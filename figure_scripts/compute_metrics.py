import os, numpy as np, xarray as xr
from make_panels import load_slopes, EXPERIMENTS, NCDIR
from scipy.stats import linregress

# ---- Kelvin-Rossby ratio from regressed U850 (equatorial 10S-10N, 40-180E) ----
print("=== Kelvin-Rossby amplitude ratio (850 hPa U regression) ===")
print("def: over 10S-10N, 40-180E; Rossby=max westerly(U>0), Kelvin=max easterly(|U<0|); ratio=Kelvin/Rossby")
for exp in EXPERIMENTS:
    suf = "mjo_U850_obs_all.850mb.nc" if exp["obs"] else "mjo_U850_model_all.850mb.80.100.nc"
    path = os.path.join(NCDIR, exp["dir"], f"{exp['modout']}.{suf}")
    sl, lon, lat, fld = load_slopes(path, ["u_timem","u_time"])
    la = (lat>=-10)&(lat<=10); lo=(lon>=40)&(lon<=180)
    sub = sl[np.ix_(la,lo)]
    west = np.nanmax(sub); east = np.nanmax(-sub)
    print(f"  {exp['label']:20s} Rossby(westerly)={west:5.2f}  Kelvin(easterly)={east:5.2f}  K/R={east/west:4.2f}")

# ---- eastward/westward MJO-band power ratio from spectra cache ----
print("\n=== MJO-band spectral power (from fig_spectra_cache.npz if present) ===")
cache = os.path.join(os.path.dirname(NCDIR),"figures","fig_spectra_cache.npz")
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True)
    print("cache keys:", list(z.keys()))
    for k in z.keys():
        arr = z[k]
        print(f"  {k}: shape {getattr(arr,'shape',None)} dtype {getattr(arr,'dtype',None)}")
else:
    print("  no spectra cache found at", cache)
