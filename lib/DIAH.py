import os
import numpy as np
import xarray as xr
import cftime

# Define constants
cp = 1004.7
R = 287.0

# Define region
latS_IO, latN_IO = -10.0, 10.0

# Define plot parameters
modout = os.getenv("modout")

# Define data directories and variables
uName, vName, wName, aName = "U", "V", "W", "T"

# Time periods
twStrt, twLast = 19790101, 19951231
cwStrtm, cwLastm = 19790101, 19901231

uNamem, aNamem, vNamem, wNamem = "U", "T", "V", "OMEGA"

# Convert time periods to datetime objects
twStrt_dt = cftime.DatetimeNoLeap(twStrt // 10000, (twStrt % 10000) // 100, twStrt % 100)
twLast_dt = cftime.DatetimeNoLeap(twLast // 10000, (twLast % 10000) // 100, twLast % 100)
cwStrtm_dt = cftime.DatetimeNoLeap(cwStrtm // 10000, (cwStrtm % 10000) // 100, cwStrtm % 100)
cwLastm_dt = cftime.DatetimeNoLeap(cwLastm // 10000, (cwLastm % 10000) // 100, cwLastm % 100)

# Constants
R = 287.0  # Specific gas constant for dry air [J/(kg·K)]
cp = 1004.0  # Specific heat at constant pressure for dry air [J/(kg·K)]

# Function to calculate Q1
# Function to calculate Q1
def calculate_q1(AIR, U, V, W, lon, lat, lev):
    print("air shape:", AIR.shape)
    
    # Calculate dTdt (K/s)
    dTdt = np.gradient(AIR, axis=0) / 86400.0  # Note: Dividing by 86400 to convert from K/day to K/s

    # Calculate dtdx (K/s)
    dlon = np.deg2rad(lon[1] - lon[0])
    dtdx = np.zeros_like(AIR)
    for nl in range(len(lat)):
        dX = 6378388.0 * np.cos(np.deg2rad(lat[nl])) * dlon
        dtdx[:, :, nl:nl+1, :] = np.gradient(AIR[:, :, nl:nl+1, :], axis=3) / dX

    # Calculate utx (K/s)
    utx = U * dtdx

    # Calculate dtdy (K/s)
    dlat = np.deg2rad(lat)
    dY = 6378388.0 * dlat
    dY = dY.reshape((1, 1, len(dlat), 1))  # Reshape dY for broadcasting
    dtdy = np.gradient(AIR, axis=2) / dY

    # Calculate vty (K/s)
    vty = V * dtdy

    # Calculate dtdp (K/s)
    dtdp = np.gradient(AIR, axis=1) / (lev[:, None, None] * 100.0)

    # Calculate wtp (K/s)
    wtp = W * dtdp

    # Calculate wt (K/s)
    wt = W * AIR * R / cp
    for i in range(len(lev)):
        wt[:, i, :, :] /= lev[i] / 100.0

    # Calculate Q1 (K/s)
    q1 = dTdt + utx + vty + wtp - wt
    q1 = q1 * 86400  # Convert back to K/day

    return q1, dTdt, utx, vty, wtp, wt



# Open the datasets using xarray
obs_data_t_path = os.getenv("obs_data_t")
obs_data_u_path = os.getenv("obs_data_u")
obs_data_v_path = os.getenv("obs_data_v")
obs_data_w_path = os.getenv("obs_data_w")
model_file_path = os.getenv("mod_data_daily")

ds_t = xr.open_dataset(obs_data_t_path)
ds_u = xr.open_dataset(obs_data_u_path)
ds_v = xr.open_dataset(obs_data_v_path)
ds_w = xr.open_dataset(obs_data_w_path)
ds_model = xr.open_dataset(model_file_path)

# Select data within the desired time range
ds_t_sel = ds_t.sel(time=slice(str(twStrt_dt)[:10], str(twLast_dt)[:10]), lat=slice(latS_IO, latN_IO))
ds_u_sel = ds_u.sel(time=slice(str(twStrt_dt)[:10], str(twLast_dt)[:10]), lat=slice(latS_IO, latN_IO))
ds_v_sel = ds_v.sel(time=slice(str(twStrt_dt)[:10], str(twLast_dt)[:10]), lat=slice(latS_IO, latN_IO))
ds_w_sel = ds_w.sel(time=slice(str(twStrt_dt)[:10], str(twLast_dt)[:10]), lat=slice(latS_IO, latN_IO))
ds_model_sel = ds_model.sel(time=slice(str(cwStrtm_dt)[:10], str(cwLastm_dt)[:10]), lat=slice(latS_IO, latN_IO))

# Extract variables from datasets
AIR = ds_t_sel[aName].values
U = ds_u_sel[uName].values
V = ds_v_sel[vName].values
W = ds_w_sel[wName].values
lon = ds_t_sel['lon'].values
lat = ds_t_sel['lat'].values
lev = ds_t_sel['level'].values
time_obs = ds_t_sel['time'].values

AIRm = ds_model_sel[aNamem].values
Um = ds_model_sel[uNamem].values
Vm = ds_model_sel[vNamem].values
Wm = ds_model_sel[wNamem].values
lonm = ds_model_sel['lon'].values
latm = ds_model_sel['lat'].values
levm = ds_model_sel['level'].values
time_model = ds_model_sel['time'].values

lev_out = np.array([925, 700, 300])

# Calculate Q1 for OBS and MODEL data
q1_obs, dTdt_obs, utx_obs, vty_obs, wtp_obs, wt_obs = calculate_q1(AIR, U, V, W, lon, lat, lev)

lev_idx = np.isin(lev, lev_out).nonzero()[0]

q1_obs_out = q1_obs[:, lev_idx, :, :]
dTdt_obs_out = dTdt_obs[:, lev_idx, :, :]
utx_obs_out = utx_obs[:, lev_idx, :, :]
vty_obs_out = vty_obs[:, lev_idx, :, :]
wtp_obs_out = wtp_obs[:, lev_idx, :, :]
wt_obs_out = wt_obs[:, lev_idx, :, :]

q1_model, dTdt_model, utx_model, vty_model, wtp_model, wt_model = calculate_q1(AIRm, Um, Vm, Wm, lonm, latm, levm)
levm_idx = np.isin(levm, lev_out).nonzero()[0]

q1_model_out = q1_model[:, levm_idx, :, :]
dTdt_model_out = dTdt_model[:, levm_idx, :, :]
utx_model_out = utx_model[:, levm_idx, :, :]
vty_model_out = vty_model[:, levm_idx, :, :]
wtp_model_out = wtp_model[:, levm_idx, :, :]
wt_model_out = wt_model[:, levm_idx, :, :]

# Create xarray Datasets for Q1 outputs with additional variables
q1_obs_ds = xr.Dataset(
    {
        "q1_ll": (["time", "level", "lat", "lon"], q1_obs_out),
        "dTdt_ll": (["time", "level", "lat", "lon"], dTdt_obs_out),
        "utx_ll": (["time", "level", "lat", "lon"], utx_obs_out),
        "vty_ll": (["time", "level", "lat", "lon"], vty_obs_out),
        "wtp_ll": (["time", "level", "lat", "lon"], wtp_obs_out),
        "wt_ll": (["time", "level", "lat", "lon"], wt_obs_out)
    },
    coords={
        "time": ds_t_sel['time'],
        "level": lev_out,
        "lat": ds_t_sel['lat'],
        "lon": ds_t_sel['lon']
    }
)

q1_model_ds = xr.Dataset(
    {
        "q1m_ll": (["time", "level", "lat", "lon"], q1_model_out),
        "dTdtm_ll": (["time", "level", "lat", "lon"], dTdt_model_out),
        "utxm_ll": (["time", "level", "lat", "lon"], utx_model_out),
        "vtym_ll": (["time", "level", "lat", "lon"], vty_model_out),
        "wtpm_ll": (["time", "level", "lat", "lon"], wtp_model_out),
        "wtm_ll": (["time", "level", "lat", "lon"], wt_model_out)
    },
    coords={
        "time": ds_model_sel['time'],
        "level": lev_out,
        "lat": ds_model_sel['lat'],
        "lon": ds_model_sel['lon']
    }
)

# Save OBS Q1 to NetCDF using xarray
obs_out_file = f"./NC_OUT/{modout}.mjo_diah2d_ll_obs_py.nc"
q1_obs_ds.to_netcdf(obs_out_file)

# Save MODEL Q1 to NetCDF using xarray
model_out_file = f"./NC_OUT/{modout}.mjo_diah2d_ll_model_py.nc"
q1_model_ds.to_netcdf(model_out_file)
