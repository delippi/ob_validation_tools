#!/usr/bin/env python3
#
# Full H(x) comparison tool with KD-tree matching (stationID aware).
#
# Matches by:
#   stationID (exact match required)  [when available in BOTH files]
#   lat/lon within tolerance
#   pressure [hPa] or height [m] within tolerance
#   datetime within tolerance
#
# Panels:
#   1. delta H(x) map
#   2. delta observation-error map
#   3. vertical-coordinate map
#   4. time-offset histogram
#   5. delta H(x) versus vertical coordinate
#   6. delta observation error versus vertical coordinate
#   7. H(x) one-to-one comparison
#   8. observation-error one-to-one comparison
#
# Pressure is used for conventional observations. Height is used for
# height-based observations such as nexrad_radialVelocity_999.
#

import argparse
import numpy as np
from netCDF4 import Dataset

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import warnings
warnings.filterwarnings("ignore")

# For KD-tree
from scipy.spatial import cKDTree

# ---------------------------------------------------------------
# Unified threshold logic for dh and d_err
# ---------------------------------------------------------------
def get_default_thresholds(varname):
    if varname in ("airTemperature", "airTemperatureAt2M"):
        return 0.5, 0.1       # (dh,K), (err,K)
    elif varname in ("windEastward", "windNorthward", "radialVelocity"):
        return 0.5, 0.1       # m/s
    elif varname in ("specificHumidity", "specificHumidityAt2M"):
        return 0.0005, 0.001  # kg/kg
    elif varname == "stationPressure":
        return 0.5, 0.5       # hPa
    else:
        return 0.5, 0.1       # default fallback values

def resolve_thresholds(varname, user_dh, user_derr):
    dh_def, derr_def = get_default_thresholds(varname)
    dh_final   = user_dh   if user_dh   is not None else dh_def
    derr_final = user_derr if user_derr is not None else derr_def
    return dh_final, derr_final


# ---------------------------------------------------------------
# Convert a netCDF variable to a floating-point NumPy array
# ---------------------------------------------------------------
def as_float_array(values):
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=float)


# ---------------------------------------------------------------
# Select the vertical coordinate used for matching and plotting
# ---------------------------------------------------------------
def read_vertical_coordinate(meta, varname, nloc):
    # NEXRAD observations are height-based, so prefer height for radialVelocity.
    if varname == "radialVelocity":
        candidates = (
            "height",
            "heightAboveSeaLevel",
            "heightAboveMeanSeaLevel",
            "pressure",
        )
    else:
        candidates = (
            "pressure",
            "height",
            "heightAboveSeaLevel",
            "heightAboveMeanSeaLevel",
        )

    for name in candidates:
        if name not in meta.variables:
            continue

        variable = meta.variables[name]
        values = as_float_array(variable[:])

        if name == "pressure":
            # IODA pressure metadata normally uses Pa. Preserve hPa input if a
            # diagnostic already declares those units.
            input_units = str(getattr(variable, "units", "Pa")).strip().lower()
            if input_units not in ("hpa", "hectopascal", "hectopascals"):
                values = values * 0.01

            return values, {
                "name": name,
                "kind": "pressure",
                "label": "Pressure (hPa)",
                "units": "hPa",
                "available": True,
                "invert_axis": True,
            }

        return values, {
            "name": name,
            "kind": "height",
            "label": "Height (m)",
            "units": "m",
            "available": True,
            "invert_axis": False,
        }

    # Keep support for observation types with no vertical coordinate. A zero
    # placeholder disables vertical discrimination but retains lat/lon/time
    # matching, which is the behavior used previously for radiances.
    return np.zeros(nloc, dtype=float), {
        "name": None,
        "kind": "none",
        "label": "Vertical coordinate unavailable",
        "units": "",
        "available": False,
        "invert_axis": False,
    }


# ---------------------------------------------------------------
# Read station IDs from either IODA strings or fixed-width char arrays
# ---------------------------------------------------------------
def read_station_ids(meta, nloc):
    if "stationIdentification" not in meta.variables:
        return np.array(["NO_STATION_ID"] * nloc, dtype=object), False

    raw = meta.variables["stationIdentification"][:]
    if np.ma.isMaskedArray(raw):
        raw = raw.filled("")
    raw = np.asarray(raw)

    # Some NEXRAD files store four-character radar IDs as an nloc x nchar
    # array instead of a one-dimensional netCDF string variable.
    if raw.ndim > 1:
        values = []
        for row in raw:
            if raw.dtype.kind == "S":
                value = b"".join(np.asarray(row).tolist()).decode(
                    "utf-8", errors="replace"
                )
            else:
                value = "".join(str(item) for item in np.asarray(row).tolist())
            values.append(value)
    else:
        values = raw.tolist()

    station_ids = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        station_ids.append(str(value).replace("\x00", "").strip())

    has_station_ids = (
        len(station_ids) == nloc
        and all(station_id for station_id in station_ids)
    )

    if not has_station_ids:
        station_ids = ["NO_STATION_ID"] * nloc

    return np.asarray(station_ids, dtype=object), has_station_ids


# ---------------------------------------------------------------
# Check valid IODA file
# ---------------------------------------------------------------
def is_ioda_diag(nc):
    return (
        "MetaData" in nc.groups
        and "_ioda_layout" in nc.ncattrs()
        and "hofx0" in nc.groups
        and "EffectiveQC2" in nc.groups
    )

# ---------------------------------------------------------------
# Load diag and group obs (optionally by stationID)
# ---------------------------------------------------------------
def load_by_station(diagfile, varname, hofx_group):
    with Dataset(diagfile, "r") as nc:

        if not is_ioda_diag(nc):
            raise RuntimeError("Not an IODA diag: " + diagfile)

        if hofx_group not in nc.groups:
            raise RuntimeError(
                f"Missing H(x) group '{hofx_group}' in {diagfile}"
            )

        meta = nc.groups["MetaData"]
        qcgrp = nc.groups["EffectiveQC2"]
        hofxgrp = nc.groups[hofx_group]
        errgrp = nc.groups["EffectiveError0"]

        for group_name, group in (
            ("EffectiveQC2", qcgrp),
            (hofx_group, hofxgrp),
            ("EffectiveError0", errgrp),
        ):
            if varname not in group.variables:
                raise RuntimeError(
                    f"Missing variable '{varname}' in group "
                    f"'{group_name}' in {diagfile}"
                )

        lat = as_float_array(meta.variables["latitude"][:])
        lon = as_float_array(meta.variables["longitude"][:])
        nloc = len(lat)

        # Select pressure for conventional data or height for NEXRAD data.
        vert, vertical_info = read_vertical_coordinate(meta, varname, nloc)

        dt = as_float_array(meta.variables["dateTime"][:])

        # read timeOffset if present
        if "timeOffset" in meta.variables:
            toff = as_float_array(meta.variables["timeOffset"][:])
        else:
            # NaN keeps missing offsets out of the histogram.
            toff = np.full(nloc, np.nan, dtype=float)

        # stationIdentification may be missing
        stn, has_station_ids = read_station_ids(meta, nloc)

        # wrap longitudes
        lon = np.where(lon > 180.0, lon - 360.0, lon)

        qc = as_float_array(qcgrp.variables[varname][:])
        hofx = as_float_array(hofxgrp.variables[varname][:])
        err = as_float_array(errgrp.variables[varname][:])
        print("")
        print("Error min/max")
        finite_err = err[np.isfinite(err)]
        if finite_err.size:
            print(np.min(finite_err), np.max(finite_err))
        else:
            print("No finite observation errors")
        print(f"{diagfile}")
        print(
            "Vertical coordinate:",
            vertical_info["name"] if vertical_info["available"] else "none",
        )
        print("")

        # Exclude locations that cannot participate in a KD-tree match.
        mask = (
            (qc == 0)
            & np.isfinite(lat)
            & np.isfinite(lon)
            & np.isfinite(vert)
            & np.isfinite(dt)
        )

        data = {}
        idx = np.where(mask)[0]

        for i in idx:
            s = stn[i]
            if s not in data:
                data[s] = {
                    "lat": [],
                    "lon": [],
                    "vert": [],
                    "time": [],
                    "timeOffset": [],
                    "hofx": [],
                    "error": []
                }
            data[s]["lat"].append(float(lat[i]))
            data[s]["lon"].append(float(lon[i]))
            data[s]["vert"].append(float(vert[i]))
            data[s]["time"].append(int(dt[i]))
            data[s]["timeOffset"].append(float(toff[i]))

            # Convert masked values or 0-D arrays into scalars
            hi = hofx[i]
            ei = err[i]
            if hasattr(hi, "filled"):
                hi = hi.filled(np.nan)   # replace masked with NaN
            if hasattr(ei, "filled"):
                ei = ei.filled(np.nan)   # replace masked with NaN
            try:
                hi = float(hi)
            except Exception:
                hi = np.nan
            try:
                ei = float(ei)
            except Exception:
                ei = np.nan
            data[s]["hofx"].append(hi)
            data[s]["error"].append(ei)

        # convert lists to numpy arrays
        for s in data:
            for key in data[s]:
                data[s][key] = np.array(data[s][key])

        return data, has_station_ids, vertical_info

# ---------------------------------------------------------------
# Load domain outline
# ---------------------------------------------------------------
def load_domain_ring(gridfile):
    ds = Dataset(gridfile, "r")

    if "grid_latt" in ds.variables:
        lat = ds.variables["grid_latt"][:]
        lon = ds.variables["grid_lont"][:]
    elif "grid_lat" in ds.variables:
        lat = ds.variables["grid_lat"][:]
        lon = ds.variables["grid_lon"][:]
    else:
        raise RuntimeError("Grid missing grid_latt/grid_lont or grid_lat/grid_lon")

    top    = np.c_[lon[0,:],        lat[0,:]]
    right  = np.c_[lon[1:,-1],      lat[1:,-1]]
    bottom = np.c_[lon[-1,-2::-1],  lat[-1,-2::-1]]
    left   = np.c_[lon[-2:0:-1,0],  lat[-2:0:-1,0]]

    ring = np.vstack([top, right, bottom, left])
    ring[:,0] = np.where(ring[:,0] > 180, ring[:,0] - 360, ring[:,0])
    ds.close()
    return ring

# ---------------------------------------------------------------
# Helper to flatten station-grouped dict to arrays
# ---------------------------------------------------------------
def flatten_obs_dict(d):
    all_lat = []
    all_lon = []
    all_vert = []
    all_time = []
    all_hofx = []
    all_err = []

    for s in d:
        all_lat.append(d[s]["lat"])
        all_lon.append(d[s]["lon"])
        all_vert.append(d[s]["vert"])
        all_time.append(d[s]["time"])
        all_hofx.append(d[s]["hofx"])
        all_err.append(d[s]["error"])

    if not all_lat:
        return (
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
            np.array([]),
        )

    return (
        np.concatenate(all_lat),
        np.concatenate(all_lon),
        np.concatenate(all_vert),
        np.concatenate(all_time),
        np.concatenate(all_hofx),
        np.concatenate(all_err),
    )

# ---------------------------------------------------------------
# Return an empty matching result with the normal tuple structure
# ---------------------------------------------------------------
def empty_match_result():
    return tuple(np.array([]) for _ in range(10))


# ---------------------------------------------------------------
# Match two flat observation subsets with a normalized KD-tree
# ---------------------------------------------------------------
def match_flat_arrays(g_arrays, j_arrays, tolerances, station_id,
                      debug_oefpc):
    g_lat, g_lon, g_vert, g_time, g_h, g_err = g_arrays
    j_lat, j_lon, j_vert, j_time, j_h, j_err = j_arrays

    if g_lat.size == 0 or j_lat.size == 0:
        return empty_match_result()

    scales = np.asarray(tolerances, dtype=float)
    if np.any(scales <= 0.0):
        raise ValueError("All matching tolerances must be greater than zero")

    # Normalize each coordinate by its own tolerance. This lets a radius-one
    # Chebyshev query represent the requested four-dimensional tolerance box.
    g_pts = np.column_stack((g_lat, g_lon, g_vert, g_time)) / scales
    j_pts = np.column_stack((j_lat, j_lon, j_vert, j_time)) / scales
    tree = cKDTree(g_pts)

    lat_list = []
    lon_list = []
    vert_list = []
    g_list = []
    j_list = []
    dh_list = []
    g_time_list = []
    j_time_list = []
    g_err_list = []
    j_err_list = []

    # Dense NEXRAD volumes can contain several nearby gates. Track used GSI
    # indices so every observation participates in at most one matched pair.
    used_g_indices = set()

    for i, point in enumerate(j_pts):
        nearby = tree.query_ball_point(point, r=1.0, p=np.inf)
        available = [k for k in nearby if k not in used_g_indices]
        if not available:
            continue

        # Select the nearest candidate in normalized coordinate space.
        candidate_indices = np.asarray(available, dtype=int)
        candidate_delta = g_pts[candidate_indices] - point
        distance_squared = np.sum(candidate_delta * candidate_delta, axis=1)
        k = int(candidate_indices[np.argmin(distance_squared)])
        used_g_indices.add(k)

        lat_list.append(j_lat[i])
        lon_list.append(j_lon[i])
        vert_list.append(j_vert[i])
        g_list.append(g_h[k])
        j_list.append(j_h[i])
        dh_list.append(j_h[i] - g_h[k])
        g_time_list.append(g_time[k])
        j_time_list.append(j_time[i])
        g_err_list.append(g_err[k])
        j_err_list.append(j_err[i])

        # -------------------------------------------------------
        # NEW DEBUG PRINT
        # -------------------------------------------------------
        if debug_oefpc:
            delta_err = j_err[i] - g_err[k]

            # print only interesting cases
            if abs(delta_err) > 0.5 or j_err[i] > 10.0:
                print("OEFPC-DEBUG:",
                      f"stnID={station_id or 'NO_STATION_ID'}",
                      f"lat={j_lat[i]:.3f}",
                      f"lon={j_lon[i]:.3f}",
                      f"vert={j_vert[i]:.1f}",
                      f"G_err={g_err[k]:.3f}",
                      f"J_err={j_err[i]:.3f}",
                      f"delta={delta_err:.3f}",
                      f"g_time={g_time[k]:.0f}",
                      f"j_time={j_time[i]:.0f}",
                      sep=" | ")

    return (np.array(lat_list), np.array(lon_list),
            np.array(vert_list), np.array(g_list),
            np.array(j_list), np.array(dh_list),
            np.array(g_time_list), np.array(j_time_list),
            np.array(g_err_list), np.array(j_err_list))


# ---------------------------------------------------------------
# Concatenate matching results from multiple station subsets
# ---------------------------------------------------------------
def concatenate_match_results(parts):
    if not parts:
        return empty_match_result()

    return tuple(
        np.concatenate([part[index] for part in parts])
        for index in range(10)
    )


# ---------------------------------------------------------------
# KD-tree stationID-based matching
# ---------------------------------------------------------------
def match_by_kdtree_station(g, j, lat_tol, lon_tol, vert_tol, time_tol,
                            debug_oefpc):
    total_g = sum([len(g[s]["lat"]) for s in g])
    total_j = sum([len(j[s]["lat"]) for s in j])

    print("")
    print("Using stationID-aware KD-tree matching.")
    print("Total GDIAG assimilated:", total_g)
    print("Total JDIAG assimilated:", total_j)

    tolerances = (lat_tol, lon_tol, vert_tol, time_tol)
    parts = []

    for s in j:
        if s not in g:
            continue

        # Extract GSI and JEDI subsets in the flat-array order expected by the
        # shared matching helper.
        g_arrays = tuple(
            g[s][key] for key in ("lat", "lon", "vert", "time", "hofx", "error")
        )
        j_arrays = tuple(
            j[s][key] for key in ("lat", "lon", "vert", "time", "hofx", "error")
        )
        parts.append(
            match_flat_arrays(
                g_arrays,
                j_arrays,
                tolerances,
                station_id=s,
                debug_oefpc=debug_oefpc,
            )
        )

    result = concatenate_match_results(parts)
    print("Matched obs (KD-tree, stationID aware):", result[0].size)
    return result


# ---------------------------------------------------------------
# KD-tree global tolerance-based matching (no stationID)
# ---------------------------------------------------------------
def match_by_kdtree_global(g, j, lat_tol, lon_tol, vert_tol, time_tol,
                           debug_oefpc):
    print("")
    print("Using global KD-tree tolerance-based matching (no stationID).")

    g_arrays = flatten_obs_dict(g)
    j_arrays = flatten_obs_dict(j)

    print("Total GDIAG assimilated:", g_arrays[0].size)
    print("Total JDIAG assimilated:", j_arrays[0].size)

    if g_arrays[0].size == 0 or j_arrays[0].size == 0:
        print("No observations in one of the inputs.")
        return empty_match_result()

    tolerances = (lat_tol, lon_tol, vert_tol, time_tol)
    result = match_flat_arrays(
        g_arrays,
        j_arrays,
        tolerances,
        station_id=None,
        debug_oefpc=debug_oefpc,
    )

    print("Matched obs (KD-tree, global tolerance):", result[0].size)
    return result


# ---------------------------------------------------------------
# Wrapper: choose station-aware or global matching
# ---------------------------------------------------------------
def match_by_kdtree(g, j, lat_tol, lon_tol, vert_tol, time_tol,
                    use_station_ids, debug_oefpc):
    if use_station_ids:
        return match_by_kdtree_station(
            g, j, lat_tol, lon_tol, vert_tol, time_tol, debug_oefpc
        )

    return match_by_kdtree_global(
        g, j, lat_tol, lon_tol, vert_tol, time_tol, debug_oefpc
    )


# ---------------------------------------------------------------
# Apply pressure-style or height-style limits to a vertical axis
# ---------------------------------------------------------------
def configure_vertical_axis(ax, values, vertical_info):
    ax.set_ylabel(vertical_info["label"])

    if not vertical_info["available"]:
        ax.set_yticks([0.0])
        ax.set_yticklabels(["N/A"])
        return

    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return

    value_min = np.min(finite_values)
    value_max = np.max(finite_values)
    span = value_max - value_min
    pad = 0.05 * span if span > 0.0 else 1.0

    if vertical_info["invert_axis"]:
        ax.set_ylim(value_max + pad, max(0.0, value_min - pad))
    else:
        ax.set_ylim(value_min - pad, value_max + pad)

# ---------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------
parser = argparse.ArgumentParser(description="KD-tree H(x) comparison tool")

parser.add_argument("-g", "--gdiag", required=True)
parser.add_argument("-j", "--jdiag", required=True)
parser.add_argument("--var", required=True)
parser.add_argument("--hofx", default="hofx0")
parser.add_argument("--grid", type=str)

parser.add_argument("--lat_tol", type=float, default=0.0005)
parser.add_argument("--lon_tol", type=float, default=0.0005)
parser.add_argument(
    "--pres_tol",
    type=float,
    default=0.1,
    help="Pressure matching tolerance in hPa.",
)
parser.add_argument(
    "--height_tol",
    type=float,
    default=1.0,
    help="Height matching tolerance in m; used by radialVelocity.",
)
parser.add_argument("--time_tol", type=float, default=10.0)

parser.add_argument("--worst", type=int, default=100)
parser.add_argument("-o", "--output", required=True)

parser.add_argument("--cenlat", type=float, default=37.5)
parser.add_argument("--cenlon", type=float, default=-97.5)
parser.add_argument("--width",  type=float, default=77.0)
parser.add_argument("--height", type=float, default=35.0)
parser.add_argument("--dh_thresh", type=float, default=None, help="Outlier threshold for delta H(x); overrides variable defaults.")
parser.add_argument("--derr_thresh", type=float, default=None, help="Outlier threshold for delta obs error; overrides variable defaults.")
parser.add_argument("--debug_oefpc", action="store_true", help="Print detailed debug info for mismatched obs errors.")


args = parser.parse_args()

print("Loading:", args.gdiag)
g, g_has_stn, g_vertical_info = load_by_station(
    args.gdiag, args.var, args.hofx
)

print("Loading:", args.jdiag)
j, j_has_stn, j_vertical_info = load_by_station(
    args.jdiag, args.var, args.hofx
)

# Matching pressure to height would create invalid pairs. Different height
# metadata names are acceptable because they share the same physical units.
if g_vertical_info["kind"] != j_vertical_info["kind"]:
    raise RuntimeError(
        "GDIAG and JDIAG use incompatible vertical coordinates: "
        f"{g_vertical_info['name']} versus {j_vertical_info['name']}"
    )

vertical_info = j_vertical_info
if vertical_info["kind"] == "height":
    vertical_tol = args.height_tol
else:
    vertical_tol = args.pres_tol

print(
    f"Vertical matching tolerance: {vertical_tol} "
    f"{vertical_info['units']}"
)

use_station_ids = g_has_stn and j_has_stn
if use_station_ids:
    print("Both files have stationIdentification; using station-aware matching.")
else:
    print("At least one file missing stationIdentification; using tolerance-based global matching.")

(lat, lon, vert, g_h, j_h, dh, g_time_match, j_time_match,
 g_err_match, j_err_match) = match_by_kdtree(
    g, j, args.lat_tol, args.lon_tol, vertical_tol, args.time_tol,
    use_station_ids=use_station_ids,
    debug_oefpc=args.debug_oefpc,
)

dh_thresh, d_err_thresh = resolve_thresholds(
    args.var,
    args.dh_thresh,
    args.derr_thresh
)

print(f"Using thresholds: dh_thresh={dh_thresh}, d_err_thresh={d_err_thresh}")
#print("dh stats:", np.nanmin(dh), np.nanmedian(dh), np.nanmax(dh))
#print("abs(dh) p95/p99:", np.nanpercentile(np.abs(dh), [95, 99]))
#print("g_h stats:", np.nanmin(g_h), np.nanmedian(g_h), np.nanmax(g_h))
#print("j_h stats:", np.nanmin(j_h), np.nanmedian(j_h), np.nanmax(j_h))


total_g = sum([len(g[s]["lat"]) for s in g])
total_j = sum([len(j[s]["lat"]) for s in j])

if dh.size == 0:
    print("No matched observations found with current tolerances. Exiting without plots.")
    import sys
    sys.exit(0)

# Convert hofx units for stationPressure (Pa -> hPa)
if args.var == "stationPressure":
    g_h = g_h * 0.01
    j_h = j_h * 0.01
    dh  = dh * 0.01

# ---------------------------------------------------------------
# Domain setup
# ---------------------------------------------------------------
half_w = args.width / 2.0
half_h = args.height / 2.0
left   = args.cenlon - half_w
right  = args.cenlon + half_w
bot    = args.cenlat - half_h
top    = args.cenlat + half_h

domain_ring = None
if args.grid:
    print("Loading domain outline")
    domain_ring = load_domain_ring(args.grid)

proj = ccrs.PlateCarree()

# ---------------------------------------------------------------
# 8-panel figure
# ---------------------------------------------------------------
print("Creating 8-panel figure...")

fig = plt.figure(figsize=(20, 18))

# Radar diagnostics are much denser than conventional diagnostics. Smaller
# markers keep the NEXRAD maps and scatter plots legible and faster to render.
map_point_size = 2 if args.var == "radialVelocity" else 6
scatter_point_size = 1 if args.var == "radialVelocity" else 4

# ----------------------------------------------------------------------
# Panel 1: delta H(x) map
# ----------------------------------------------------------------------
mask_out_dh = np.abs(dh) > dh_thresh

# pick values
lon_out = lon[mask_out_dh]
lat_out = lat[mask_out_dh]
dh_out  = dh[mask_out_dh]

ax1 = fig.add_subplot(4,2,1, projection=proj)
ax1.set_extent([left, right, bot, top], crs=proj)
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

if domain_ring is not None:
    ax1.plot(domain_ring[:,0], domain_ring[:,1], color="gray",
             linewidth=0.8, transform=proj)

vmax = 1

s1 = ax1.scatter(lon, lat, c=dh, s=map_point_size,
                 cmap="coolwarm", vmin=-vmax, vmax=vmax,
                 transform=proj)
s1 = ax1.scatter(lon_out, lat_out, c=dh_out, s=map_point_size,
                 cmap="coolwarm", vmin=-vmax, vmax=vmax,
                 transform=proj)
plt.colorbar(s1, ax=ax1, shrink=0.7, label="delta H(x)")
ax1.set_title("delta H(x) = JEDI - GSI")

# ----------------------------------------------------------------------
# Panel 2: delta obs error map
# ----------------------------------------------------------------------
d_err = j_err_match - g_err_match
mask_out_err = np.abs(d_err) > d_err_thresh

# Outliers
lon_out = lon[mask_out_err]
lat_out = lat[mask_out_err]
d_err_out = d_err[mask_out_err]

ax2 = fig.add_subplot(4, 2, 2, projection=proj)
ax2.set_extent([left, right, bot, top], crs=proj)
ax2.add_feature(cfeature.COASTLINE)
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)

if domain_ring is not None:
    ax2.plot(domain_ring[:,0], domain_ring[:,1], color="gray",
             linewidth=0.8, transform=proj)

vmax_err=1

s2 = ax2.scatter(lon, lat, c=d_err, s=map_point_size,
                 cmap="coolwarm", vmin=-vmax_err, vmax=vmax_err,
                 transform=proj)
s2 = ax2.scatter(lon_out, lat_out, c=d_err_out, s=map_point_size,
                 cmap="coolwarm", vmin=-vmax_err, vmax=vmax_err,
                 transform=proj)
plt.colorbar(s2, ax=ax2, shrink=0.7, label="delta obs error")
ax2.set_title("delta obs error map (JEDI - GSI)")

# ----------------------------------------------------------------------
# Panel 3: vertical-coordinate-colored map
# ----------------------------------------------------------------------
ax3 = fig.add_subplot(4,2,3, projection=proj)
ax3.set_extent([left, right, bot, top], crs=proj)
ax3.add_feature(cfeature.COASTLINE)
ax3.add_feature(cfeature.BORDERS)
ax3.add_feature(cfeature.STATES)

if domain_ring is not None:
    ax3.plot(domain_ring[:,0], domain_ring[:,1], color="gray",
             linewidth=0.8, transform=proj)

if vertical_info["available"]:
    finite_vert = vert[np.isfinite(vert)]
    if vertical_info["kind"] == "pressure":
        color_min = 100.0
        color_max = 1050.0
    elif finite_vert.size:
        color_min = np.min(finite_vert)
        color_max = np.max(finite_vert)
        if color_min == color_max:
            color_min -= 0.5
            color_max += 0.5
    else:
        color_min = 0.0
        color_max = 1.0

    s3 = ax3.scatter(lon, lat, c=vert, s=map_point_size,
                     cmap="viridis", vmin=color_min, vmax=color_max,
                     transform=proj)
    plt.colorbar(
        s3, ax=ax3, shrink=0.7, label=vertical_info["label"]
    )
    ax3.set_title(f"{vertical_info['label']}-colored scatter")
else:
    # Observation types without a vertical coordinate still receive a useful
    # location map without a misleading pressure color scale.
    ax3.scatter(lon, lat, s=map_point_size, c="gray", transform=proj)
    ax3.set_title("Observation locations (no vertical coordinate)")

# ----------------------------------------------------------------------
# Panel 4: time hist
# ----------------------------------------------------------------------
ax4 = fig.add_subplot(4, 2, 4)
g_off = []
j_off = []

for s in g:
    g_off.extend(list(g[s]["timeOffset"]))
for s in j:
    j_off.extend(list(j[s]["timeOffset"]))

g_off = np.array(g_off, dtype=float)
j_off = np.array(j_off, dtype=float)

# Keep only finite values before plotting.
g_off = g_off[np.isfinite(g_off)]
j_off = j_off[np.isfinite(j_off)]

bins = np.linspace(-5400, 5400, 61)

ax4.hist(j_off, bins=bins, alpha=0.8, label="JDIAG timeOffset", color="orange")
ax4.hist(g_off, bins=bins, alpha=0.5, label="GDIAG timeOffset", color="blue", zorder=10)

ax4.axvline(-5400, color="k", linestyle="--", label="twindow min (-5400)")
ax4.axvline( 5400, color="k", linestyle="--", label="twindow max (+5400)")

ax4.set_xlim(-6000, 6000)
ax4.set_xlabel("timeOffset (seconds)")
ax4.set_ylabel("Count")
ax4.set_title("Time-distribution comparison (timeOffset): GDIAG vs JDIAG")
ax4.legend()

# ----------------------------------------------------------------------
# Panel 5: delta H(x) vs vertical coordinate
# ----------------------------------------------------------------------
ax5 = fig.add_subplot(4,2,5)
ax5.scatter(dh, vert, s=scatter_point_size, c="k")
ax5.set_xlabel("delta H(x)")
configure_vertical_axis(ax5, vert, vertical_info)
ax5.set_title(f"delta H(x) vs {vertical_info['label']}")
ax5.set_xlim(-4, 4)

# ----------------------------------------------------------------------
# Panel 6: delta obs error vs vertical coordinate
# ----------------------------------------------------------------------
ax6 = fig.add_subplot(4,2,6)
d_err = j_err_match - g_err_match
#print("d_err stats:", np.nanmin(d_err), np.nanmedian(d_err), np.nanmax(d_err))
#print("g_err stats:", np.nanmin(g_err_match), np.nanmedian(g_err_match), np.nanmax(g_err_match))
#print("j_err stats:", np.nanmin(j_err_match), np.nanmedian(j_err_match), np.nanmax(j_err_match))
#print("len(d_err):", len(d_err))
#print("len(vert):", len(vert))
#print("outliers:", np.sum(np.abs(d_err) > d_err_thresh))

n_exceed = np.sum(np.abs(d_err) > d_err_thresh)
pct_exceed = 100.0 * n_exceed / d_err.size if d_err.size > 0 else 0.0
ax6.scatter(d_err, vert, s=scatter_point_size, c="k")
#mask_in  = np.abs(d_err) <= d_err_thresh
#mask_out = np.abs(d_err) >  d_err_thresh
#
#ax6.scatter(d_err[mask_in],  vert[mask_in],  s=4, c="k", alpha=0.4)
#ax6.scatter(d_err[mask_out], vert[mask_out], s=6, c="red")

ax6.set_xlabel("delta obs error (JEDI - GSI)")
configure_vertical_axis(ax6, vert, vertical_info)
ax6.set_title(f"delta obs error vs {vertical_info['label']}")
# Compute xlim safely
if d_err.size == 0:
    xlim = 1.0
else:
    xlim = np.nanmax(np.abs(d_err))  # <- SAFE even when all NaN
    if not np.isfinite(xlim) or xlim == 0:
        xlim = 1.0
    if xlim >= 500.0:
        xlim = 10.0

ax6.set_xlim(-1*xlim, xlim)
# Add info text similar to Panel 6
info6 = (
    f"delta error threshold: {d_err_thresh}\n"
    f"Count exceeding: {n_exceed}\n"
    f"Percent: {pct_exceed:.1f}%\n"
    f"N matched: {d_err.size}"
)

ax6.text(
    0.02,
    0.98,
    info6,
    transform=ax6.transAxes,
    ha="left",
    va="top",
    fontsize=9,
    family="monospace",
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)

# ----------------------------------------------------------------------
# Panel 7: 1-to-1 comparison
# ----------------------------------------------------------------------
ax7 = fig.add_subplot(4,2,7)
ax7.scatter(g_h, j_h, s=scatter_point_size, c="gray")

# Only finite values can define the identity and regression lines.
mask = np.isfinite(g_h) & np.isfinite(j_h)
g_fit = g_h[mask]
j_fit = j_h[mask]
if g_fit.size:
    lo = min(g_fit.min(), j_fit.min())
    hi = max(g_fit.max(), j_fit.max())
    if lo == hi:
        lo -= 0.5
        hi += 0.5
else:
    lo = -1.0
    hi = 1.0

ax7.plot([lo, hi], [lo, hi], "k--")
ax7.set_xlabel("GSI hofx")
ax7.set_ylabel("JEDI hofx")
ax7.set_title("JEDI/GSI H(x) comparison")
# Annotate assimilated counts
match_mode = "stationID-aware" if use_station_ids else "global-tolerance"

percent_g = 100.0 * dh.size / total_g if total_g > 0 else 0.0
percent_j = 100.0 * dh.size / total_j if total_j > 0 else 0.0

# --- Linear trendline (least squares) ---
if g_fit.size > 1:
    # Fit slope and intercept: j = a*g + b
    a, b = np.polyfit(g_fit, j_fit, 1)

    # Compute R^2
    j_pred = a * g_fit + b
    ss_res = np.sum((j_fit - j_pred)**2)
    ss_tot = np.sum((j_fit - np.mean(j_fit))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan

    # Plot trendline over same range as 1:1 line
    xx = np.array([lo, hi])
    yy = a * xx + b
    ax7.plot(xx, yy, "r-", linewidth=1.2, label="trendline")

    # Add equation text
    eqn = (
        f"trend: J = {a:.3f} * G + {b:.3f}\n"
        f"R^2 = {r2:.3f}"
    )

    ax7.text(
        0.70,
        0.02,
        eqn,
        transform=ax7.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        family="monospace",
        bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
    )


info = (
    f"GDIAG assimilated: {total_g}\n"
    f"JDIAG assimilated: {total_j}\n"
    f"Matched pairs: {dh.size}\n"
    f"Match % (GDIAG): {percent_g:.1f}%\n"
    f"Match % (JDIAG): {percent_j:.1f}%\n"
    f"Mode: {match_mode}\n"
    f"lat_tol: {args.lat_tol}\n"
    f"lon_tol: {args.lon_tol}\n"
    f"{vertical_info['kind']}_tol: {vertical_tol} "
    f"{vertical_info['units']}\n"
    f"time_tol: {args.time_tol}"
)

ax7.text(
    0.02,
    0.98,
    info,
    transform=ax7.transAxes,
    ha="left",
    va="top",
    fontsize=9,
    family="monospace",
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)

# ----------------------------------------------------------------------
# Panel 8: Obs error one-to-one comparison (with inlier vs outlier analysis)
# ----------------------------------------------------------------------
ax8 = fig.add_subplot(4, 2, 8)

# Raw scatter
ax8.scatter(
    g_err_match,
    j_err_match,
    s=scatter_point_size,
    c="gray",
    label="all obs",
    alpha=0.5,
)

# --------------- Define inliers vs outliers -----------------
threshold = d_err_thresh
d_err = j_err_match - g_err_match
sentinel = 1.0e10  # or detect from metadata
mask_phys = (g_err_match < sentinel) & (j_err_match < sentinel)
mask_all = np.isfinite(g_err_match) & np.isfinite(j_err_match) & mask_phys

#mask_all = np.isfinite(g_err_match) & np.isfinite(j_err_match)
mask_in  = mask_all & (np.abs(d_err) <= threshold)
mask_out = mask_all & (np.abs(d_err) > threshold)

g_in,  j_in  = g_err_match[mask_in],  j_err_match[mask_in]
g_out, j_out = g_err_match[mask_out], j_err_match[mask_out]

# Highlight outliers
ax8.scatter(g_out, j_out, s=2, c="red", label=f"outliers ({g_out.size})", alpha=0.8)

# Axis limits
valid_error_values = np.concatenate([
    g_err_match[mask_all],
    j_err_match[mask_all],
])
if valid_error_values.size:
    absmax = np.max(np.abs(valid_error_values))
else:
    absmax = 1.0
if absmax == 0 or not np.isfinite(absmax):
    absmax = 1.0
if absmax >= 500.0:
    absmax = 10.0
pad = absmax * 0.05
ax8.set_xlim(-(absmax + pad), absmax + pad)
ax8.set_ylim(-(absmax + pad), absmax + pad)

# --------------- Regression helper -----------------
def compute_regression(x, y):
    if x.size < 2:
        return np.nan, np.nan, np.nan
    a, b = np.polyfit(x, y, 1)
    y_pred = a * x + b
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
    return a, b, r2

# Compute regressions
a_in,  b_in,  r2_in  = compute_regression(g_in,  j_in)
a_all, b_all, r2_all = compute_regression(g_err_match[mask_all], j_err_match[mask_all])

# --------------- Plot regression lines ---------------
xx = np.array([-(absmax + pad), absmax + pad])

# Inlier trendline (BLUE)
if np.isfinite(a_in):
    ax8.plot(xx, a_in*xx + b_in, "b-", linewidth=1.2, label="inlier trend")

# All-obs trendline (RED)
if np.isfinite(a_all):
    ax8.plot(xx, a_all*xx + b_all, "r--", linewidth=1.2, label="all-points trend")

# Identity line
ax8.plot(xx, xx, "k-", linewidth=1.0, label="1:1 line")

# Labels
ax8.set_xlabel("GSI obs error")
ax8.set_ylabel("JEDI obs error")
ax8.set_title("JEDI/GSI obs error comparison")

# Info box
info = (
    f"Total matched obs: {mask_all.sum()}\n"
    f"Inliers (|$\\Delta$ error| <= {threshold}): {g_in.size}\n"
    f"Outliers: {g_out.size}\n"
    f"\n"
    f"Inlier trend: J = {a_in:.3f}*G + {b_in:.3f}\n"
    f"Inlier R$^2$ = {r2_in:.3f}\n"
    f"\n"
    f"All-points trend: J = {a_all:.3f}*G + {b_all:.3f}\n"
    f"All-points R$^2$ = {r2_all:.3f}"
)

ax8.text(
    0.02, 0.02, info,
    transform=ax8.transAxes,
    ha="left", va="bottom",
    fontsize=8, family="monospace",
    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
)

ax8.legend(loc="upper left", fontsize=8)

# ---------------------------------------------------------------
# Simple supertitle
# ---------------------------------------------------------------

# Extract from filename: gdiag_adpupa_airTemperature_120.nc
fname = args.gdiag.split("/")[-1]
short = fname.replace("gdiag_", "").replace(".nc", "")

fig.suptitle(short, fontsize=16, fontweight="bold", y=0.98)

plt.tight_layout()
plt.savefig(args.output + "_8panel.png", dpi=180)
print("Saved:", args.output + "_8panel.png")

print("Done.")

