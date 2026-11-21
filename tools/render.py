#!/usr/bin/env python3
"""render.py -- PNG panels + contact sheet from the giants' VTK slices (P1 of the 2026-11-20 viz proposal).

Usage:
  python3 -I tools/render.py RUN [RUN_B] [--out DIR] [--kinds zonal,radial] [--ranges FILE] [--every N] [--jobs N]
                              [--levels I,I,...] [--zoom LAT0,LAT1,LON0,LON1]

Kinds: zonal, radial, panorama, proczm, cloudview (default all that are present).
cloudview (2026-11-21): <Planet>_cloudview_<it>.vtk -- every column array (X_above, the column above a layer) as a lat x lon map
  for layer 0 (= the whole column) and for each --levels layer, plus PhotoRGB (the writer's colour recipe, per-file scale) as an image.

--levels: panorama levels (i, 0 = bottom) for the wind / temperature maps (zonal wind w, meridional wind v,
  the departure from the zonal mean -- speed with arrows -- and T minus its zonal mean); default the middle level.
--zoom: crop the lat x lon maps (radial slices, panorama) to the window, e.g. --zoom=-40,0,150,210 for the GRS (with "=", a leading minus would read as an option)
  at 22 S 180 E; the colour ranges are then taken inside the window. Output goes to png_zoom/ by default.

RUN is an arm directory (the one holding output*/) or the output directory itself.
With RUN_B the same panels are drawn for B and as A-B differences at matching iterations.

Output (default RUN/png, or RUN/png_vs_<B> for pairs):
  <kind>/<panel>_<iter>.png    one PNG per panel and iteration (A_, B_, D_ prefixes in pair mode)
  ranges.json                  the colour ranges used -- fixed over ALL iterations and BOTH runs;
                               pass it back with --ranges to keep later renders comparable
  index.html                   contact sheet, rows = iterations, columns = panels

Conventions taken from the writers (ParaView_<P>.cpp, ParaViewWriter.h open_slice):
  zonal slice: DIMENSIONS jm im 1, point index i*jm + j, j = 0 at 90 N; 'height' field = level height in km.
  radial slice: DIMENSIONS km jm 1, point index j*km + k.
  Temperature: Saturn writes t*t_ref/10 (= K/10); Jupiter / Neptune / Uranus write t*t_ref - 273.15 (C).
  Species: x 1e3 of kg/m3 = g/m3; NH4SH (and w_/j_ nh4sh) x 1e9 = ug/m3, shown here in g/m3.
  Precipitation P_*: written as the model's arrays (kg/m2/s); Q_*: W/m3.
"""
import argparse, glob, json, os, re, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np

# ------------------------------------------------------------------ reading
def read_slice(path):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n
    r = vtk.vtkStructuredGridReader(); r.SetFileName(path)
    r.ReadAllScalarsOn(); r.ReadAllVectorsOn(); r.Update()
    g = r.GetOutput(); nf, ns, _ = g.GetDimensions(); pd = g.GetPointData(); d = {}
    for a in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(a); x = v2n(arr)
        if x.ndim == 1 and arr.GetName() not in d:     # rho_mix is written twice; keep the first
            d[arr.GetName()] = x.reshape(ns, nf).astype(np.float64)
    return d

PLANET = {"Saturn": "sat", "Jupiter": "jup", "Neptune": "nept", "Uranus": "uran"}

def temp_K(planet, x):
    return x * 10.0 if planet == "Saturn" else x + 273.15

# ------------------------------------------------------------------ panels
_KIND = ["zonal"]          # kind of the slice being built (set in main), for builders that depend on it

def t_anomaly(t):
    """T minus the level mean. Zonal slices (rows = levels, cols = latitudes, j = 0 at 90 N): the
    AREA mean of the level, weighted by cos(latitude) as TB_DIAG weights it (sin(colatitude)), so the
    zero is the planet's average at that height (2026-11-20; a plain mean over the rows before).
    Radial slices (rows = latitudes, cols = longitudes): the mean along the latitude circle, unweighted."""
    if _KIND[0] in ("zonal", "proczm"):
        w = np.cos(np.radians(90.0 - np.arange(t.shape[1]) * 180.0 / (t.shape[1] - 1)))
        return t - (t * w).sum(axis=1, keepdims=True) / w.sum()
    return t - t.mean(axis=1, keepdims=True)

# name, title, unit, builder(d, planet) -> 2D array or None, scale ('lin' | 'sym' | 'log'), overlay contours
def g(d, *names, f=1.0):
    if not all(n in d for n in names): return None
    return sum(d[n] for n in names) * f

PANELS = [
  ("T",        "Temperature",                 "K",      lambda d, p: None if "Temperature" not in d else temp_K(p, d["Temperature"]), "lin", None),
  ("Tanom",    "T minus level mean",          "K",      lambda d, p: None if "Temperature" not in d else t_anomaly(temp_K(p, d["Temperature"])), "sym", None),
  ("v",        "v (meridional)",              "m/s",    lambda d, p: g(d, "v-Component"), "sym", None),
  ("w",        "w (zonal)",                   "m/s",    lambda d, p: g(d, "w-Component"), "sym", None),
  ("u",        "u (radial)",                  "m/s",    lambda d, p: g(d, "u-Component"), "sym", None),
  ("H2O",      "H2O vapour; cloud / ice",     "g/m3",   lambda d, p: g(d, "H2O"), "log", ("H2OCloud", "H2OIce")),
  ("NH3",      "NH3 vapour; cloud / ice",     "g/m3",   lambda d, p: g(d, "NH3"), "log", ("NH3Cloud", "NH3Ice")),
  ("H2S",      "H2S vapour; cloud / ice",     "g/m3",   lambda d, p: g(d, "H2S"), "log", ("H2SCloud", "H2SIce")),
  ("NH4SH",    "NH4SH",                       "g/m3",   lambda d, p: g(d, "NH4SH", f=1e-6), "log", None),
  ("CH4cond",  "CH4 cloud + ice",             "g/m3",   lambda d, p: g(d, "CH4Cloud", "CH4Ice"), "log", None),
  ("H2Ocond",  "H2O cloud + ice",             "g/m3",   lambda d, p: g(d, "H2OCloud", "H2OIce"), "log", None),
  ("NH3cond",  "NH3 cloud + ice",             "g/m3",   lambda d, p: g(d, "NH3Cloud", "NH3Ice"), "log", None),
  ("Chrom",    "Chromophore",                 "mg/m3",  lambda d, p: g(d, "Chromophore_mgm3"), "log", None),
  ("QLat",     "Latent heating Q_Latent",     "W/m3",   lambda d, p: g(d, "Q_Latent"), "sym", None),
  ("Qprecip",  "Precip latent heat Q_precip", "W/m3",   lambda d, p: g(d, "Q_precip"), "sym", None),
  ("Qrad",     "Radiative heating Q_rad",     "W/m3",   lambda d, p: g(d, "Q_rad"), "sym", None),
  ("Ph2o",     "H2O precip flux (rain+snow+graupel)", "kg/m2/s", lambda d, p: g(d, "P_rain", "P_snow", "P_graupel"), "log", None),
  ("Pnh3",     "NH3 precip flux (rain+snow)", "kg/m2/s",lambda d, p: g(d, "P_nh3_rain", "P_nh3_snow"), "log", None),
  ("Pch4",     "CH4 precip flux (rain+snow)", "kg/m2/s",lambda d, p: g(d, "P_ch4_rain", "P_ch4_snow"), "log", None),
  ("Pnh4sh",   "NH4SH settling flux",         "as written", lambda d, p: g(d, "P_nh4sh"), "log", None),
  ("Buoy",     "Buoyancy force",              "N/m3",   lambda d, p: g(d, "BuoyancyForce"), "sym", None),
  ("nuet",     "Eddy viscosity nu_t",         "m2/s",   lambda d, p: g(d, "nue_t"), "log", None),
]
# process fields (ATSAT_PROC_DIAG, ProcSat.h): zonal + time means written by the model as <Planet>_proczm_<it>.vtk (until 2026-11-20 <Planet>_proczm_0_<it>.vtk; both are read)
def _pg(*n):
    return lambda d, p: g(d, *n)
PROC_PANELS = [
  ("p_dT_satadj",  "SatAdj heating",                       "K/s",     _pg("dT_satadj"), "sym", None),
  ("p_cond_H2O",   "H2O condensation (+) / evaporation (-)", "g/m3/s", _pg("cond_H2O", "evap_H2O"), "sym", None),
  ("p_cond_NH3",   "NH3 condensation (+) / evaporation (-)", "g/m3/s", _pg("cond_NH3", "evap_NH3"), "sym", None),
  ("p_cond_CH4",   "CH4 condensation (+) / evaporation (-)", "g/m3/s", _pg("cond_CH4", "evap_CH4"), "sym", None),
  ("p_cond_H2S",   "H2S condensation (+) / evaporation (-)", "g/m3/s", _pg("cond_H2S", "evap_H2S"), "sym", None),
  ("p_nh4sh_form", "NH4SH formation",                      "g/m3/s",  _pg("nh4sh_form"), "lin", None),
  ("p_nh4sh_dec",  "NH4SH decomposition",                  "g/m3/s",  lambda d, p: None if "nh4sh_decomp" not in d else -d["nh4sh_decomp"], "lin", None),
  ("p_nh4sh_net",  "NH4SH net formation",                  "g/m3/s",  _pg("nh4sh_form", "nh4sh_decomp"), "sym", None),
  ("p_dT_chem",    "NH4SH heating",                        "K/s",     _pg("dT_chem"), "sym", None),
  ("p_conv_H2O",   "H2O condensate -> precip (conv)",      "g/m3/s",  _pg("conv_H2O"), "sym", None),
  ("p_conv_NH3",   "NH3 condensate -> precip (conv)",      "g/m3/s",  _pg("conv_NH3"), "sym", None),
  ("p_revap_H2O",  "H2O re-evaporation of precip",         "g/m3/s",  _pg("revap_H2O"), "sym", None),
  ("p_revap_NH3",  "NH3 re-evaporation of precip",         "g/m3/s",  _pg("revap_NH3"), "sym", None),
  ("p_conv_H2S",   "H2S condensate -> precip (conv)",      "g/m3/s",  _pg("conv_H2S"), "sym", None),
  ("p_revap_H2S",  "H2S re-evaporation of precip",         "g/m3/s",  _pg("revap_H2S"), "sym", None),
  ("p_conv_CH4",   "CH4 condensate -> precip (conv)",      "g/m3/s",  _pg("conv_CH4"), "sym", None),
  ("p_dT_precip",  "Precip heating",                       "K/s",     _pg("dT_precip"), "sym", None),
  ("p_Pflux_H2O",  "H2O precip flux (mean)",               "kg/m2/s", _pg("Pflux_H2O"), "lin", None),
  ("p_Pflux_NH3",  "NH3 precip flux (mean)",               "kg/m2/s", _pg("Pflux_NH3"), "lin", None),
  ("p_Pflux_H2S",  "H2S precip flux (mean)",               "kg/m2/s", _pg("Pflux_H2S"), "lin", None),
  ("p_Pflux_CH4",  "CH4 precip flux (mean)",               "kg/m2/s", _pg("Pflux_CH4"), "lin", None),
  ("p_Pflux_NH4SH","NH4SH settling flux (mean)",           "kg/m2/s", _pg("Pflux_NH4SH"), "lin", None),
  ("p_mca_freq",   "MCA firing frequency",                 "share of calls", _pg("mca_freq"), "lin", None),
  ("p_dT_mca",     "MCA heating",                          "K/s",     _pg("dT_mca"), "sym", None),
  ("p_dq_H2O_mca", "MCA H2O redistribution",               "g/m3/s",  _pg("dq_H2O_mca"), "sym", None),
  ("p_dq_NH3_mca", "MCA NH3 redistribution",               "g/m3/s",  _pg("dq_NH3_mca"), "sym", None),
  ("p_dq_H2S_mca", "MCA H2S redistribution",               "g/m3/s",  _pg("dq_H2S_mca"), "sym", None),
  ("p_dq_CH4_mca", "MCA CH4 redistribution",               "g/m3/s",  _pg("dq_CH4_mca"), "sym", None),
  ("p_dNH4SH_mca", "MCA NH4SH redistribution",             "g/m3/s",  _pg("dNH4SH_mca"), "sym", None),
]
PANELS = PANELS + PROC_PANELS
PROC_KEYS = {p[0] for p in PROC_PANELS}
# radial slices (lat x lon at one level): a subset that reads well as a map
RADIAL_KEYS = ["T", "Tanom", "v", "w", "u", "H2Ocond", "NH3cond", "NH4SH", "QLat", "Ph2o", "Pnh3"]

# ------------------------------------------------------------------ files
def _has(o): return any(glob.glob(os.path.join(o, p)) for p in ("*_zonal_*.vtk", "*_radial_*.vtk", "*_panorama_*.vts", "*_proczm_*.vtk", "*_cloudview_*.vtk"))
def find_output(run):
    if _has(run): return run
    for o in sorted(glob.glob(os.path.join(run, "output*"))):
        if os.path.isdir(o) and _has(o): return o
    sys.exit(f"no *_zonal_*.vtk / *_radial_*.vtk / *_panorama_*.vts under {run}")

PAT = re.compile(r"(?P<planet>[A-Za-z]+)_(?P<kind>zonal|radial|longal|proczm)_(?P<pos>\d+)_(?P<it>\d+)\.vtk$")
PANO = re.compile(r"(?P<planet>[A-Za-z]+)_panorama_(?P<it>\d+)\.vts$")
PROCZ = re.compile(r"(?P<planet>[A-Za-z]+)_proczm_(?P<it>\d+)\.vtk$")          # <Planet>_proczm_<it>.vtk (since 2026-11-20)
CLOUDV = re.compile(r"(?P<planet>[A-Za-z]+)_cloudview_(?P<it>\d+)\.vtk$")      # <Planet>_cloudview_<it>.vtk
def list_files(outdir, kinds):
    fs = {}
    for f in glob.glob(os.path.join(outdir, "*.vtk")):
        m = PAT.search(os.path.basename(f))
        if m and m["kind"] in kinds: fs[(m["kind"], int(m["pos"]), int(m["it"]))] = (f, m["planet"]); continue
        m = PROCZ.search(os.path.basename(f))
        if m and "proczm" in kinds: fs[("proczm", 0, int(m["it"]))] = (f, m["planet"]); continue
        m = CLOUDV.search(os.path.basename(f))
        if m and "cloudview" in kinds: fs[("cloudview", 0, int(m["it"]))] = (f, m["planet"])
    if "panorama" in kinds:
        for f in glob.glob(os.path.join(outdir, "*_panorama_*.vts")):
            m = PANO.search(os.path.basename(f))
            if m: fs[("panorama", 0, int(m["it"]))] = (f, m["planet"])
    return fs

# ------------------------------------------------------------------ panorama -> lat x lon maps
PANO_NEED = ["Temperature", "H2OCloud", "H2OIce", "NH3Cloud", "NH3Ice", "CH4Cloud", "CH4Ice", "NH4SH", "Chromophore_mgm3",
             "P_rain_mmd", "P_snow_mmd", "P_graupel_mmd", "P_nh3_rain_mmd", "P_nh3_snow_mmd", "P_ch4_rain_mmd", "P_ch4_snow_mmd",
             "P_nh4sh_mmd", "Radiation", "Q_Latent", "u-component", "v-component", "w-component"]
LEVELS = []                # panorama levels for the wind / T maps (set in main)
ZOOM = None                # (lat0, lat1, lon0, lon1) or None (set in main)
QUIVER = {}                # map name -> function(d) -> (U, V) arrows, for the eddy-wind maps

def zoom_idx(nj, nk, lat=None):
    """row / column indices of the zoom window on a lat x lon map (j = 0 at 90 N, k = 0 at 0 E; lat given for the cloudview)."""
    if lat is None: lat = 90.0 - np.arange(nj) * 180.0 / (nj - 1)
    lon = np.arange(nk) * 360.0 / (nk - 1)
    la0, la1, lo0, lo1 = ZOOM
    return (np.where((lat >= min(la0, la1)) & (lat <= max(la0, la1)))[0],
            np.where((lon >= min(lo0, lo1)) & (lon <= max(lo0, lo1)))[0])
def crop(Z, lat=None):
    if ZOOM is None or Z is None: return Z
    jj, kk = zoom_idx(Z.shape[0], Z.shape[1], lat); return Z[np.ix_(jj, kk)]
def read_pano(path):
    """{name: array (j, k, i)} -- the panorama is i-fastest, extent (im, jm, km); j = 0 at 90 N."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n
    r = vtk.vtkXMLStructuredGridReader(); r.SetFileName(path)
    r.Update(); g = r.GetOutput(); im, jm, km = g.GetDimensions(); pd = g.GetPointData(); d = {}
    for n in PANO_NEED:
        a = pd.GetArray(n)
        if a is not None: d[n] = v2n(a).reshape(km, jm, im).transpose(1, 0, 2).astype(np.float64)
    return d

def read_cloudview(path):
    """<Planet>_cloudview_<it>.vtk -> ({name: array (i, j, k)}, rgb (i, j, k, 3) or None, lat [deg] per row).
    Point order i, j, k (k fastest), DIMENSIONS km nj im; the spherical default leaves the pole rows out, so the
    latitude is taken from the points (shells: asin(z / |p|); flat stack: x = 0.1 j over the full 90 N..90 S)."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n
    r = vtk.vtkStructuredGridReader(); r.SetFileName(path)
    r.ReadAllScalarsOn(); r.ReadAllColorScalarsOn(); r.Update()
    g = r.GetOutput(); km, nj, im = g.GetDimensions(); pd = g.GetPointData(); d = {}; rgb = None
    for a in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(a); x = v2n(arr)
        if arr.GetName() == "PhotoRGB":
            x = x.astype(np.float64); rgb = (x / 255.0 if x.max() > 1.0 else x).reshape(im, nj, km, 3)
        elif x.ndim == 1: d[arr.GetName()] = x.reshape(im, nj, km).astype(np.float64)
    pts = v2n(g.GetPoints().GetData()).reshape(im, nj, km, 3)[0, :, 0, :]
    rr = np.linalg.norm(pts, axis=1)
    if rr.min() > 0.5 and np.ptp(rr) < 1e-3: lat = np.degrees(np.arcsin(np.clip(pts[:, 2] / rr, -1, 1)))
    else: lat = 90.0 - np.round(pts[:, 0] / 0.1) * 180.0 / (nj - 1)
    return d, rgb, lat

def cv_maps(d, rgb, hkm):
    """cloudview lat x lon maps: (name, title, unit, array, scale); layer 0 and the LEVELS layers."""
    M = []; im = next(iter(d.values())).shape[0]
    for L in [0] + [l for l in LEVELS if l != 0]:
        if not 0 <= L < im - 1: continue
        hl = "whole column (layer 0)" if L == 0 else f"above layer {L}" + (f" ({hkm[L]:.0f} km)" if hkm is not None else "")
        sfx = "" if L == 0 else f"_L{L}"
        if rgb is not None: M.append((f"cvPhoto{sfx}", "PhotoRGB (writer's recipe, per-file scale), " + hl, "", rgb[L], "rgb"))
        for n, x in d.items():
            unit = n.rsplit("_", 1)[-1].replace("gm2", "g/m2").replace("mgm2", "mg/m2") if "_" in n else "as written"
            unit = "mg/m2" if n.endswith("mgm2") else unit
            M.append((f"cv{n.split('_')[0]}{sfx}", f"{n}, " + hl, unit, x[L], "lin"))
    return M

def level_heights_km(outdir, im):
    """level heights from any zonal slice in the same output directory ('height' field), else None."""
    for f in sorted(glob.glob(os.path.join(outdir, "*_zonal_*.vtk"))):
        z = read_slice(f)
        if "height" in z and z["height"].shape[0] == im: return z["height"][:, 0]
    # fallback: uniform levels over L_atm [km] from the arm's config.xml (exact on ATSAT, whose levels are uniform)
    for c in (os.path.join(outdir, "config.xml"), os.path.join(os.path.dirname(os.path.normpath(outdir)), "config.xml")):
        if os.path.isfile(c):
            m = re.search(r"<L_atm>\s*([0-9.eE+-]+)\s*</L_atm>", open(c).read())
            if m:
                print(f"  level heights: uniform over L_atm = {m[1]} km from {c} (no zonal slice)"); return np.linspace(0, float(m[1]), im)
    return None

def pano_maps(d, hkm, planet=""):
    """lat x lon maps: (name, title, unit, array, scale)."""
    dz = np.gradient(hkm) * 1e3 if hkm is not None else None      # m per level
    def col(*ns, f=1.0):
        if not all(n in d for n in ns): return None
        x = sum(d[n] for n in ns) * f
        return (x * dz).sum(axis=2) if dz is not None else x.max(axis=2)
    cu = "g/m2" if dz is not None else "g/m3 (column max, no heights)"
    def top(*ns):
        if hkm is None or not all(n in d for n in ns): return None
        x = sum(d[n] for n in ns); thr = 1e-3 * max(np.nanmax(x), 1e-30)
        idx = np.where(x > thr, np.arange(x.shape[2])[None, None, :], -1).max(axis=2)
        return np.where(idx >= 0, hkm[np.clip(idx, 0, None)], np.nan)
    def pmax(*ns):
        if not all(n in d for n in ns): return None
        return sum(d[n] for n in ns).max(axis=2)
    M = [("colH2O",   "H2O cloud + ice column",     cu,  col("H2OCloud", "H2OIce"), "lin"),
         ("colNH3",   "NH3 cloud + ice column",     cu,  col("NH3Cloud", "NH3Ice"), "lin"),
         ("colCH4",   "CH4 cloud + ice column",     cu,  col("CH4Cloud", "CH4Ice"), "lin"),
         ("colNH4SH", "NH4SH column",               cu,  col("NH4SH", f=1e-6), "lin"),
         ("colChrom", "Chromophore column",         "mg/m2" if dz is not None else "mg/m3 (max)", col("Chromophore_mgm3"), "lin"),
         ("topH2O",   "H2O condensate top (1e-3 of max)", "km", top("H2OCloud", "H2OIce"), "lin"),
         ("topNH3",   "NH3 condensate top (1e-3 of max)", "km", top("NH3Cloud", "NH3Ice"), "lin"),
         ("PH2O",     "H2O precip, column max",      "mm/d", pmax("P_rain_mmd", "P_snow_mmd", "P_graupel_mmd"), "lin"),
         ("PNH3",     "NH3 precip, column max",      "mm/d", pmax("P_nh3_rain_mmd", "P_nh3_snow_mmd"), "lin"),
         ("PCH4",     "CH4 precip, column max",      "mm/d", pmax("P_ch4_rain_mmd", "P_ch4_snow_mmd"), "lin"),
         ("PNH4SH",   "NH4SH settling, column max",  "mm/d", pmax("P_nh4sh_mmd"), "lin"),
         ("QLatcol",  "Latent heating, column max",  "W/m3", None if "Q_Latent" not in d else d["Q_Latent"].max(axis=2), "lin"),
         ("Rlid",     "Radiation at the lid (log's OLR quantity)", "as written", None if "Radiation" not in d else d["Radiation"][:, :, -1], "lin")]
    im = next(iter(d.values())).shape[2]
    for L in LEVELS:
        if not 0 <= L < im: continue
        hl = f" level {L}" + (f" ({hkm[L]:.0f} km)" if hkm is not None else "")
        if "w-component" in d and "v-component" in d:
            w = d["w-component"][:, :, L]; v = d["v-component"][:, :, L]
            we = w - w.mean(axis=1, keepdims=True); ve = v - v.mean(axis=1, keepdims=True)
            M.append((f"w_L{L}", "zonal wind w," + hl, "m/s", w, "sym"))
            M.append((f"v_L{L}", "meridional wind v (+ = southward)," + hl, "m/s", v, "sym"))
            M.append((f"eddy_L{L}", "wind minus its zonal mean (speed, arrows)," + hl, "m/s", np.hypot(we, ve), "lin"))
            QUIVER[f"eddy_L{L}"] = (we, -ve)      # arrow east = +w, north = -v (v points to growing colatitude)
        if "Temperature" in d:
            t = temp_K(planet, d["Temperature"][:, :, L])      # Saturn's panorama writes K/10, the others C
            M.append((f"Teddy_L{L}", "T minus its zonal mean," + hl, "K", t - t.mean(axis=1, keepdims=True), "sym"))
    return [m for m in M if m[3] is not None]

# ------------------------------------------------------------------ geometry
def axes_for(kind, d):
    ns, nf = next(iter(d.values())).shape
    if kind in ("zonal", "proczm"):         # rows = levels i, cols = j
        lat = 90.0 - np.arange(nf) * 180.0 / (nf - 1)
        h = d["height"][:, 0] if "height" in d else np.arange(ns, dtype=float)
        p = d["PressureStat"].mean(axis=1) if "PressureStat" in d else None
        return lat, h, p
    if kind in ("radial", "panorama", "cloudview"):      # rows = j, cols = k
        lon = np.arange(nf) * 360.0 / (nf - 1)
        lat = d["_lat"] if "_lat" in d else 90.0 - np.arange(ns) * 180.0 / (ns - 1)
        return lon, lat, None
    return None

# ------------------------------------------------------------------ one PNG
def draw(job):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    (png, kind, title, unit, Z, scale, rng, X, Y, P, overlays, quiv) = job
    if scale == "rgb":           # cloudview PhotoRGB: an image, no colour bar
        fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=110, constrained_layout=True)
        ax.imshow(np.clip(Z, 0, 1), extent=[X.min(), X.max(), Y.min(), Y.max()], origin="upper", aspect="auto", interpolation="nearest")
        ax.set_xlabel("longitude [deg]", fontsize=8); ax.set_ylabel("latitude [deg]", fontsize=8)
        ax.tick_params(labelsize=7); ax.set_title(title, fontsize=9); fig.savefig(png); plt.close(fig)
        return png
    if kind in ("panorama", "cloudview"):
        fig, (ax, axp) = plt.subplots(1, 2, figsize=(9.6, 3.6), dpi=110, constrained_layout=True,
                                      gridspec_kw={"width_ratios": [3, 1]})
    else:
        fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=110, constrained_layout=True)
    if scale == "log":
        lo, hi = rng; Zp = np.where(Z > lo, Z, np.nan)
        norm = LogNorm(vmin=lo, vmax=hi); cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("#e8e8e8")
    elif scale == "sym":
        lo, hi = rng; Zp = Z; norm = Normalize(lo, hi); cmap = plt.get_cmap("RdBu_r")
    else:
        lo, hi = rng; Zp = Z; norm = Normalize(lo, hi); cmap = plt.get_cmap("inferno")
    im = ax.pcolormesh(X, Y, Zp, shading="nearest", norm=norm, cmap=cmap)
    for (O, col, lab) in overlays or []:
        if O is not None and np.nanmax(O) > 0:
            lev = [max(np.nanmax(O) * 1e-3, 1e-6)]
            ax.contour(X, Y, O, levels=lev, colors=col, linewidths=0.9)
            ax.plot([], [], color=col, lw=0.9, label=lab)
    if overlays and any(o[0] is not None and np.nanmax(o[0]) > 0 for o in overlays):
        ax.legend(loc="upper right", fontsize=7, framealpha=0.7)
    if quiv is not None:
        U, V = quiv; sj = max(1, len(Y) // 14); sk = max(1, len(X) // 20)
        smax = float(np.nanmax(np.hypot(U, V))) or 1.0
        ax.quiver(X[::sk], Y[::sj], U[::sj, ::sk], V[::sj, ::sk], color="#222", width=0.003, pivot="mid",
                  scale=smax * 14, scale_units="width")      # the longest arrow ~ 1/14 of the plot width
    cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label(unit, fontsize=8); cb.ax.tick_params(labelsize=7)
    if kind in ("zonal", "proczm"):
        ax.set_xlim(90, -90); ax.set_xlabel("latitude [deg]", fontsize=8); ax.set_ylabel("height [km]", fontsize=8)
        ax.set_xticks([90, 60, 30, 0, -30, -60, -90])
        if P is not None:
            sec = ax.twinx(); sec.set_ylim(ax.get_ylim()); idx = np.linspace(0, len(Y) - 1, 6).astype(int)
            sec.set_yticks(Y[idx]); sec.set_yticklabels([f"{P[i]:.3g}" for i in idx], fontsize=7)
            sec.set_ylabel("level-mean p [bar]", fontsize=8)
    else:
        ax.set_xlabel("longitude [deg]", fontsize=8); ax.set_ylabel("latitude [deg]", fontsize=8)
        ax.set_ylim(Y.min(), Y.max()) if ZOOM is not None else ax.set_ylim(-90, 90)
    ax.tick_params(labelsize=7); ax.set_title(title, fontsize=9)
    if kind in ("panorama", "cloudview"):       # zonal mean with the longitude min / max as a band (the model runs zonally symmetric)
        zm = np.nanmean(Z, axis=1); lo_, hi_ = np.nanmin(Z, axis=1), np.nanmax(Z, axis=1)
        axp.fill_betweenx(Y, lo_, hi_, color="#9ab", alpha=0.5, lw=0); axp.plot(zm, Y, color="#123", lw=1.2)
        axp.set_ylim(Y.min(), Y.max()) if ZOOM is not None else axp.set_ylim(-90, 90)
        axp.set_xlim(rng[0], rng[1]); axp.set_xlabel(unit, fontsize=8)
        axp.set_title("zonal mean (band: lon min-max)", fontsize=8); axp.tick_params(labelsize=7); axp.grid(alpha=0.3)
    fig.savefig(png); plt.close(fig)
    return png

# ------------------------------------------------------------------ ranges
def joint_range(arrs, scale):
    v = np.concatenate([a[np.isfinite(a)].ravel() for a in arrs if a is not None]) if arrs else np.array([])
    if v.size == 0: return None
    if scale == "sym":
        m = np.percentile(np.abs(v), 99.5); m = m if m > 0 else (np.abs(v).max() or 1.0); return [-float(m), float(m)]
    if scale == "log":
        pos = v[v > 0]
        if pos.size == 0: return None
        hi = float(np.percentile(pos, 99.9)); return [hi * 1e-4, hi]
    return [float(np.percentile(v, 0.5)), float(np.percentile(v, 99.5))]

# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run"); ap.add_argument("run_b", nargs="?")
    ap.add_argument("--out"); ap.add_argument("--kinds", default="zonal,radial,panorama,proczm,cloudview")
    ap.add_argument("--ranges"); ap.add_argument("--every", type=int, default=1, help="keep iterations that are multiples of N")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--levels", help="panorama levels for wind / T maps, e.g. 20,28 (default: the middle level)")
    ap.add_argument("--zoom", help="LAT0,LAT1,LON0,LON1: crop the lat x lon maps to this window")
    a = ap.parse_args()
    global ZOOM
    if a.zoom: ZOOM = tuple(float(x) for x in a.zoom.split(","))
    if a.levels: LEVELS.extend(int(x) for x in a.levels.split(","))
    kinds = a.kinds.split(",")
    oa = find_output(a.run); fa = list_files(oa, kinds)
    fb = list_files(find_output(a.run_b), kinds) if a.run_b else {}
    keys = sorted(k for k in fa if k[2] % a.every == 0 and (not a.run_b or k in fb))
    if not keys: sys.exit("no matching files")
    out = a.out or os.path.join(a.run, ("png" if not a.run_b else "png_vs_" + os.path.basename(os.path.normpath(a.run_b))) + ("_zoom" if ZOOM else ""))
    os.makedirs(out, exist_ok=True)
    ranges = json.load(open(a.ranges)) if a.ranges else {}

    # pass 1: build every panel array (A, B, D) and collect joint ranges
    data = {}; built = {}; pmeta = {}; quivers = {}
    hA = hB = None
    for k in keys:
        f, planet = fa[k]
        if k[0] == "cloudview":
            cA, rA, latA = read_cloudview(f)
            cB, rB, latB = read_cloudview(fb[k][0]) if a.run_b else (None, None, None)
            imc = next(iter(cA.values())).shape[0]
            if hA is None: hA = level_heights_km(oa, imc)
            mB = {m[0]: m for m in cv_maps(cB, rB, hA)} if cB is not None else {}
            for (name, title, unit, A, scale) in cv_maps(cA, rA, hA):
                pmeta[name] = (name, title, unit, None, scale, None)
                built[(k, name)] = (crop(A, latA), crop(mB[name][3], latA) if name in mB else None)
            j_, k_ = next(iter(cA.values())).shape[1:]
            data[k] = (planet, {"_shape": np.zeros((j_, k_)), "_lat": latA}, None); del cA, cB, rA, rB
            continue
        if k[0] == "panorama":
            pA = read_pano(f); pB = read_pano(fb[k][0]) if a.run_b else None
            if not LEVELS: LEVELS.append(next(iter(pA.values())).shape[2] // 2)
            if hA is None: hA = level_heights_km(oa, next(iter(pA.values())).shape[2])
            if a.run_b and hB is None: hB = level_heights_km(find_output(a.run_b), next(iter(pB.values())).shape[2])
            mB = {m[0]: m for m in pano_maps(pB, hB, planet)} if pB is not None else {}
            for (name, title, unit, A, scale) in pano_maps(pA, hA, planet):
                pmeta[name] = (name, title, unit, None, scale, None)
                built[(k, name)] = (crop(A), crop(mB[name][3]) if name in mB else None)
                if name in QUIVER: quivers[(k, name)] = tuple(crop(x) for x in QUIVER[name])
            j_, k_ = next(iter(pA.values())).shape[:2]
            data[k] = (planet, {"_shape": np.zeros((j_, k_))}, None); del pA, pB
            continue
        dA = read_slice(f); dB = read_slice(fb[k][0]) if a.run_b else None
        data[k] = (planet, dA, dB)
    for k, (planet, dA, dB) in data.items():
        kind = k[0]; _KIND[0] = kind
        if kind in ("panorama", "cloudview"): continue
        for (name, title, unit, fn, scale, ov) in PANELS:
            if kind == "radial" and name not in RADIAL_KEYS: continue
            if (kind == "proczm") != (name in PROC_KEYS): continue
            A = fn(dA, planet); B = fn(dB, planet) if dB is not None else None
            if A is None: continue
            if kind == "radial": A, B = crop(A), crop(B)
            built[(k, name)] = (A, B)
    for (name, title, unit, fn, scale, ov) in PANELS + list(pmeta.values()):
        for kind in kinds:
            key = f"{kind}:{name}"
            vals = [x for (k, n), x in built.items() if n == name and k[0] == kind]
            if not vals: continue
            if scale == "rgb": ranges[key] = [0.0, 1.0]; continue
            if key not in ranges:
                r = joint_range([v[0] for v in vals] + [v[1] for v in vals if v[1] is not None], scale)
                if r: ranges[key] = r
            if a.run_b and f"{kind}:D_{name}" not in ranges:
                r = joint_range([v[0] - v[1] for v in vals if v[1] is not None], "sym")
                if r: ranges[f"{kind}:D_{name}"] = r
    json.dump(ranges, open(os.path.join(out, "ranges.json"), "w"), indent=1)

    # pass 2: jobs
    jobs = []; sheet = {}
    meta = {p[0]: p for p in PANELS}; meta.update(pmeta)
    for (k, name), (A, B) in built.items():
        kind, pos, it = k; planet, dA, dB = data[k]
        _, title, unit, _, scale, ov = meta[name]
        if f"{kind}:{name}" not in ranges: continue
        X, Y, P = axes_for(kind, dA)
        if ZOOM is not None and kind in ("radial", "panorama", "cloudview"):
            jj, kk = zoom_idx(len(Y), len(X), Y if kind == "cloudview" else None); X, Y = X[kk], Y[jj]
        os.makedirs(os.path.join(out, kind), exist_ok=True)
        def overlays(d):
            if not ov or kind != "zonal": return None
            return [(d.get(ov[0]), "white", ov[0]), (d.get(ov[1]), "cyan", ov[1])]
        tag = f"{planet} {kind} {pos}  it {it}"
        variants = [("A_" if a.run_b else "", A, dA, os.path.basename(os.path.normpath(a.run)), scale, ranges[f"{kind}:{name}"])]
        if a.run_b and B is not None:
            variants.append(("B_", B, dB, os.path.basename(os.path.normpath(a.run_b)), scale, ranges[f"{kind}:{name}"]))
            if f"{kind}:D_{name}" in ranges:
                variants.append(("D_", A - B, None, "A-B", "sym", ranges[f"{kind}:D_{name}"]))
        for (pre, Z, d, lab, sc, rng) in variants:
            png = os.path.join(out, kind, f"{pre}{name}_{pos}_{it:06d}.png")
            jobs.append((png, kind, f"{title}  [{lab}]  {tag}", unit, Z, sc, rng, X, Y, P, overlays(d) if d is not None else None,
                         quivers.get((k, name)) if pre in ("", "A_") else None))
            sheet.setdefault((pre, kind, pos), {}).setdefault(it, {})[name] = os.path.relpath(png, out)
    with ProcessPoolExecutor(a.jobs) as ex:
        for _ in ex.map(draw, jobs, chunksize=4): pass

    # contact sheet
    order = [p[0] for p in PANELS] + list(pmeta)
    h = ["<!doctype html><meta charset=utf-8><title>Render sheet</title>",
         "<style>body{font:13px system-ui,sans-serif;margin:16px;background:#fff;color:#111}"
         "table{border-collapse:collapse}td,th{padding:2px;border:1px solid #ddd;vertical-align:top}"
         "th{position:sticky;top:0;background:#f4f4f4;font-weight:600}img{width:260px;display:block}"
         "h2{margin:24px 0 6px}nav a{margin-right:12px}</style>",
         f"<h1>{os.path.basename(os.path.normpath(a.run))}" + (f" vs {os.path.basename(os.path.normpath(a.run_b))}" if a.run_b else "") + "</h1>",
         "<p>Colour ranges are fixed over all iterations" + (" and both runs" if a.run_b else "") + " (ranges.json). Click a thumbnail for the full PNG.</p><nav>"]
    secs = sorted(sheet)
    for s in secs: h.append(f"<a href='#{''.join(map(str, s))}'>{({'A_':'A ','B_':'B ','D_':'A-B ','':''})[s[0]]}{s[1]} {s[2]}</a>")
    h.append("</nav>")
    for s in secs:
        rows = sheet[s]; cols = [c for c in order if any(c in r for r in rows.values())]
        h.append(f"<h2 id='{''.join(map(str, s))}'>{({'A_':'A: ','B_':'B: ','D_':'A-B: ','':''})[s[0]]}{s[1]} {s[2]}</h2><table><tr><th>it</th>"
                 + "".join(f"<th>{c}</th>" for c in cols) + "</tr>")
        for it in sorted(rows):
            h.append(f"<tr><th>{it}</th>" + "".join(
                f"<td><a href='{rows[it][c]}'><img loading=lazy src='{rows[it][c]}'></a></td>" if c in rows[it] else "<td></td>" for c in cols) + "</tr>")
        h.append("</table>")
    open(os.path.join(out, "index.html"), "w").write("\n".join(h))
    print(f"{len(jobs)} PNGs, {len(keys)} slice files -> {out}/index.html")

if __name__ == "__main__":
    main()
