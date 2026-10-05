# -*- coding: utf-8 -*-
"""
Decadal Shoreline Change Tool
=============================
ArcGIS Pro Python Toolbox (.pyt)

Extracts shorelines from two dates of aerial/satellite imagery (e.g., NAIP),
measures erosion/accretion along perpendicular transects, and tags the results
with restoration-relevant planning layers.

Tools (run in this order):
  1. Extract Shoreline          imagery -> water raster -> water polygon -> shoreline line
  2. Generate Transects         perpendicular transects along the baseline (Date 1) shoreline
  3. Calculate Shoreline Change erosion/accretion distance and rate (m/yr) per transect
  4. Restoration Overlay        flag transects that touch habitat, ownership, project layers

Sign convention (same as USGS DSAS):
  negative = erosion   (shoreline moved landward)
  positive = accretion (shoreline moved waterward)

Requirements
  - ArcGIS Pro with arcpy
  - Spatial Analyst extension (Tool 1 only)
  - All inputs in a PROJECTED coordinate system (UTM, State Plane, etc.)
  - Distances in the tool dialogs are entered in meters

The workflow avoids Advanced-license tools (Polygon To Line, Eliminate Polygon Part,
Erase) by doing that geometry work directly in arcpy, so it runs at a Basic license.

License: MIT
"""

import os
import arcpy


# =============================================================================
# Shared helpers
# =============================================================================

def _fail(message):
    """Report an error in the GP messages and stop the tool."""
    arcpy.AddError(message)
    raise arcpy.ExecuteError(message)


def _scratch_fc(name):
    return arcpy.CreateUniqueName(name, arcpy.env.scratchGDB)


def _scratch_tif(name):
    return arcpy.CreateUniqueName(name + ".tif", arcpy.env.scratchFolder)


def _delete(*items):
    for item in items:
        try:
            if item and arcpy.Exists(item):
                arcpy.management.Delete(item)
        except Exception:
            pass


def _projected_sr(dataset, label):
    """Return the dataset's spatial reference, or stop if it isn't projected."""
    sr = arcpy.Describe(dataset).spatialReference
    if sr is None or sr.type != "Projected":
        _fail(
            f"{label} must be in a projected coordinate system (e.g., UTM or "
            "State Plane) so distances are meaningful. Project it and rerun."
        )
    return sr


def _is_shapefile(path):
    return path.lower().endswith(".shp")


def _create_fc(path, geom_type, sr, fields=()):
    """Create a feature class with fields given as (name, type[, length])."""
    ws, name = os.path.split(path)
    if arcpy.Exists(path):
        arcpy.management.Delete(path)
    arcpy.management.CreateFeatureclass(ws, name, geom_type, spatial_reference=sr)
    for f in fields:
        if len(f) > 2:
            arcpy.management.AddField(path, f[0], f[1], field_length=f[2])
        else:
            arcpy.management.AddField(path, f[0], f[1])


def _dissolve(features, sr):
    """Dissolve all features into a single scratch feature class in `sr`."""
    out = _scratch_fc("dslv")
    with arcpy.EnvManager(outputCoordinateSystem=sr):
        arcpy.management.Dissolve(features, out)
    return out


def _read_union(fc):
    """Read every geometry in a feature class and union them into one."""
    geom = None
    with arcpy.da.SearchCursor(fc, ["SHAPE@"]) as cur:
        for (g,) in cur:
            if g is None:
                continue
            geom = g if geom is None else geom.union(g)
    return geom


def _dissolve_to_geometry(features, sr):
    tmp = _dissolve(features, sr)
    geom = _read_union(tmp)
    _delete(tmp)
    return geom


_band_layer_counter = [0]


def _band(raster_path, band_index):
    """Return one band of a multiband raster as a float Raster object."""
    _band_layer_counter[0] += 1
    lyr = f"dsct_band{band_index}_{_band_layer_counter[0]}"
    arcpy.management.MakeRasterLayer(raster_path, lyr, band_index=str(band_index))
    return arcpy.sa.Float(lyr)


def _split_rings(part):
    """Split one polygon part (Array of Points with None separators) into rings."""
    rings, ring = [], []
    for pt in part:
        if pt is None:
            if ring:
                rings.append(ring)
            ring = []
        else:
            ring.append(pt)
    if ring:
        rings.append(ring)
    return rings


def _fill_small_holes(poly, min_hole_area, sr):
    """
    Remove interior rings (land 'holes' inside the water polygon) smaller than
    min_hole_area (map units squared). Larger holes, e.g. real islands, are kept.
    """
    if min_hole_area <= 0:
        return poly

    kept_rings = []
    removed = 0
    for part in poly:
        rings = _split_rings(part)
        if not rings:
            continue
        kept_rings.append(rings[0])  # exterior ring
        for hole in rings[1:]:
            if arcpy.Polygon(arcpy.Array(hole), sr).area >= min_hole_area:
                kept_rings.append(hole)
            else:
                removed += 1

    if removed == 0:
        return poly
    # arcpy simplifies the rebuilt polygon: rings nested inside an exterior
    # ring become holes again.
    return arcpy.Polygon(arcpy.Array([arcpy.Array(r) for r in kept_rings]), sr)


def _point_from_multipoint(mp, i):
    p = mp.getPart(i)
    if isinstance(p, arcpy.Array):
        p = p[0]
    return p


# =============================================================================
# Core workflow functions (usable from a standalone script as well)
# =============================================================================

def extract_shoreline(
    in_raster,
    out_shoreline,
    method="NDWI",
    threshold=0.0,
    green_band=2,
    nir_band=4,
    red_band=1,
    blue_band=3,
    reference_line=None,
    buffer_m=300.0,
    analysis_area=None,
    keep_largest=True,
    min_water_m2=10000.0,
    min_island_m2=500.0,
    majority_filter=True,
    smooth_m=10.0,
    min_length_m=50.0,
    out_water=None,
):
    """
    Imagery -> water raster -> water polygon -> cleaned shoreline line.

    method: "NDWI" (water where NDWI > threshold) or
            "BRIGHTNESS" (water where mean RGB brightness < threshold)
    NAIP band order is 1=Red, 2=Green, 3=Blue, 4=NIR.
    """
    arcpy.CheckOutExtension("Spatial")
    scratch = []

    desc = arcpy.Describe(in_raster)
    sr = _projected_sr(in_raster, "Input imagery")
    mpu = sr.metersPerUnit
    raster_path = desc.catalogPath
    ras = arcpy.Raster(raster_path)
    cell = ras.meanCellWidth

    needed = [green_band, nir_band] if method == "NDWI" else [red_band, green_band, blue_band]
    for b in needed:
        if b < 1 or b > ras.bandCount:
            _fail(f"Band {b} doesn't exist; this image has {ras.bandCount} band(s).")

    try:
        with arcpy.EnvManager(
            outputCoordinateSystem=sr, snapRaster=raster_path, cellSize=raster_path
        ):
            # --- 1. Analysis area ------------------------------------------------
            mask_fc = None
            if analysis_area:
                arcpy.AddMessage("Using the supplied analysis area polygon.")
                mask_fc = _dissolve(analysis_area, sr)
                scratch.append(mask_fc)
            elif reference_line:
                arcpy.AddMessage(f"Buffering reference shoreline by {buffer_m} m.")
                buf = _scratch_fc("buffer")
                scratch.append(buf)
                arcpy.analysis.Buffer(
                    reference_line, buf, f"{buffer_m} Meters", dissolve_option="ALL"
                )
                mask_fc = buf
            else:
                arcpy.AddWarning(
                    "No reference shoreline or analysis area given; processing the "
                    "full image extent. This can be slow for large mosaics."
                )

            if mask_fc:
                mask_geom = _read_union(mask_fc)
                ext = mask_geom.extent
                clipped = _scratch_tif("clip")
                scratch.append(clipped)
                arcpy.AddMessage("Clipping imagery to the analysis area...")
                arcpy.management.Clip(
                    raster_path,
                    f"{ext.XMin} {ext.YMin} {ext.XMax} {ext.YMax}",
                    clipped,
                    mask_fc,
                    "",
                    "ClippingGeometry",
                    "NO_MAINTAIN_EXTENT",
                )
                work_raster = clipped
            else:
                e = ras.extent
                mask_geom = arcpy.Polygon(
                    arcpy.Array(
                        [
                            arcpy.Point(e.XMin, e.YMin),
                            arcpy.Point(e.XMin, e.YMax),
                            arcpy.Point(e.XMax, e.YMax),
                            arcpy.Point(e.XMax, e.YMin),
                        ]
                    ),
                    sr,
                )
                work_raster = raster_path

            # --- 2. Water classification -----------------------------------------
            if method == "NDWI":
                arcpy.AddMessage(f"Calculating NDWI; water where NDWI > {threshold}.")
                g = _band(work_raster, green_band)
                n = _band(work_raster, nir_band)
                denom = g + n
                index = arcpy.sa.Con(denom == 0, 0.0, (g - n) / denom)
                is_water = index > threshold
            else:
                arcpy.AddMessage(f"Calculating brightness; water where brightness < {threshold}.")
                r = _band(work_raster, red_band)
                g = _band(work_raster, green_band)
                b = _band(work_raster, blue_band)
                index = (r + g + b) / 3.0
                is_water = index < threshold

            water = arcpy.sa.Con(is_water, 1, 0)
            if majority_filter:
                water = arcpy.sa.MajorityFilter(water, "EIGHT", "MAJORITY")
            if mask_fc:
                # Guarantees cells outside the analysis area are NoData, not 0.
                water = arcpy.sa.ExtractByMask(water, mask_fc)
            water = arcpy.sa.SetNull(water != 1, 1)

            water_tif = _scratch_tif("water")
            scratch.append(water_tif)
            water.save(water_tif)

            # --- 3. Water polygon --------------------------------------------------
            arcpy.AddMessage("Converting water raster to polygons...")
            poly_fc = _scratch_fc("waterpoly")
            scratch.append(poly_fc)
            arcpy.conversion.RasterToPolygon(
                water_tif, poly_fc, "SIMPLIFY", "Value", "SINGLE_OUTER_PART"
            )
            polys = [row[0] for row in arcpy.da.SearchCursor(poly_fc, ["SHAPE@"])]
            if not polys:
                _fail(
                    "No water was detected. Try adjusting the threshold "
                    "(NDWI: lower it; brightness: raise it)."
                )

            if keep_largest:
                polys = [max(polys, key=lambda p: p.area)]
                arcpy.AddMessage("Kept the largest water body.")
            else:
                min_water = min_water_m2 / mpu ** 2
                polys = [p for p in polys if p.area >= min_water]
                arcpy.AddMessage(f"Kept {len(polys)} water bodies >= {min_water_m2} m².")
                if not polys:
                    _fail("No water bodies are larger than the minimum water area.")

            min_hole = min_island_m2 / mpu ** 2
            polys = [_fill_small_holes(p, min_hole, sr) for p in polys]

            if out_water:
                _create_fc(out_water, "POLYGON", sr, [("AREA_M2", "DOUBLE")])
                with arcpy.da.InsertCursor(out_water, ["SHAPE@", "AREA_M2"]) as cur:
                    for p in polys:
                        cur.insertRow([p, p.area * mpu ** 2])

            # --- 4. Shoreline line --------------------------------------------------
            # Drop the edges where water meets the analysis-area boundary;
            # those aren't real shoreline.
            inner = mask_geom.buffer(-2 * cell)
            if inner is None or inner.area <= 0:
                inner = mask_geom

            raw_lines = _scratch_fc("rawline")
            scratch.append(raw_lines)
            _create_fc(raw_lines, "POLYLINE", sr)
            with arcpy.da.InsertCursor(raw_lines, ["SHAPE@"]) as cur:
                for p in polys:
                    line = p.boundary().intersect(inner, 2)
                    if line is not None and line.length > 0:
                        cur.insertRow([line])

            # --- 5. Smoothing -------------------------------------------------------
            smooth_fc = raw_lines
            if smooth_m and smooth_m > 0:
                arcpy.AddMessage(f"Smoothing shoreline (PAEK, {smooth_m} m tolerance)...")
                smooth_fc = _scratch_fc("smoothline")
                scratch.append(smooth_fc)
                try:
                    arcpy.cartography.SmoothLine(
                        raw_lines, smooth_fc, "PAEK", f"{smooth_m} Meters"
                    )
                except arcpy.ExecuteError:
                    arcpy.AddWarning("Smooth Line failed; keeping the unsmoothed shoreline.")
                    smooth_fc = raw_lines

            # --- 6. Explode to single parts and drop short fragments ---------------
            _create_fc(
                out_shoreline,
                "POLYLINE",
                sr,
                [
                    ("LENGTH_M", "DOUBLE"),
                    ("SOURCE", "TEXT", 255),
                    ("METHOD", "TEXT", 20),
                    ("THRESHOLD", "DOUBLE"),
                ],
            )
            min_len = min_length_m / mpu
            source = os.path.basename(raster_path)[:255]
            kept = 0
            with arcpy.da.SearchCursor(smooth_fc, ["SHAPE@"]) as s_cur, arcpy.da.InsertCursor(
                out_shoreline, ["SHAPE@", "LENGTH_M", "SOURCE", "METHOD", "THRESHOLD"]
            ) as i_cur:
                for (geom,) in s_cur:
                    if geom is None:
                        continue
                    for i in range(geom.partCount):
                        part = arcpy.Polyline(geom.getPart(i), sr)
                        if part.length >= min_len:
                            i_cur.insertRow(
                                [part, part.length * mpu, source, method, threshold]
                            )
                            kept += 1

            if kept == 0:
                arcpy.AddWarning(
                    "No shoreline segments survived the minimum length filter. "
                    "Check the threshold or lower the minimum segment length."
                )
            else:
                arcpy.AddMessage(f"Wrote {kept} shoreline segment(s) to {out_shoreline}.")
    finally:
        _delete(*scratch)


def generate_transects(baseline, out_transects, spacing_m=50.0, length_m=300.0, include_ends=False):
    """
    Perpendicular transects along the baseline (earlier) shoreline.
    Transects are centered on the shoreline, so a 300 m transect reaches
    150 m landward and 150 m waterward.
    """
    _projected_sr(baseline, "Baseline shoreline")
    arcpy.AddMessage(
        f"Generating {length_m} m transects every {spacing_m} m along the baseline..."
    )
    arcpy.management.GenerateTransectsAlongLines(
        baseline,
        out_transects,
        f"{spacing_m} Meters",
        f"{length_m} Meters",
        "END_POINTS" if include_ends else "NO_END_POINTS",
    )
    arcpy.management.AddField(out_transects, "TRANSECTID", "LONG")
    count = 0
    with arcpy.da.UpdateCursor(out_transects, ["TRANSECTID"]) as cur:
        for count, _row in enumerate(cur, start=1):
            cur.updateRow([count])
    arcpy.AddMessage(f"Created {count} transects.")


def calculate_change(
    transects,
    baseline_water,
    comparison_shoreline,
    date1,
    date2,
    out_fc,
    stable_m=1.0,
    out_vectors=None,
):
    """
    For each transect:
      - the baseline position is the transect midpoint (where it crosses the Date 1 shoreline)
      - the Date 2 position is the nearest intersection with the comparison shoreline
      - if the Date 2 position falls inside the Date 1 water polygon, land moved
        into former water -> accretion (+); otherwise -> erosion (-)
    """
    sr = _projected_sr(transects, "Transects")
    mpu = sr.metersPerUnit

    years = (date2 - date1).days / 365.25
    if years <= 0:
        _fail("Date 2 must be later than Date 1.")
    arcpy.AddMessage(f"Time span: {years:.2f} years.")

    arcpy.AddMessage("Preparing baseline water polygon and comparison shoreline...")
    water_geom = _dissolve_to_geometry(baseline_water, sr)
    comp_geom = _dissolve_to_geometry(comparison_shoreline, sr)
    if water_geom is None or comp_geom is None:
        _fail("The baseline water polygon or comparison shoreline is empty.")

    arcpy.management.CopyFeatures(transects, out_fc)
    existing = [f.name.upper() for f in arcpy.ListFields(out_fc)]
    id_field = "TRANSECTID" if "TRANSECTID" in existing else arcpy.Describe(out_fc).OIDFieldName

    new_fields = [
        ["DATE1", "DATE"],
        ["DATE2", "DATE"],
        ["YEARS", "DOUBLE"],
        ["DIST_M", "DOUBLE"],
        ["ABS_DIST_M", "DOUBLE"],
        ["RATE_M_YR", "DOUBLE"],
        ["CHANGE", "TEXT", "", 20],
        ["N_INTERSX", "SHORT"],
    ]
    arcpy.management.AddFields(out_fc, new_fields)

    if out_vectors:
        _create_fc(
            out_vectors,
            "POLYLINE",
            sr,
            [
                ("TRANSECTID", "LONG"),
                ("DIST_M", "DOUBLE"),
                ("RATE_M_YR", "DOUBLE"),
                ("CHANGE", "TEXT", 20),
            ],
        )
        vec_cur = arcpy.da.InsertCursor(
            out_vectors, ["SHAPE@", "TRANSECTID", "DIST_M", "RATE_M_YR", "CHANGE"]
        )
    else:
        vec_cur = None

    total = int(arcpy.management.GetCount(out_fc)[0])
    arcpy.SetProgressor("step", "Measuring shoreline change...", 0, max(total, 1), 1)

    stats = {"Erosion": 0, "Accretion": 0, "Stable": 0, "No Intersect": 0}
    rates = []
    fields = [
        "SHAPE@", id_field, "DATE1", "DATE2", "YEARS", "DIST_M",
        "ABS_DIST_M", "RATE_M_YR", "CHANGE", "N_INTERSX",
    ]

    try:
        with arcpy.da.UpdateCursor(out_fc, fields) as cur:
            for row in cur:
                t = row[0]
                row[2], row[3], row[4] = date1, date2, years

                if t is None or t.length == 0:
                    row[5:10] = [None, None, None, "No Intersect", 0]
                    stats["No Intersect"] += 1
                    cur.updateRow(row)
                    arcpy.SetProgressorPosition()
                    continue

                anchor = t.positionAlongLine(0.5, True)
                hits = t.intersect(comp_geom, 1)
                n_hits = hits.pointCount if hits is not None else 0

                if n_hits == 0:
                    row[5:10] = [None, None, None, "No Intersect", 0]
                    stats["No Intersect"] += 1
                else:
                    pts = [
                        arcpy.PointGeometry(_point_from_multipoint(hits, i), sr)
                        for i in range(n_hits)
                    ]
                    nearest = min(pts, key=anchor.distanceTo)
                    abs_m = anchor.distanceTo(nearest) * mpu
                    sign = 1.0 if water_geom.contains(nearest) else -1.0
                    dist_m = sign * abs_m
                    rate = dist_m / years

                    if abs_m < stable_m:
                        change = "Stable"
                    else:
                        change = "Accretion" if sign > 0 else "Erosion"

                    row[5:10] = [dist_m, abs_m, rate, change, n_hits]
                    stats[change] += 1
                    rates.append(rate)

                    if vec_cur is not None:
                        vec = arcpy.Polyline(
                            arcpy.Array([anchor.firstPoint, nearest.firstPoint]), sr
                        )
                        vec_cur.insertRow([vec, row[1], dist_m, rate, change])

                cur.updateRow(row)
                arcpy.SetProgressorPosition()
    finally:
        if vec_cur is not None:
            del vec_cur
        arcpy.ResetProgressor()

    arcpy.AddMessage("---- Shoreline change summary ----")
    for k, v in stats.items():
        arcpy.AddMessage(f"  {k}: {v}")
    if rates:
        arcpy.AddMessage(f"  Mean rate: {sum(rates) / len(rates):.2f} m/yr")
        arcpy.AddMessage(f"  Max erosion rate: {min(rates):.2f} m/yr")
        arcpy.AddMessage(f"  Max accretion rate: {max(rates):.2f} m/yr")
    if stats["No Intersect"]:
        arcpy.AddWarning(
            f"{stats['No Intersect']} transect(s) didn't cross the Date 2 shoreline. "
            "Longer transects or a better Date 2 extraction may help."
        )


def restoration_overlay(change_fc, overlay_layers, out_fc, search_m=0.0, erosion_threshold=0.5):
    """
    Adds a 0/1 flag field per overlay layer, an overlay count, and a simple
    PRIORITY class:
      High     - erosion faster than the threshold AND touches >= 1 overlay layer
      Moderate - erosion faster than the threshold, OR any erosion that touches an overlay
      Low      - everything else
    The scoring is meant as a starting point; adjust it to your planning criteria.
    """
    arcpy.management.CopyFeatures(change_fc, out_fc)
    names = [f.name.upper() for f in arcpy.ListFields(out_fc)]
    if "RATE_M_YR" not in names or "CHANGE" not in names:
        _fail("Input must be the output of Calculate Shoreline Change (needs RATE_M_YR and CHANGE).")

    ws = os.path.dirname(out_fc)
    max_len = 10 if _is_shapefile(out_fc) else 30
    lyr = "dsct_overlay_lyr"
    arcpy.management.MakeFeatureLayer(out_fc, lyr)
    distance = f"{search_m} Meters" if search_m and search_m > 0 else ""

    flag_fields = []
    try:
        for layer in overlay_layers:
            layer_name = arcpy.Describe(layer).name
            base = os.path.splitext(layer_name)[0]
            fname = arcpy.ValidateFieldName(("OV_" + base)[:max_len], ws)
            suffix = 1
            while fname.upper() in names:
                tail = f"_{suffix}"
                fname = arcpy.ValidateFieldName(("OV_" + base)[: max_len - len(tail)] + tail, ws)
                suffix += 1
            names.append(fname.upper())
            flag_fields.append(fname)

            arcpy.management.AddField(out_fc, fname, "SHORT", field_alias=f"Touches {base}")
            arcpy.management.SelectLayerByLocation(lyr, "INTERSECT", layer, distance, "NEW_SELECTION")
            hit_ids = {r[0] for r in arcpy.da.SearchCursor(lyr, ["OID@"])}
            arcpy.management.SelectLayerByAttribute(lyr, "CLEAR_SELECTION")

            with arcpy.da.UpdateCursor(out_fc, ["OID@", fname]) as cur:
                for oid, _ in cur:
                    cur.updateRow([oid, 1 if oid in hit_ids else 0])
            arcpy.AddMessage(f"{base}: {len(hit_ids)} transect(s) flagged as {fname}.")
    finally:
        _delete(lyr)

    arcpy.management.AddFields(out_fc, [["OVERLAY_CT", "SHORT"], ["PRIORITY", "TEXT", "", 10]])

    counts = {"High": 0, "Moderate": 0, "Low": 0}
    with arcpy.da.UpdateCursor(out_fc, ["RATE_M_YR", "CHANGE", "OVERLAY_CT", "PRIORITY"] + flag_fields) as cur:
        for row in cur:
            rate, change = row[0], row[1]
            ct = sum(v or 0 for v in row[4:])
            fast_erosion = rate is not None and rate <= -abs(erosion_threshold)
            if fast_erosion and ct > 0:
                priority = "High"
            elif fast_erosion or (change == "Erosion" and ct > 0):
                priority = "Moderate"
            else:
                priority = "Low"
            row[2], row[3] = ct, priority
            counts[priority] += 1
            cur.updateRow(row)

    arcpy.AddMessage(
        f"Priority: High {counts['High']}, Moderate {counts['Moderate']}, Low {counts['Low']}"
    )


# =============================================================================
# Toolbox definition
# =============================================================================

def _param(name, label, dtype, ptype="Required", direction="Input", value=None, multi=False, category=None):
    p = arcpy.Parameter(
        displayName=label,
        name=name,
        datatype=dtype,
        parameterType=ptype,
        direction=direction,
        multiValue=multi,
    )
    if value is not None:
        p.value = value
    if category:
        p.category = category
    return p


class Toolbox(object):
    def __init__(self):
        self.label = "Decadal Shoreline Change"
        self.alias = "shorelinechange"
        self.tools = [ExtractShoreline, GenerateTransects, CalculateShorelineChange, RestorationOverlay]


class ExtractShoreline(object):
    METHOD_NDWI = "NDWI (Green/NIR)"
    METHOD_BRIGHT = "Brightness (RGB only)"
    KEEP_LARGEST = "Largest water body only"
    KEEP_ALL = "All water bodies above minimum area"

    def __init__(self):
        self.label = "1. Extract Shoreline"
        self.description = (
            "Classifies water vs. land from imagery (NDWI or RGB brightness threshold) "
            "and converts the water edge into a cleaned shoreline line. Run once per imagery date."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        cleanup = "Shoreline Cleanup"
        bands = "Band Assignment (NAIP: R=1, G=2, B=3, NIR=4)"

        p_ref = _param("reference_line", "Reference Shoreline / Lagoon Edge (for buffer)", "GPFeatureLayer", "Optional")
        p_ref.filter.list = ["Polyline"]
        p_area = _param("analysis_area", "Analysis Area Polygon (overrides buffer)", "GPFeatureLayer", "Optional")
        p_area.filter.list = ["Polygon"]

        p_method = _param("method", "Extraction Method", "GPString", value=self.METHOD_NDWI)
        p_method.filter.type = "ValueList"
        p_method.filter.list = [self.METHOD_NDWI, self.METHOD_BRIGHT]

        p_keep = _param("keep", "Water Bodies to Keep", "GPString", value=self.KEEP_LARGEST, category=cleanup)
        p_keep.filter.type = "ValueList"
        p_keep.filter.list = [self.KEEP_LARGEST, self.KEEP_ALL]

        return [
            _param("in_raster", "Input Imagery", "GPRasterLayer"),                                  # 0
            _param("out_shoreline", "Output Shoreline", "DEFeatureClass", direction="Output"),       # 1
            p_method,                                                                                 # 2
            _param("threshold", "Water Threshold (NDWI > value, or Brightness < value)", "GPDouble", value=0.0),  # 3
            _param("green_band", "Green Band", "GPLong", value=2, category=bands),                    # 4
            _param("nir_band", "NIR Band", "GPLong", value=4, category=bands),                        # 5
            _param("red_band", "Red Band", "GPLong", value=1, category=bands),                        # 6
            _param("blue_band", "Blue Band", "GPLong", value=3, category=bands),                      # 7
            p_ref,                                                                                    # 8
            _param("buffer_m", "Buffer Distance (m)", "GPDouble", value=300.0),                       # 9
            p_area,                                                                                   # 10
            p_keep,                                                                                   # 11
            _param("min_water_m2", "Minimum Water Body Area (m²)", "GPDouble", value=10000.0, category=cleanup),   # 12
            _param("min_island_m2", "Fill Land Patches Smaller Than (m²)", "GPDouble", value=500.0, category=cleanup),  # 13
            _param("majority", "Apply Majority Filter (removes speckle)", "GPBoolean", value=True, category=cleanup),  # 14
            _param("smooth_m", "Smoothing Tolerance (m, 0 = none)", "GPDouble", value=10.0, category=cleanup),       # 15
            _param("min_length_m", "Minimum Shoreline Segment Length (m)", "GPDouble", value=50.0, category=cleanup),  # 16
            _param("out_water", "Output Water Polygon (needed for the baseline date)", "DEFeatureClass", "Optional", "Output"),  # 17
        ]

    def isLicensed(self):
        return arcpy.CheckExtension("Spatial") == "Available"

    def updateParameters(self, parameters):
        ndwi = (parameters[2].valueAsText or "").startswith("NDWI")
        parameters[5].enabled = ndwi
        parameters[6].enabled = not ndwi
        parameters[7].enabled = not ndwi
        if parameters[2].altered and not parameters[2].hasBeenValidated:
            parameters[3].value = 0.0 if ndwi else 30.0
        parameters[12].enabled = parameters[11].valueAsText == self.KEEP_ALL
        parameters[9].enabled = parameters[10].value is None

    def updateMessages(self, parameters):
        if parameters[0].value:
            try:
                count = arcpy.Describe(parameters[0].valueAsText).bandCount
                ndwi = (parameters[2].valueAsText or "").startswith("NDWI")
                idx = [4, 5] if ndwi else [6, 4, 7]
                for i in idx:
                    if parameters[i].value and parameters[i].value > count:
                        parameters[i].setErrorMessage(f"The input imagery only has {count} band(s).")
                if ndwi and count < 4:
                    parameters[2].setWarningMessage(
                        "This imagery has fewer than 4 bands. NDWI needs a NIR band; "
                        "consider the Brightness (RGB only) method."
                    )
            except Exception:
                pass

    def execute(self, parameters, messages):
        p = parameters
        method = "NDWI" if p[2].valueAsText.startswith("NDWI") else "BRIGHTNESS"
        extract_shoreline(
            in_raster=p[0].valueAsText,
            out_shoreline=p[1].valueAsText,
            method=method,
            threshold=float(p[3].value),
            green_band=int(p[4].value),
            nir_band=int(p[5].value),
            red_band=int(p[6].value),
            blue_band=int(p[7].value),
            reference_line=p[8].valueAsText,
            buffer_m=float(p[9].value or 300.0),
            analysis_area=p[10].valueAsText,
            keep_largest=p[11].valueAsText == self.KEEP_LARGEST,
            min_water_m2=float(p[12].value or 0),
            min_island_m2=float(p[13].value or 0),
            majority_filter=bool(p[14].value),
            smooth_m=float(p[15].value or 0),
            min_length_m=float(p[16].value or 0),
            out_water=p[17].valueAsText,
        )


class GenerateTransects(object):
    def __init__(self):
        self.label = "2. Generate Transects"
        self.description = (
            "Creates evenly spaced transects perpendicular to the baseline (earlier) shoreline. "
            "Transects are centered on the shoreline."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_base = _param("baseline", "Baseline (Date 1) Shoreline", "GPFeatureLayer")
        p_base.filter.list = ["Polyline"]
        return [
            p_base,
            _param("out_transects", "Output Transects", "DEFeatureClass", direction="Output"),
            _param("spacing_m", "Transect Spacing (m)", "GPDouble", value=50.0),
            _param("length_m", "Transect Length (m, centered on shoreline)", "GPDouble", value=300.0),
            _param("include_ends", "Include Transects at Line Ends", "GPBoolean", value=False),
        ]

    def updateMessages(self, parameters):
        for i in (2, 3):
            if parameters[i].value is not None and parameters[i].value <= 0:
                parameters[i].setErrorMessage("Must be greater than 0.")

    def execute(self, parameters, messages):
        generate_transects(
            parameters[0].valueAsText,
            parameters[1].valueAsText,
            float(parameters[2].value),
            float(parameters[3].value),
            bool(parameters[4].value),
        )


class CalculateShorelineChange(object):
    def __init__(self):
        self.label = "3. Calculate Shoreline Change"
        self.description = (
            "Measures distance and annualized rate of shoreline change along each transect "
            "and classifies it as erosion, accretion, or stable."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_tr = _param("transects", "Transects", "GPFeatureLayer")
        p_tr.filter.list = ["Polyline"]
        p_water = _param("baseline_water", "Baseline (Date 1) Water Polygon", "GPFeatureLayer")
        p_water.filter.list = ["Polygon"]
        p_comp = _param("comparison", "Comparison (Date 2) Shoreline", "GPFeatureLayer")
        p_comp.filter.list = ["Polyline"]
        return [
            p_tr,                                                                            # 0
            p_water,                                                                         # 1
            p_comp,                                                                          # 2
            _param("date1", "Date 1 (baseline imagery date)", "GPDate"),                     # 3
            _param("date2", "Date 2 (comparison imagery date)", "GPDate"),                   # 4
            _param("out_fc", "Output Change Transects", "DEFeatureClass", direction="Output"),  # 5
            _param("stable_m", "Stable Threshold (m) - changes smaller than this are 'Stable'", "GPDouble", value=1.0),  # 6
            _param("out_vectors", "Output Change Vectors (optional)", "DEFeatureClass", "Optional", "Output"),  # 7
        ]

    def updateMessages(self, parameters):
        d1, d2 = parameters[3].value, parameters[4].value
        if d1 and d2 and d2 <= d1:
            parameters[4].setErrorMessage("Date 2 must be later than Date 1.")

    def execute(self, parameters, messages):
        p = parameters
        calculate_change(
            transects=p[0].valueAsText,
            baseline_water=p[1].valueAsText,
            comparison_shoreline=p[2].valueAsText,
            date1=p[3].value,
            date2=p[4].value,
            out_fc=p[5].valueAsText,
            stable_m=float(p[6].value or 0),
            out_vectors=p[7].valueAsText,
        )


class RestorationOverlay(object):
    def __init__(self):
        self.label = "4. Restoration Overlay"
        self.description = (
            "Flags change transects that intersect habitat, ownership, project, or other "
            "planning layers and assigns a simple restoration priority class."
        )
        self.canRunInBackground = False

    def getParameterInfo(self):
        p_in = _param("change_fc", "Change Transects (from Tool 3)", "GPFeatureLayer")
        p_in.filter.list = ["Polyline"]
        return [
            p_in,
            _param("overlays", "Overlay Layers (ownership, habitat, projects...)", "GPFeatureLayer", multi=True),
            _param("out_fc", "Output Prioritized Transects", "DEFeatureClass", direction="Output"),
            _param("search_m", "Search Distance (m)", "GPDouble", value=0.0),
            _param("erosion_threshold", "Significant Erosion Rate (m/yr, positive number)", "GPDouble", value=0.5),
        ]

    def execute(self, parameters, messages):
        restoration_overlay(
            change_fc=parameters[0].valueAsText,
            overlay_layers=parameters[1].values,
            out_fc=parameters[2].valueAsText,
            search_m=float(parameters[3].value or 0),
            erosion_threshold=float(parameters[4].value or 0),
        )
