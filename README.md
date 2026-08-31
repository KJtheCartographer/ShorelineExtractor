# Decadal Shoreline Change Tool

An ArcGIS Pro Python toolbox that extracts shorelines from aerial/satellite imagery (including free NAIP imagery) and quantifies lateral shoreline movement — erosion or accretion — over time. Built to support coastal restoration planning by linking shoreline-change results to habitat, ownership, and planning layers.

## Overview

Rather than simply showing *where* a shoreline has moved, this tool is designed to help answer a more useful planning question:

> **Where might shoreline change create an opportunity or need for restoration?**

It does this by extracting shoreline position from two imagery dates, measuring change along evenly spaced transects, calculating an annualized rate of change, and enabling that output to be overlaid with restoration-relevant spatial layers.

## Methodology

**1. Imagery**
- Uses two dates of aerial/satellite imagery (e.g., April 2016 and April 2025 in the current workflow).
- Imagery is processed natively in ArcGIS Pro.

**2. Analysis Area**
- Creates a 300-meter inland buffer from the shoreline/lagoon edge.
- Clips both imagery dates to this buffer to reduce processing time and keep the analysis focused on the shoreline zone.

**3. Shoreline Extraction**

*Standard workflow (multispectral imagery with NIR):*
- Calculates NDWI: `NDWI = (Green − NIR) / (Green + NIR)`
- Applies a threshold (approximately `NDWI > 0`) to classify water vs. land.
- Converts the water raster to polygons, isolates the primary waterbody, and converts it to a shoreline line.
- Cleans and smooths the shoreline to remove small artifacts.

*High-resolution / RGB-only workflow (no usable NIR band):*
- Calculates a brightness raster: `Brightness = (Red + Green + Blue) / 3`
- Applies an empirically adjusted threshold (e.g., `Con("Brightness" < 30, 1, 0)`)
- Follows the same downstream path: Water Raster → Water Polygon → Shoreline Line → Transects → Shoreline Change
- Iso Cluster unsupervised classification and supervised classification (Random Trees, SVM) are noted as alternative approaches when a simple brightness threshold isn't sufficient.

**4. Transect Generation**
- Generates transects spaced 50 meters apart along the baseline (earlier) shoreline.
- Transects are oriented perpendicular to the shoreline.
- Transects are approximately 300 meters long — long enough to intersect both shoreline positions.

**5. Shoreline Change Calculation**
- Intersects each transect with both shoreline positions.
- Measures the distance between the two positions along each transect.
- Classifies direction of change:
  - **Erosion** — shoreline moved landward
  - **Accretion** — shoreline moved waterward
- Converts total distance to an annualized rate:
  `Shoreline Change Rate = Distance Changed ÷ Number of Years`
  (output in meters/year)

**6. Restoration & Planning Application**
- Shoreline-change results can be overlaid with:
  - Public ownership boundaries
  - Existing habitat
  - Active or planned restoration projects
  - Other environmental/planning layers
- This allows results to move from raw change detection to actionable restoration prioritization.

### Workflow Summary

```
Imagery (Date 1 & Date 2)
        ↓
Shoreline Extraction (NDWI or Brightness Threshold)
        ↓
Shoreline Cleanup/Smoothing
        ↓
50-m Transect Generation
        ↓
Transect–Shoreline Intersection
        ↓
Erosion/Accretion Distance & Rate (m/yr)
        ↓
Restoration Prioritization Overlay
```

## Requirements

- ArcGIS Pro (with `arcpy`)
- Two dates of aerial/satellite imagery covering the analysis area (e.g., NAIP)
- Spatial Analyst extension (for raster thresholding and Con expressions)

## Inputs

- Imagery for two time periods (multispectral with NIR, or high-res RGB)
- Approximate shoreline/lagoon edge (for buffer generation, if not already delineated)

## Outputs

- Extracted shoreline lines for each imagery date
- Transects with erosion/accretion classification
- Shoreline change rate (meters/year) per transect
- Optional overlay-ready layers for restoration planning

## Notes

- NDWI threshold and brightness threshold values may need empirical adjustment depending on imagery quality, water turbidity, and lighting conditions.
- Transect spacing and length are configurable; defaults reflect a 50-m spacing / 300-m length used in the current workflow.

## License

*(Add your preferred license here, e.g., MIT)*
