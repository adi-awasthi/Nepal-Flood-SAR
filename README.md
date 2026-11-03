# Multi-Temporal SAR Flood Mapping — August 2026 Nepal Flash Floods

> **Bhotekoshi-Trishuli-Narayani Corridor | Google Earth Engine | Sentinel-1 SAR**

A Google Earth Engine (GEE) script for detecting and mapping flood inundation during the August 2026 Nepal flash floods using Sentinel-1 Synthetic Aperture Radar (SAR) change detection. The project maps flood progression across **two distinct terrain types** — mountain valleys and downstream Terai plains — with adaptive thresholds for each.

---

## Event Background

| Detail | Info |
|---|---|
| **Event** | Flash flood along Bhotekoshi river from Nepal-Tibet border |
| **Date** | 26 August 2026 (~9 AM local) |
| **Cause** | Glacial lake outburst / extreme rainfall → Bhotekoshi river surge |
| **Affected Districts** | Rasuwa (worst), Nuwakot, Dhading, downstream Chitwan & Nawalparasi |
| **Casualties** | 165+ dead, 800+ missing |
| **Flood Path** | Bhotekoshi → Trishuli → Narayani river system |

---

## Study Areas

Two separate ROIs account for the dramatically different terrain:

| ROI | Region | Bounding Box | Terrain | Slope Threshold |
|---|---|---|---|---|
| **Mountain Corridor** | Rasuwa, Nuwakot, Dhading | `[85.05, 27.75, 85.55, 28.25]` | Steep valleys, 1500-4000m | 15 degrees |
| **Terai Plains** | Chitwan, Nawalparasi | `[83.90, 27.35, 84.85, 27.75]` | Flat floodplain, 100-300m | 5 degrees |

> SAR flood detection works best on flat terrain. The Terai ROI produces cleaner results; the Mountain ROI captures the flood origin but has more noise from SAR geometric distortions (layover, shadow).

---

## Temporal Windows

| Phase | Start | End | Description |
|---|---|---|---|
| **Baseline** | 2026-03-01 | 2026-06-15 | Dry season reference (pre-monsoon) |
| **Pre-Event** | 2026-08-14 | 2026-08-26 | Just before the flood |
| **Peak Flood** | 2026-08-26 | 2026-09-08 | Flood event + immediate aftermath |
| **Receding** | 2026-09-08 | 2026-09-22 | Post-flood recession |

---

## Data Sources

| Dataset | Description | Use |
|---|---|---|
| `COPERNICUS/S1_GRD` | Sentinel-1 Ground Range Detected SAR | Primary flood detection |
| `JRC/GSW1_4/GlobalSurfaceWater` | JRC Global Surface Water | Permanent water masking |
| `USGS/SRTMGL1_003` | SRTM 30m Digital Elevation Model | Slope masking & terrain visualization |

### Sentinel-1 Parameters
- **Instrument Mode:** IW (Interferometric Wide swath)
- **Orbit Pass:** Both ascending and descending (maximizes scene count)
- **Polarization:** VV + VH (dual-pol)
- **Resolution:** 10 m
- **Satellites:** S1A + S1C (operational in 2026, ~6-day revisit)

---

## Methodology

```
Sentinel-1 GRD (VV + VH)
      |
      v
  Temporal Mosaicking (per window, per ROI)
      |
      v
  Convert dB --> Linear Power (natural units)
      |
      v
  Refined Lee Speckle Filter (7x7, directional, edge-preserving)
      |
      v
  Change Detection (flood / baseline backscatter ratio)
   |-- INCREASE: VV ratio > 1.25 AND VH ratio > 1.30
   |   (double-bounce scattering in flooded vegetation/urban)
   |
   |-- DECREASE: VV ratio < 0.55 AND VH ratio < 0.55
       (specular reflection on open water = backscatter drops)
      |
      v
  Post-Classification Masks
   |-- Permanent Water Mask (JRC GSW seasonality >= 5 months)
   |-- Slope Mask (adaptive: 15 deg mountain / 5 deg Terai)
      |
      v
  Connected Pixel Filter (sieve: min 8 connected pixels)
      |
      v
  Flood Inundation Map + Area Estimation (km^2)
```

### Key Improvements Over Standard Approach

1. **Dual ROI with Adaptive Thresholds** — Different slope thresholds for mountain (15 deg) vs. flat terrain (5 deg), acknowledging that steep valleys need relaxed filtering while plains need strict filtering.

2. **Bidirectional Change Detection** — Detects both backscatter *increase* (double-bounce in flooded vegetation/urban areas) and *decrease* (specular reflection on open water), capturing the full range of flood signatures.

3. **Multi-Temporal Progression** — Four temporal windows track the flood from pre-event through peak to recession, enabling temporal analysis of flood dynamics.

4. **No Orbit-Pass Filter** — Uses both ascending and descending passes to maximize scene count, critical for rapid-response flood mapping.

---

## File Structure

```
Nepal-Flood-SAR/
|
|-- nepal_flood_2026_sar_analysis.js   # Main GEE script (analysis + export)
|-- README.md                          # Project documentation
|-- results/                           # (optional) screenshots and exported maps
    |-- terai_flood_map.png
    |-- mountain_flood_map.png
    |-- console_output.png
    |-- vv_difference.png
```

---

## How to Run

1. Open [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
2. Create a new repository and file, paste the contents of `nepal_flood_2026_sar_analysis.js`
3. Click **Run**
4. Check the **Console** tab for:
   - Scene counts per temporal window (verify data availability)
   - Flooded area estimates in km^2
5. Toggle map layers on/off to inspect:
   - Flood maps (red = peak flood, yellow = pre-event, blue = receding)
   - Reference layers (elevation, slope, permanent water, VV backscatter)
6. Go to **Tasks** tab and click **Run** next to each export task to download GeoTIFFs to Google Drive

### Map Layers

| Layer | Color | Default |
|---|---|---|
| Mountain - Pre-Event | Yellow `#ffcc00` | Hidden |
| Mountain - Peak Flood | Red `#cc0000` | Hidden |
| Mountain - Receding | Blue `#3388ff` | Hidden |
| **Terai - Peak Flood** | **Red `#cc0000`** | **Visible** |
| Terai - Pre-Event | Yellow `#ffcc00` | Hidden |
| Terai - Receding | Blue `#3388ff` | Hidden |
| Permanent Water (JRC) | Dark blue | Hidden |
| Elevation (SRTM) | Green-Red gradient | Hidden |
| Slope | White-Red gradient | Hidden |
| Baseline VV (dB) | Grayscale | Hidden |
| Flood VV (dB) | Grayscale | Hidden |
| VV Difference | Blue-White-Red | Hidden |

### Exported GeoTIFFs (Google Drive)

All exports go to folder `Nepal_Flood_SAR/` at 10m resolution:

| File | Description |
|---|---|
| `Nepal_Flood_Terai_PeakFlood.tif` | Binary flood map - Terai plains |
| `Nepal_Flood_Mountain_PeakFlood.tif` | Binary flood map - Mountain corridor |
| `Nepal_Baseline_VV_Terai.tif` | Pre-flood VV backscatter |
| `Nepal_Flood_VV_Terai.tif` | Flood-period VV backscatter |
| `Nepal_VV_Difference_Terai.tif` | VV change (flood minus baseline) |

---

## Technical Notes

- **Mountain terrain challenges:** SAR imaging in steep terrain produces geometric distortions — layover (slopes facing the sensor appear compressed), foreshortening, and radar shadow (slopes facing away get no signal). These cause false positives in flood detection. The relaxed 15-degree slope mask helps, but mountain results should be interpreted with caution.

- **Specular vs. double-bounce:** Open water acts as a mirror (specular reflection), sending the SAR signal away from the sensor = lower backscatter. Water under vegetation or between buildings creates double-bounce = higher backscatter. Detecting both is essential for complete flood mapping.

- **Permanent water masking:** JRC Global Surface Water `seasonality >= 5` removes pixels that are wet for 5+ months/year (rivers, lakes), so only *new* flooding is mapped.

- **Connected pixel filter:** Removes isolated 1-7 pixel clusters that are likely speckle noise rather than real flood extent.

---

## Dependencies

- [Google Earth Engine](https://earthengine.google.com/) account (free for non-commercial use)
- No external libraries — uses built-in GEE JavaScript API

---

## References

- Sentinel-1 SAR data: [ESA Copernicus](https://sentinel.esa.int/web/sentinel/missions/sentinel-1)
- JRC Global Surface Water: [Pekel et al., 2016 - Nature](https://www.nature.com/articles/nature20584)
- Refined Lee Filter: Quegan & Yu, 2001; adapted implementation for GEE
- SRTM DEM: [USGS / NASA](https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1)
- NASA ARSET SAR Training: [Applied Remote Sensing Training](https://appliedsciences.nasa.gov/get-involved/training/english/arset-sar-flood-mapping)

---

## Author

**Aditya**
M.Tech, Centre of Studies in Resources Engineering (CSRE)
IIT Bombay

---

## License

For academic and research use only. SAR data courtesy of ESA Copernicus Open Access Hub.
