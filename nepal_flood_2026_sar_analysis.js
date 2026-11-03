/*
==========================================================================================
  Multi-Temporal SAR Flood Mapping — August 2026 Nepal Flash Floods
  Bhotekoshi–Trishuli–Narayani corridor (Rasuwa → Nuwakot → Chitwan)

  Event: 26 Aug 2026 flash flood along Bhotekoshi river from Nepal-Tibet border
  Affected districts: Rasuwa (worst), Nuwakot, Dhading
  Deaths: 165+  |  Missing: 800+

  Author: Aditya | IIT Bombay | CSRE
==========================================================================================
*/

// ---------------------------------------------------------------------------------------
// 1. STUDY AREAS
//    Two ROIs: (A) Mountain corridor where flood originated
//              (B) Downstream Terai plains where flood spread and is easier to detect
//    SAR flood detection works best on flat terrain; the Terai ROI gives cleaner results
// ---------------------------------------------------------------------------------------

// ROI-A: Mountain corridor — Rasuwa, Nuwakot, Dhading (Bhotekoshi–Trishuli valley)
var roiMountain = ee.Geometry.Rectangle([85.05, 27.75, 85.55, 28.25]);

// ROI-B: Downstream Terai — Chitwan, Nawalparasi (Narayani river floodplain)
var roiTerai = ee.Geometry.Rectangle([83.90, 27.35, 84.85, 27.75]);

// Combined for overview
var roiFull = roiMountain.bounds().union(roiTerai.bounds(), 100);

Map.centerObject(roiFull, 9);
Map.addLayer(roiMountain, {color: 'orange'}, 'ROI: Mountain Corridor', true, 0.3);
Map.addLayer(roiTerai, {color: 'cyan'}, 'ROI: Terai Plains', true, 0.3);

// ---------------------------------------------------------------------------------------
// 2. TEMPORAL WINDOWS
//    Flash flood struck 26 Aug 2026 (~9 AM)
//    S1A + S1C operational → ~6-day revisit in 2026
// ---------------------------------------------------------------------------------------
var windows = {
  baseline:  {start: '2026-03-01', end: '2026-06-15'},   // dry season reference
  preEvent:  {start: '2026-08-14', end: '2026-08-26'},   // just before the flood
  peakFlood: {start: '2026-08-26', end: '2026-09-08'},   // flood event + immediate aftermath
  receding:  {start: '2026-09-08', end: '2026-09-22'}    // post-flood recession
};

// ---------------------------------------------------------------------------------------
// 3. SENTINEL-1 COLLECTION HELPER
// ---------------------------------------------------------------------------------------
function getS1Mosaic(start, end, roi) {
  var collection = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(roi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV', 'VH']);

  print('Scenes (' + start + ' to ' + end + '):', collection.size());
  return collection.mosaic().clip(roi);
}

// ---------------------------------------------------------------------------------------
// 4. REFINED LEE SPECKLE FILTER
// ---------------------------------------------------------------------------------------
function toNatural(img) { return ee.Image(10.0).pow(img.divide(10.0)); }

function refinedLeeFilter(image) {
  var bandNames = image.bandNames();

  var filterOneBand = function (bandName) {
    var img = image.select([bandName]);

    var weights3 = ee.List.repeat(ee.List.repeat(1, 3), 3);
    var kernel3 = ee.Kernel.fixed(3, 3, weights3, 1, 1, false);
    var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
    var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

    var sampleWeights = ee.List([
      [0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0],
      [0, 1, 0, 1, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1, 0],
      [0, 0, 0, 0, 0, 0, 0]
    ]);
    var sampleKernel = ee.Kernel.fixed(7, 7, sampleWeights, 3, 3, false);
    var sampleMean = mean3.neighborhoodToBands(sampleKernel);
    var sampleVar  = variance3.neighborhoodToBands(sampleKernel);

    var gradients = sampleMean.select(1).subtract(sampleMean.select(7)).abs();
    gradients = gradients.addBands(sampleMean.select(6).subtract(sampleMean.select(2)).abs());
    gradients = gradients.addBands(sampleMean.select(3).subtract(sampleMean.select(5)).abs());
    gradients = gradients.addBands(sampleMean.select(0).subtract(sampleMean.select(8)).abs());
    var maxGradient = gradients.reduce(ee.Reducer.max());
    var gradMask = gradients.eq(maxGradient);
    gradMask = gradMask.addBands(gradMask);

    var directions = sampleMean.select(1).subtract(sampleMean.select(4))
      .gt(sampleMean.select(4).subtract(sampleMean.select(7))).multiply(1);
    directions = directions.addBands(sampleMean.select(6).subtract(sampleMean.select(4))
      .gt(sampleMean.select(4).subtract(sampleMean.select(2))).multiply(2));
    directions = directions.addBands(sampleMean.select(3).subtract(sampleMean.select(4))
      .gt(sampleMean.select(4).subtract(sampleMean.select(5))).multiply(3));
    directions = directions.addBands(sampleMean.select(0).subtract(sampleMean.select(4))
      .gt(sampleMean.select(4).subtract(sampleMean.select(8))).multiply(4));
    directions = directions.addBands(directions.select(0).not().multiply(5));
    directions = directions.addBands(directions.select(1).not().multiply(6));
    directions = directions.addBands(directions.select(2).not().multiply(7));
    directions = directions.addBands(directions.select(3).not().multiply(8));
    directions = directions.updateMask(gradMask).reduce(ee.Reducer.sum());

    var sampleStats = sampleVar.divide(sampleMean.multiply(sampleMean));
    var sigmaV = sampleStats.toArray().arraySort().arraySlice(0, 0, 5)
      .arrayReduce(ee.Reducer.mean(), [0]);

    var rectWeights = ee.List.repeat(ee.List.repeat(0, 7), 3)
      .cat(ee.List.repeat(ee.List.repeat(1, 7), 4));
    var diagWeights = ee.List([
      [1, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0],
      [1, 1, 1, 1, 0, 0, 0], [1, 1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1, 0],
      [1, 1, 1, 1, 1, 1, 1]
    ]);
    var rectKernel = ee.Kernel.fixed(7, 7, rectWeights, 3, 3, false);
    var diagKernel = ee.Kernel.fixed(7, 7, diagWeights, 3, 3, false);

    var dirMean = img.reduceNeighborhood(ee.Reducer.mean(), rectKernel).updateMask(directions.eq(1));
    var dirVar  = img.reduceNeighborhood(ee.Reducer.variance(), rectKernel).updateMask(directions.eq(1));
    dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diagKernel).updateMask(directions.eq(2)));
    dirVar  = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diagKernel).updateMask(directions.eq(2)));

    for (var i = 1; i < 4; i++) {
      dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rectKernel.rotate(i)).updateMask(directions.eq(2 * i + 1)));
      dirVar  = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rectKernel.rotate(i)).updateMask(directions.eq(2 * i + 1)));
      dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diagKernel.rotate(i)).updateMask(directions.eq(2 * i + 2)));
      dirVar  = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diagKernel.rotate(i)).updateMask(directions.eq(2 * i + 2)));
    }

    dirMean = dirMean.reduce(ee.Reducer.sum());
    dirVar  = dirVar.reduce(ee.Reducer.sum());
    var varX = dirVar.subtract(dirMean.multiply(dirMean).multiply(sigmaV)).divide(sigmaV.add(1.0));
    var b = varX.divide(dirVar);
    var result = dirMean.add(b.multiply(img.subtract(dirMean)));

    return result.arrayFlatten([['sum']]).rename([bandName]);
  };

  var filtered = bandNames.map(filterOneBand);
  return ee.ImageCollection(filtered).toBands().rename(bandNames);
}

// ---------------------------------------------------------------------------------------
// 5. FLOOD DETECTION
// ---------------------------------------------------------------------------------------
var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
var srtm = ee.Image('USGS/SRTMGL1_003');

function detectFlood(baselineFiltered, waveFiltered, roi, slopeThreshold) {
  // INCREASE detection: double-bounce in urban/vegetated flooded areas
  var ratioVV = waveFiltered.select('VV').divide(baselineFiltered.select('VV'));
  var ratioVH = waveFiltered.select('VH').divide(baselineFiltered.select('VH'));
  var floodIncrease = ratioVV.gt(1.25).and(ratioVH.gt(1.30));

  // DECREASE detection: specular reflection on open water = backscatter drops
  var floodDecrease = ratioVV.lt(0.55).and(ratioVH.lt(0.55));

  // Combined
  var floodedRaw = floodIncrease.or(floodDecrease);

  // Mask permanent water
  var permanentWaterMask = gsw.select('seasonality').gte(5).unmask(0).not();

  // Slope mask — parameter allows different thresholds for mountain vs Terai
  var slopeMask = ee.Terrain.slope(srtm).lt(slopeThreshold);

  var floodedMasked = floodedRaw.updateMask(permanentWaterMask).updateMask(slopeMask);

  // Sieve: minimum 8 connected flood pixels
  var connections = floodedMasked.connectedPixelCount(25, true);
  return floodedMasked.updateMask(connections.gte(8)).selfMask().clip(roi);
}

// ---------------------------------------------------------------------------------------
// 6. RUN ANALYSIS FOR BOTH ROIs
// ---------------------------------------------------------------------------------------
var analysisConfigs = [
  {name: 'Mountain', roi: roiMountain, slopeThresh: 15},  // relaxed for valleys
  {name: 'Terai',    roi: roiTerai,    slopeThresh: 5}    // strict for plains
];

var wavePhases = [
  {key: 'preEvent',  label: 'Pre-Event (14-26 Aug)',  color: '#ffcc00', visible: false},
  {key: 'peakFlood', label: 'Peak Flood (26 Aug-8 Sep)', color: '#cc0000', visible: true},
  {key: 'receding',  label: 'Receding (8-22 Sep)',    color: '#3388ff', visible: false}
];

analysisConfigs.forEach(function (config) {
  print('========== ' + config.name + ' Region ==========');

  var baseRaw = getS1Mosaic(windows.baseline.start, windows.baseline.end, config.roi);
  var baseFiltered = refinedLeeFilter(toNatural(baseRaw));

  wavePhases.forEach(function (phase) {
    var w = windows[phase.key];
    var waveRaw = getS1Mosaic(w.start, w.end, config.roi);
    var waveFiltered = refinedLeeFilter(toNatural(waveRaw));

    var flood = detectFlood(baseFiltered, waveFiltered, config.roi, config.slopeThresh);

    var layerName = config.name + ' — ' + phase.label;
    Map.addLayer(flood, {palette: [phase.color]}, layerName,
                 phase.visible && config.name === 'Terai');  // Terai peak visible by default

    var areaKm2 = flood
      .multiply(ee.Image.pixelArea().divide(1e6))
      .reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry: config.roi,
        scale: 10,
        maxPixels: 1e10,
        tileScale: 4
      });

    print(config.name + ' | ' + phase.label + ' | flooded area (km²):', areaKm2);
  });
});

// ---------------------------------------------------------------------------------------
// 7. REFERENCE LAYERS
// ---------------------------------------------------------------------------------------

// Permanent water
var permWater = gsw.select('seasonality').gte(8).selfMask();
Map.addLayer(permWater, {palette: ['#0044cc']}, 'Permanent Water (JRC)', false);

// Elevation
Map.addLayer(srtm.clip(roiFull),
  {min: 100, max: 4000, palette: ['#1a9850', '#91cf60', '#d9ef8b', '#fee08b', '#fc8d59', '#d73027']},
  'Elevation (SRTM)', false);

// Slope
Map.addLayer(ee.Terrain.slope(srtm).clip(roiFull),
  {min: 0, max: 45, palette: ['white', 'yellow', 'orange', 'red']},
  'Slope (degrees)', false);

// ---------------------------------------------------------------------------------------
// 8. RAW BACKSCATTER COMPARISON (visual inspection)
// ---------------------------------------------------------------------------------------
var baseVV = getS1Mosaic(windows.baseline.start, windows.baseline.end, roiTerai).select('VV');
var floodVV = getS1Mosaic(windows.peakFlood.start, windows.peakFlood.end, roiTerai).select('VV');

Map.addLayer(baseVV, {min: -25, max: 0}, 'Terai Baseline VV (dB)', false);
Map.addLayer(floodVV, {min: -25, max: 0}, 'Terai Flood VV (dB)', false);

// Difference image for visual change detection
var diffVV = floodVV.subtract(baseVV);
Map.addLayer(diffVV, {min: -10, max: 10, palette: ['blue', 'white', 'red']},
  'VV Difference (flood - baseline)', false);

// ---------------------------------------------------------------------------------------
print('=== Nepal Flash Flood SAR Analysis Complete ===');
print('Event: 26 August 2026 Bhotekoshi flash flood');
print('Affected: Rasuwa, Nuwakot, Dhading → downstream Chitwan/Nawalparasi');
print('Data: Sentinel-1 GRD (S1A + S1C), IW mode, VV+VH');
print('Note: Mountain results have more noise due to SAR layover/shadow in steep terrain');
print('      Terai results are cleaner — flat terrain is ideal for SAR flood detection');

// ---------------------------------------------------------------------------------------
// 9. EXPORT RESULTS TO GOOGLE DRIVE
//    Run the script, then go to Tasks tab → click "Run" next to each export task
//    All GeoTIFFs will appear in Google Drive folder "Nepal_Flood_SAR"
// ---------------------------------------------------------------------------------------

// Recompute flood layers for export
var baseTerai_exp = refinedLeeFilter(toNatural(getS1Mosaic(windows.baseline.start, windows.baseline.end, roiTerai)));
var floodTerai_exp = refinedLeeFilter(toNatural(getS1Mosaic(windows.peakFlood.start, windows.peakFlood.end, roiTerai)));
var floodMapTerai = detectFlood(baseTerai_exp, floodTerai_exp, roiTerai, 5);

var baseMtn_exp = refinedLeeFilter(toNatural(getS1Mosaic(windows.baseline.start, windows.baseline.end, roiMountain)));
var floodMtn_exp = refinedLeeFilter(toNatural(getS1Mosaic(windows.peakFlood.start, windows.peakFlood.end, roiMountain)));
var floodMapMountain = detectFlood(baseMtn_exp, floodMtn_exp, roiMountain, 15);

// 1. Terai Flood Map (main result)
Export.image.toDrive({
  image: floodMapTerai.toFloat(),
  description: 'Nepal_Flood_Terai_PeakFlood',
  folder: 'Nepal_Flood_SAR',
  region: roiTerai,
  scale: 10,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});

// 2. Mountain Flood Map
Export.image.toDrive({
  image: floodMapMountain.toFloat(),
  description: 'Nepal_Flood_Mountain_PeakFlood',
  folder: 'Nepal_Flood_SAR',
  region: roiMountain,
  scale: 10,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});

// 3. Baseline VV backscatter (Terai) — before/after comparison
Export.image.toDrive({
  image: getS1Mosaic(windows.baseline.start, windows.baseline.end, roiTerai).select('VV'),
  description: 'Nepal_Baseline_VV_Terai',
  folder: 'Nepal_Flood_SAR',
  region: roiTerai,
  scale: 10,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});

// 4. Flood-period VV backscatter (Terai) — before/after comparison
Export.image.toDrive({
  image: getS1Mosaic(windows.peakFlood.start, windows.peakFlood.end, roiTerai).select('VV'),
  description: 'Nepal_Flood_VV_Terai',
  folder: 'Nepal_Flood_SAR',
  region: roiTerai,
  scale: 10,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});

// 5. VV Difference image (Terai) — change detection visualization
var baseVV_exp = getS1Mosaic(windows.baseline.start, windows.baseline.end, roiTerai).select('VV');
var floodVV_exp = getS1Mosaic(windows.peakFlood.start, windows.peakFlood.end, roiTerai).select('VV');
Export.image.toDrive({
  image: floodVV_exp.subtract(baseVV_exp).rename('VV_difference'),
  description: 'Nepal_VV_Difference_Terai',
  folder: 'Nepal_Flood_SAR',
  region: roiTerai,
  scale: 10,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF'
});
