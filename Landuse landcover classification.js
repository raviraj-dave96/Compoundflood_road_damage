// -----------------------------------------------------------------------
// Sentinel-1 and Sentinel-2 Data Preprocessing and Classification
// -----------------------------------------------------------------------

// Define Region of Interest (ROI) — change 'roi' to your defined area
var roi = ee.Geometry.Polygon([[[longitude1, latitude1], [longitude2, latitude2], ...]]); // Define your area here

// -----------------------------------------------------------------------
// Sentinel-1 C-band SAR Ground Range Collection (Log scale, VV, Descending)
// -----------------------------------------------------------------------
var collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filterMetadata('resolution_meters', 'equals', 10)
  .filterBounds(roi)
  .select('VV');
print(collectionVV, 'Sentinel-1 VV Collection');

// -----------------------------------------------------------------------
// Sentinel-1 C-band SAR Ground Range Collection (Log scale, VH, Descending)
// -----------------------------------------------------------------------
var collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filterMetadata('resolution_meters', 'equals', 10)
  .filterBounds(roi)
  .select('VH');
print(collectionVH, 'Sentinel-1 VH Collection');

// -----------------------------------------------------------------------
// Filter by Date Range and Create Mosaic for Both VV and VH
// -----------------------------------------------------------------------
var SARVV = collectionVV.filterDate('2019-01-01', '2019-01-10').mosaic();
var SARVH = collectionVH.filterDate('2019-01-01', '2019-01-10').mosaic();

// -----------------------------------------------------------------------
// Cloud Masking Functions for Sentinel-2 and Landsat
// -----------------------------------------------------------------------
function maskL8sr(image) {
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var qa = image.select('pixel_qa');
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0).and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask).divide(10000).select("B[0-9]*").copyProperties(image, ["system:time_start"]);
}

function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

// -----------------------------------------------------------------------
// Sentinel-2 Data Preprocessing and Cloud Masking
// -----------------------------------------------------------------------
var dataset = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate('2019-01-01', '2019-01-10')
  .filterBounds(roi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2clouds);
print(dataset, 'Sentinel-2 Data');

// -----------------------------------------------------------------------
// Calculate NDVI and Create Composite Image
// -----------------------------------------------------------------------
var comp = dataset.mean();
var ndvi = comp.normalizedDifference(['B8', 'B4']).rename('NDVI');
var composite = ee.Image.cat(comp, ndvi);

// -----------------------------------------------------------------------
// Apply Smoothing Filter to SAR Data
// -----------------------------------------------------------------------
var SMOOTHING_RADIUS = 30;
var SARVV_filtered = SARVV.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');
var SARVH_filtered = SARVH.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');

// -----------------------------------------------------------------------
// Merge Feature Collections (Define the classes for land cover)
// -----------------------------------------------------------------------
var newfc = open_water.merge(baren_land).merge(Vegetation1).merge(urban).merge(Road); // Merge your feature collections here

// -----------------------------------------------------------------------
// SAR Classification
// -----------------------------------------------------------------------
var final = ee.Image.cat(SARVV_filtered, SARVH_filtered);
var bands = ['VH', 'VV'];
var training = final.select(bands).sampleRegions({
  collection: newfc,
  properties: ['landcover'],
  scale: 10
});

// Train the classifier (Random Forest)
var classifier = ee.Classifier.smileRandomForest(20).train({
  features: training,
  classProperty: 'landcover',
  inputProperties: bands
});

// Run Classification
var classified = final.select(bands).classify(classifier);
print('RF-SAR error matrix: ', classifier.confusionMatrix());
print('RF-SAR accuracy: ', classifier.confusionMatrix().accuracy());

// -----------------------------------------------------------------------
// Optical Classification
// -----------------------------------------------------------------------
var bandss2 = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B11', 'B12', 'NDVI'];
var trainings2 = composite.select(bandss2).sampleRegions({
  collection: newfc,
  properties: ['landcover'],
  scale: 10
});

var classifiers2 = ee.Classifier.smileRandomForest(20).train({
  features: trainings2,
  classProperty: 'landcover',
  inputProperties: bandss2
});

var classifieds2 = composite.select(bandss2).classify(classifiers2);
print('RF-S2 training error matrix: ', classifiers2.confusionMatrix());
print('RF-S2 training accuracy: ', classifiers2.confusionMatrix().accuracy());

// -----------------------------------------------------------------------
// Validation for Optical Classification
// -----------------------------------------------------------------------
var valfc = ee.FeatureCollection(test_sample1); // Replace with your validation sample
var validationS2 = composite.select(bandss2).sampleRegions({
  collection: valfc,
  properties: ['Code'],
  scale: 10
});

validationS2 = validationS2.classify(classifiers2);
var validationAccuracyS2 = validationS2.errorMatrix('Code', 'classification');
print('Validation optical error matrix', validationAccuracyS2);
print('Validation optical accuracy', validationAccuracyS2.accuracy());

// -----------------------------------------------------------------------
// Combined Optical + SAR Classification
// -----------------------------------------------------------------------
var opt_sar = ee.Image.cat(composite, SARVV_filtered, SARVH_filtered);
var bands_opt_sar = ['VH', 'VV', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B11', 'B12', 'NDVI'];
var training_opt_sar = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: newfc,
  properties: ['landcover'],
  scale: 10
});

var classifier_opt_sar = ee.Classifier.smileRandomForest(20).train({
  features: training_opt_sar,
  classProperty: 'landcover',
  inputProperties: bands_opt_sar
});

var classifiedboth = opt_sar.select(bands_opt_sar).classify(classifier_opt_sar);
print('RF-Opt/SAR error matrix: ', classifier_opt_sar.confusionMatrix());
print('RF-Opt/SAR accuracy: ', classifier_opt_sar.confusionMatrix().accuracy());

// -----------------------------------------------------------------------
// Validation for Optical + SAR Classification
// -----------------------------------------------------------------------
var validation_opt_sar = opt_sar.select(bands_opt_sar).sampleRegions({
  collection: valfc,
  properties: ['Code'],
  scale: 10
});

validation_opt_sar = validation_opt_sar.classify(classifier_opt_sar);
var validationAccuracyOS = validation_opt_sar.errorMatrix('Code', 'classification');
print('Validation optical/SAR error matrix', validationAccuracyOS);
print('Validation optical/SAR accuracy', validationAccuracyOS.accuracy());
print("Consumer's optical/SAR accuracy", validationAccuracyOS.consumersAccuracy());
print("Producer's accuracy", validationAccuracyOS.producersAccuracy());
print('Validation optical/SAR Kappa coefficient: ', validationAccuracyOS.kappa());

// -----------------------------------------------------------------------
// Final Classification and Export
// -----------------------------------------------------------------------
var aoi = ee.FeatureCollection(roi); // Define your AOI
Map.addLayer(aoi, {}, 'ROI');
Map.centerObject(aoi, 8);

// Clip classification to AOI
var LUData_clip = classifiedboth.clip(aoi);

// Display the Classification
Map.addLayer(LUData_clip, {min: 1, max: 5, palette: ['1667fa', 'c9270d', 'cf7b68', 'ee9a1c', '146d0e']}, 'Optical/SAR Classification_clipped');

// Export the Classification Result
Export.image.toDrive({
  image: LUData_clip, 
  description: 'Optical_Radar2',
  fileNamePrefix: 'LUclass2',
  scale: 10,
  region: aoi, 
  maxPixels: 1e10
});
