// Region of interest (ROI) selection
var roi = ee.FeatureCollection(table); // Replace 'table' with your own feature collection variable
Map.addLayer(roi, {}, 'ROI');
Map.centerObject(roi, 8);

// Image selection and filtering from Sentinel-1 SAR imagery
var R1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi) // Filters the images based on the ROI
  .filterDate('2018-08-12', '2018-08-17') // Filters based on the specific date range
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')); // Filters based on VH polarization

// Display the first image of the collection
Map.addLayer(R1.first(), {bands: 'VH', min: -20, max: 0}, 'SAR image');

// Function for speckle filtering to reduce noise effect
var filter_Speckles = function(img) {
  var vh = img.select('VH'); // Select the VH polarization band
  var vh_smoothed = vh.focal_median(100, 'circle', 'meters').rename('VH_Filtered'); // Apply a focal median filter
  return img.addBands(vh_smoothed); // Add filtered VH band to the original image
};

// Apply speckle filter to the entire collection
R1 = R1.map(filter_Speckles);

// Display the first filtered image
Map.addLayer(R1.first(), {bands: 'VH_Filtered', min: -20, max: 0}, 'Filtered SAR image');

// Function to classify waterbody (flood) using Otsu's threshold
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
   return means.sort(bss).get([-1]);
};

// Function to classify water using Otsu's threshold
var classify_Water_Otsu = function(img) {
  var vh = img.select('VH_Filtered');
  // Create histogram for Otsu's method
  var histogram = vh.reduceRegion({
    reducer: ee.Reducer.histogram(255, 2).combine('mean', null, true).combine('variance', null, true), 
    geometry: roi,
    scale: 10,
    bestEffort: true
  });
  
  // Apply Otsu thresholding
  var threshold = otsu(histogram.get('VH_histogram'));
  
  // Classify water using threshold
  var water = vh.lt(threshold).rename('Water');
  water = water.updateMask(water); // Mask out non-water pixels
  return img.addBands(water); // Add water classification as a new band
};

// Apply water classification using Otsu's method to the entire image collection
R1 = R1.map(classify_Water_Otsu);
print(R1); // Print the image collection to the console for inspection

// Select the water classification band for the first image
var classification = ee.Image(R1.first()).clip(roi).select('Water');
var SARimage = ee.Image(R1.first()); // Select the first image of the collection

// Display SAR image with VH band
var R1_Layer = ui.Map.Layer(SARimage, {
  bands: ['VH'],
  max: 0,
  min: -20
});
Map.layers().reset([R1_Layer]);

// Visualization parameters for water classification
var visParams = {
  min: 0,
  max: 1,
  palette: ['#FFFFFF', '#0000FF'] // Water is displayed in blue
};

// Add the water classification layer to the map
Map.addLayer(classification, visParams, 'Water');

// Export flood extent as raster to Google Drive
Export.image.toDrive({
  image: classification, 
  description: 'Flood_extent_raster',
  fileNamePrefix: 'flooded',
  scale: 10,
  region: roi, // Specify region of interest
  maxPixels: 1e10 // Increase the pixel limit if necessary
});

// Export the flooded area as shapefile (for further analysis in QGIS or ArcGIS)
// Convert the flood raster to polygons
var classification_vec = classification.reduceToVectors({
  scale: 10,
  geometryType: 'polygon',
  geometry: roi,
  eightConnected: false, // Not considering diagonals for connectivity
  bestEffort: true, // Maximize the effort to achieve a good result
  tileScale: 2, // Adjust to manage large areas
});

// Export flood polygons as shapefile to Google Drive
Export.table.toDrive({
  collection: classification_vec,
  description: 'Flood_extent_vector',
  fileFormat: 'SHP',
  fileNamePrefix: 'flooded_vec'
});
