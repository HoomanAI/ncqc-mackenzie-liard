/**
 * gee_fig1a_sentinel2.js
 * ======================
 * Exports the Sentinel-2 true-colour basemap used by Figure 1a of the NCQC
 * manuscript: the Mackenzie-Liard confluence with all three gauges in frame.
 *
 * Paste into the Google Earth Engine Code Editor (https://code.earthengine.google.com),
 * press Run, then open the Tasks tab and click Run on the export task.
 *
 * !! DO NOT CHANGE `REGION` !!  04_make_figures.py georeferences the exported
 * image by assuming exactly this rectangle, so that no GIS library is needed
 * to place the gauge markers. If you must change it, change EXTENT in
 * 04_make_figures.py to match.
 */

// lon_min, lat_min, lon_max, lat_max  -- 56.3 km x 22.1 km at 61.8 N
var REGION = ee.Geometry.Rectangle([-121.62, 61.705, -120.55, 61.905], null, false);

// ---- cloud / shadow / snow mask from the L2A scene classification ----------
function maskS2(img) {
  var scl = img.select('SCL');
  var bad = scl.eq(3)                     // cloud shadow
    .or(scl.eq(8))                        // cloud, medium probability
    .or(scl.eq(9))                        // cloud, high probability
    .or(scl.eq(10))                       // thin cirrus
    .or(scl.eq(11));                      // snow / ice
  return img.updateMask(bad.not()).divide(10000);
}

// ---- open-water-season composite ------------------------------------------
// July-September only: the figure should show the channel ice-free, and the
// NCQC open-water regime is the one the planform is meant to illustrate.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(REGION)
  .filterDate('2021-07-01', '2025-09-30')
  .filter(ee.Filter.calendarRange(7, 9, 'month'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 25))
  .map(maskS2);

print('scenes in composite:', s2.size());

// visualize() bakes the stretch into an 8-bit 3-band image, so the export is a
// plain RGB raster that PIL can read without rasterio/GDAL.
var rgb = s2.select(['B4', 'B3', 'B2']).median()
            .visualize({min: 0.015, max: 0.17, gamma: 1.25});

Map.centerObject(REGION, 10);
Map.addLayer(rgb, {}, 'Sentinel-2 true colour');
Map.addLayer(ee.Image().paint(REGION, 0, 2), {palette: ['ffcc00']}, 'export region');

// ---- gauge locations, for a visual check before exporting ------------------
var gauges = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([-120.7917, 61.8165]), {code: '10FB006'}),  // Strong Point
  ee.Feature(ee.Geometry.Point([-121.2278, 61.7426]), {code: '10ED002'}),  // Liard
  ee.Feature(ee.Geometry.Point([-121.3589, 61.8685]), {code: '10GC001'})   // Fort Simpson
]);
Map.addLayer(gauges.style({color: 'ff0000', pointSize: 7}), {}, 'gauges');

// ---- export ----------------------------------------------------------------
// EPSG:4326 keeps the raster in plate carree, so pixel -> lon/lat is linear and
// the figure code can place markers from REGION alone.
Export.image.toDrive({
  image: rgb,
  description: 'MackenzieLiard_S2_truecolour',
  fileNamePrefix: 'MackenzieLiard_S2_truecolour',
  region: REGION,
  scale: 10,              // 10 m -> about 5630 x 2210 px; use 20 if the export is slow
  crs: 'EPSG:4326',
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});

// Attribution required in the figure caption:
//   "Contains modified Copernicus Sentinel data 2021-2025."
