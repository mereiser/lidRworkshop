# Introduction to lidR workshop 
# Ecological Restoration Institute 
# Day 2 - First Order Products 
# 8/23/23

###################################################################################
## We will process a catalog to build a DEM and pitless canopy height model 
## (CHM), as well as DSM, mean intensity rasters (as a rasterBrick, which is more 
## memory effecient). RasterStack and RasterBrick are very similar, the difference 
## being in the virtual characteristic of the RasterStack. While a RasterBrick 
## has to refer to one multi-layer file or is in itself a multi-layer object 
## with data loaded in memory, a RasterStack may "virtually" connect several 
## raster objects written to different files or in memory. Processing will be
## more efficient for a RasterBrick than for a RasterStack, but RasterStack 
## has the advantage of facilitating pixel based calculations on separate 
## raster layers.
###################################################################################

# Load libraries
library(lidR)
library(terra)        # 1.5-21 
library(sf)           # 1.0-7
library(leaflet)      # 2.1.1
library(parallel)
library(dplyr)
library(mapview)
library(raster)
library(rasterVis)

# Set threads 
plan(multisession, workers = (parallel::detectCores(logical=FALSE)-2))
set_lidr_threads(1L)

# Location where you unzipped the files
dir = "C:/Users/mr2988//lidR_workshop/"

# Build a project
project = catalog(paste0(dir, "input/las/"))
st_crs(project) = 26912

# visualize it 
plot(project, map=TRUE)

# If you want to use a subset of the catalog, you can intersect it with a ROI that has the same projection
# and this will remove the "extra" las tiles. It will not reshape the tiles though... see catalog_retile for that
roi = sf::read_sf(paste0(dir, "input/vector/", "FtValleyExROI_Subset.shp"))
roi = sf::st_transform(roi, sf::st_crs(project))
project = catalog_intersect(project, roi)

# Set some catalog options
opt_chunk_buffer(project) = 10        # I typically use the maximum expected radius for a huge tree, say ~15m 
opt_chunk_size(project) = 200         # This is an art, but I typically use between 150-300m. This will be a funciton of your computing environment
opt_select(project) = "xyzicr"        # we only need the coordinates, height, intensity and classification, and return number information
opt_stop_early(project) = FALSE       # Use at your peril
opt_wall_to_wall(project) = TRUE
opt_chunk_alignment(project) = c(0,0)
opt_progress(project) = TRUE


# Run some diagnostics of the catalog
las_check(project) # Performs a deep inspection of the catalog and prints a report
summary(project) # Prints a summary overview of the catalog

# the las/laz aren't indexed using lax, so have a look at help("writelax")
# if you want to index all the las/laz files in a folder you can use
# files=list.files(paste0(dir, "input/las/"), full.names = TRUE)
# for (i in files) writelax(i)
# for more info, see https://cran.r-project.org/web/packages/lidR/vignettes/lidR-computation-speed-LAScatalog.html

# Interpolate the ground classified points and creates a digital elevation model (DEM) raster
dtm = rasterize_terrain(project, res=1, algorithm = knnidw(k = 8L, p = 2))

# Plot it
plot(dtm, col=height.colors(100))
plot(terra::mask(dtm, vect(roi)), col=height.colors(100))

# We can also plot the hillshade. When we do this from a catalog we get a seamless image 
# compared to stitching together las tiles where we get lines between tiles 
slope = terra::terrain(dtm, "slope", unit="radians")
aspect = terra::terrain(dtm, "aspect", unit="radians")
hillshade = terra::shade(slope, aspect, 40, 270)
plot(hillshade, col = gray(0:50/50), legend = FALSE)

# Now let's construct a pit-free digital surface model (DSM) based on the computation of a set of
# triangulations at different heights (Khosravipour et al. 2014). Other methods exist.
# Also, this approach assumes you have no outliers or noise issues to resolve. This dataset
# has both 
dsm = rasterize_canopy(project, res = 1, pitfree(c(0,2,5,10,15), c(0, 1.5)))
plot(terra::mask(dsm, vect(roi)), col=height.colors(100))
lidR::plot_dtm3d(dsm)

# We'll just bypass the creation of a DSM and create a pit-free canopy height model (CHM), also 1m and we'll 
# do our best to remove the building and powerline points prior to processing. FYI - our best isn't that good!
# My data contains outliers so we have to remove them (filter) before we
# build the chm and thus we need a outlier/noise removal function. Here we do it by normalizing
# the point cloud and removing point above a threshold (50m in this case) and then creating the CHM. Pretty simple...
processed_chm = function(cluster, threshold)
{
  las = readLAS(cluster)
  if (is.empty(las)) return(NULL) # No need to process the las object if it's empty
  
  las = normalize_height(las, knnidw(k = 8, p = 2)) # Normalize it
  las = filter_poi(las, Z < threshold) # get rid of points above 50m hieght above the ground (i.e., ain't no tree northern AZ tree > 164 feet!)
  
  if (is.null(las)) return(NULL)
  chm = rasterize_canopy(las, res=1, pitfree(thresholds = c(0, 2, 5, 10, 15), max_edge = c(0, 1.5), subcircle = 0.5))
  # More info - https://github.com/Jean-Romain/lidR/wiki/Rasterizing-perfect-canopy-height-models
  
  bbox = lidR::extent(cluster) # Define the extent of the original raster
  chm = terra::crop(chm, bbox) #and crop it
  return(chm)
}
chm = catalog_apply(project, processed_chm, threshold = 50)
chm = do.call(terra::merge, chm)
chm = terra::mask(chm, vect(as.spatial(project)))

# Plot it, multiple ways
plot(terra::mask(chm, vect(roi)), col=height.colors(100))
rasterVis::plot3D(as(chm, "Raster"), zfac=1, col=height.colors(100))

# We could just add the two rasters to get a pitless digtal surface model (DSM)  
dsm = terra::crop(dtm, chm) + chm
# Plot it
plot(dsm, col=height.colors(50))

# Set all values below dbh to NA and those above to 1 to create a "forest canopy" mask
canopy = terra::classify(chm, matrix(c(-Inf,1.37,NA,  1.37,Inf,1), ncol=3, byrow=TRUE)) 
# Check it, by masking the chm by he canopy mask, and then we will use it later
plot(terra::mask(chm,canopy), col=height.colors(50))
plot(terra::mask(terra::mask(chm,canopy),vect(roi)), col=height.colors(50))

# and for fun! Why do the edges look like that? Can you think of a way to avoid those edge effects?
rasterVis::plot3D(as(dtm, "Raster"), zfac=1, col=height.colors(100), drape=raster::raster(canopy))

# Build a mean intensity raster, removing the outliers 
filtered_int = function(cluster)
{
  las = readLAS(cluster)
  if (is.empty(las)) return(NULL)
  las = filter_poi(las, Classification != LASNOISE) # Filter out the potential noise two ways
  mean_int = pixel_metrics(las, ~mean(Intensity), res = 1)
  bbox  = lidR::extent(cluster)
  mean_int = terra::crop(mean_int, bbox)
  return(mean_int)
}
mean_int = catalog_apply(project, filtered_int)
mean_int = do.call(terra::merge, mean_int)

# Plot it
cuts = quantile(values(mean_int), probs = seq(0, 1, 0.1), na.rm = TRUE) #set breaks
pal = colorRampPalette(c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'))
plot(mean_int, breaks=cuts, col = rev(pal(10)))  #plot with defined breaks

# Create a RasterBrick object of the four layers
firstOrderProd = c(dtm, dsm, chm, mean_int)

# Write the brick to an GeoTiff this saves it as an ESRI world file 
writeRaster(firstOrderProd, filename = paste0(dir, "output/fstOrderProds.tif"), 
            gdal=c("COMPRESS=NONE", "TFW=YES","of=COG"), options=c('TFW=YES'), overwrite = TRUE)
writeRaster(chm, filename = paste0(dir, "output/CHM_1m.tif"), 
            gdal=c("COMPRESS=NONE", "TFW=YES","of=COG"), options=c('TFW=YES'), overwrite = TRUE)
writeRaster(dtm, filename = paste0(dir, "output/DEM_1m.tif"), 
            gdal=c("COMPRESS=NONE", "TFW=YES","of=COG"), options=c('TFW=YES'), overwrite = TRUE)
writeRaster(canopy, filename = paste0(dir, "output/CC_1m.tif"), 
            gdal=c("COMPRESS=NONE", "TFW=YES","of=COG"), options=c('TFW=YES'), overwrite = TRUE)

