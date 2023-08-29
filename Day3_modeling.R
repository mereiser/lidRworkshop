#Modeling 
# I removed the voxel and tree metrics from the lidRattR package. That package offers additional metrics 
# for trees and voxels which can help improve the model. Highly reccommend checking it out if you're 
# looking to make wall-to-wall predictions with lidar 

#load packages
library(lidR)
library(terra)
library(sf)
library(mapview)
library(stats)
library(readxl)
library(dplyr)
library(corrplot)
library(randomForest) 
library(rfUtilities)  # remotes::install_github("jeffreyevans/rfUtilities", force = TRUE)
library(Boruta)       

dir = "C:/Users/mr2988/OneDrive - Northern Arizona University/lidar/lidR_workshop/"

# Set threads 
plan(multisession, workers = (parallel::detectCores(logical=FALSE)/2)-1)
set_lidr_threads(2L)

# Build a project
project = catalog(paste0(dir, "input/las/"))
crs(project) = sf::st_crs("epsg:26912")

# Set catalog options
opt_chunk_buffer(project) = 15        
opt_chunk_size(project) = 260         
opt_select(project) = "xyzicra"        
opt_stop_early(project) = FALSE       
opt_wall_to_wall(project) = TRUE
opt_progress(project) = TRUE
opt_filter(project) <- "-thin_with_voxel 0.1"

# Read in the shapefile of the inventory plot locations and buffer to the plot radius 
pointList = read_sf(paste0(dir, "input/vector/", "FUWI_Plots_WGS84_4326_Subset.shp"))
pointList = st_transform(pointList, sf::st_crs(project))
plotList = st_buffer(pointList, dist=11.34970)

# Clip the catalog to the plot boundaries and and normalize the point clouds 
plot_ctg <- clip_roi(project, plotList)
plot_ctg = lapply(plot_ctg, function(x) normalize_height(x, knnidw()))
plot_ctg = lapply(plot_ctg, function(x) filter_poi(x, Z > 0.1))
#plot(new_ctg[[1]])

# Calculate metrics for each plot '
plots_metrics <- lapply(plot_ctg, cloud_metrics, func = .stdmetrics)
plots_metrics <- data.table::rbindlist(plots_metrics)
plots = cbind(plotList, plots_metrics)

# Read in our FVS output so we can model biomass across the study area 
FVS_plots = read_xls(paste0(dir, "input/data/FVSOut.xls"), sheet = "FVS_Compute")
biomass = FVS_plots[c("StandID", "STBIOMSS")]
plots = inner_join(biomass, plots, join_by(StandID == Comment))
#plots = plots[-c(3:24)]
# Clean up our data so we can plot a correlation matrix to select our predictor variables 
plots_numeric = plots %>% select(where(is.numeric)) %>% select_if(~ !any(is.na(.)))
plots_numeric = plots_numeric[, colSums(plots_numeric != 0) > 0]
plots_numeric = plots_numeric[, sapply(plots_numeric, function(x) { sd(x) != 0} )]

# And now we can look at our correlation matrix. 
cor_mat = cor(as.matrix(plots_numeric))
corrplot(cor_mat)
# We'll also list out the 5 strongest correlations (positive and negative)
as.data.frame(cor_mat) %>% arrange(desc(STBIOMSS)) %>% select("STBIOMSS") %>% head(5)
as.data.frame(cor_mat) %>% arrange(STBIOMSS) %>% select("STBIOMSS") %>% head(5)

#Build a few linear models from a few highly correlated variables and pick the best one 
m1 = lm(STBIOMSS ~ itot, data = plots)
summary(m1)
m2 = lm(STBIOMSS ~ itot + n, data = plots)
summary(m2)
m3 = lm(STBIOMSS ~ itot + ipcumzq70, data = plots)
summary(m3)

#Plot our observed vs predicted biomass for m3
plot(plots$STBIOMSS, predict(m3))
abline(0,1)

#Now we can apply this plot-based model for wall-to-wall coverage across our study area
clean_las = function(cluster, threshold)
{
  las = readLAS(cluster)
  if (is.empty(las)) return(NULL) # No need to process the las object if it's empty
  
  las = normalize_height(las, knnidw()) # Normalize it
  las = filter_poi(las, Z < threshold) 
  las = filter_poi(las, Z > 0.1)
  if (is.null(las)) return(NULL)
  metrics = grid_metrics(las, .stdmetrics, res = 11.34970*2)
  #terra::crs(metrics) = "epsg:26912"
  bbox = lidR::extent(cluster)
  metrics = terra::crop(metrics, bbox)
  return(metrics)
}

metrics = catalog_apply(project, clean_las, threshold = 50)
metrics_w2w = do.call(terra::merge, metrics)
names(metrics_w2w) = names(metrics[[1]])

#Plot our predictor variables 
plot(metrics_w2w$itot)
plot(metrics_w2w$ipcumzq90)

#Make our predictions and plot!
STBIOMSS_pred = coef(m3)[1] + metrics_w2w$p2th*coef(m3)[2] + metrics_w2w$ipcumzq90*coef(m3)[3]
plot(STBIOMSS_pred)

################################################################################################
#Now let's kick it up a notch and try modeling without having to manually select our variables
#and as a random forest model instead of a linear one. We'll use the package boruta to do our 
#variable selection. We already have our project and plots loaded so we won't do that again

# Define a custom set of las metrics to extract for each plot (this is basically .stdmetrics, 
# repeated to give you an example of a custom metrics function)
source(paste0(dir, "scripts/", "myMetricsScript.R"))

#Custom metrics function 
myMetrics = function(x, y, z, i, rn, class, dz, th, minht, above)
{
  S  = stdmetrics(x,y,z,i,rn,class,dz,th)
  ShpM = stdshapemetrics(x,y,z)
  mSM = mystdMetrics(z,i,rn,minht,above)
  metrics = c(S,mSM,ShpM)
  return(metrics)
}

# calculate our custom myMetrics for each plot and store the results in a dataframe with the associated plot ID
overall_a_metrics = lapply(plot_ctg, function(x) cloud_metrics(x, myMetrics(X, Y, Z, Intensity, ReturnNumber, Classification, 
                                                                           dz = 1, th = 2, minht = 1.37, above = 2)))

overall_a_metrics = data.table::rbindlist(overall_a_metrics)
names(overall_a_metrics)=make.names(trimws(names(overall_a_metrics),"b"),unique = TRUE)

# Check your variable for NAs
xvars = overall_a_metrics %>% select_if(~ !any(is.na(.)))
xvars$StandID = plotList$Comment

# Write the las metrics outputs to a .csv
write.csv(xvars, paste0(dir, "output/lidR_Metrics.csv"))

# Read in our plot data (basal area, volume, etc.) processed from FVS, filter it for the 
# inventory year, and merge it with the predictor variables from above into a new dataframe
plot_attr= read_xls(paste0(dir, "input/data/", "FVSOut.xls"), sheet = "FVS_Compute")
plot_attr = plot_attr %>% dplyr::group_by(StandID) %>% slice(which.min(Year))

df=inner_join(plot_attr, xvars, by = c("StandID" = "StandID"))


# Construct a random forest model
# Caution: this step takes fairly long but can be shortened by setting importance=FALSE
# First we might need to remove columns with NAs
df = df[ , apply(df, 2, function(x) !any(is.na(x)))]
# And then we use bortuta to select the predictive features
boruta <- Boruta(x=df[ , -which(names(df) %in% names(plot_attr))], y=df$SBA/10.117136203, doTrace = 0, pValue = 0.05)
# Explore the importance of all relevant features 
plot(boruta)
# Subset the features to only those selected by Boruta
finalvars = getSelectedAttributes(boruta, withTentative = F)
# Fit the randomForest model
(modelRF_BA = randomForest(x=df[ , finalvars], y=df$SBA/10.117136203, # This get's the response 
                           # variable to be BA per 20m pixel
                           importance = TRUE))

# Model Inspection
# You have to pick a model to view, so I picked TPA in the 6" class and larger, modelRF_TPA.GT.5IN
modelRF = modelRF_BA
# 9.a Check error convergence
plot(modelRF)
# 9.b Plot mean decrease in variable importance and node purity
varImpPlot(modelRF)

# Preform a permutation test cross-validation for the random forest model
( modelRF_BA.cv = rf.crossValidation(modelRF_BA, df[ , finalvars], p=0.10, n=99, ntree=501) )
# This will also write the cross validation results to a text file for later reference
sink(paste0(dir, "output/","modelRF_BA.cv.txt"))
print("modelRF_BA permutation test cross-validation")
print(modelRF_BA.cv)
sink()

# Now we build a brick of the predictor variables for the entire catalog (for 
# prediction) using a custom myMetrics wrapper function and catalog_apply
myMetrics_Catalog = function(cluster, res)
{
  las = readLAS(cluster)
  if (is.empty(las)) return(NULL)
  las = normalize_height(las, knnidw(k = 8, p = 2))
  
  las = decimate_points(las, homogenize(1,24))
  
  area_metrics = pixel_metrics(las, myMetrics(X, Y, Z, Intensity, ReturnNumber, Classification, 
                                              dz = 1, th = 2, minht = 1.37, above = 2), res=20)
  
  #metrics = c(area_metrics)
  
  terra::crs(area_metrics) = "epsg:26912"
  bbox = lidR::extent(cluster)
  metrics = terra::crop(area_metrics, bbox) # This removes the points in the buffered area
  
  return(metrics)
}

output_myMetrics = catalog_apply(project, myMetrics_Catalog)
names = output_myMetrics[[1]]@ptr[["names"]]
output_myMetrics = do.call(terra::merge, output_myMetrics)
names(output_myMetrics) = names

# We only need the layers that the model needs, and we should check for NA's 
output_myMetrics = subset(output_myMetrics, names(modelRF$importance[,1]))

# Make predictions and plot it
predRF = terra::predict(object=output_myMetrics, model=modelRF, type="response", progress="window", na.rm=TRUE)

terra::crs(predRF) = "epsg:26912"
plot(predRF, col=rev(heat.colors(10)))

# Write the Basal Area raster to GeoTiff
writeRaster(predRF, filename=paste0(dir, "output/", "BA.GT.5IN.tif"), 
            gdal=c("COMPRESS=NONE", "TFW=YES","of=COG"), datatype='INT1U', overwrite=TRUE)

# Make a more fancy plot
library(leaflet) # leaflet still depends on raster - https://github.com/e-sensing/sits/issues/435
predRF = terra::project(predRF, "+proj=longlat +datum=WGS84")
# 13.a If you have large areas that are treeless or water or bare ground, you might want to crop the 
# prediction raster to the canopy layer we made earlier
# canopy = terra::project(canopy, "+proj=longlat +datum=WGS84")
# predBA = terra::crop(predBA,canopy)

# Create a palette
pal_ba = colorNumeric(colorRampPalette(c("#000004FF","#51127CFF","#B63679FF","#FB8861FF","#FCFDBFFF"))(100), terra::values(predRF),
                      na.color = "transparent")
m = leaflet(predRF) %>% addTiles %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Esri.WorldImagery") %>%
  addProviderTiles(providers$OpenTopoMap, group = "OpenTopoMap") %>%
  addRasterImage(raster::raster(predRF), colors = pal_ba, opacity = 0.75, group = "Basal Area", maxBytes = 30*1024*1024) %>%
  addLayersControl("topleft",
                   baseGroups = c("Esri.WorldImagery", "OpenTopoMap"),
                   overlayGroups = c("Basal Area"),
                   options = layersControlOptions(collapsed = TRUE)) %>% 
  addLegend(pal = pal_ba, values = terra::values(predRF),
            title = "Basal area per pixel (20m)")
m

