#Modeling 

library(lidR)
library(mapview)
library(stats)
library(readxl)
library(dplyr)
library(corrplot)
library(randomForest) 
library(rfUtilities)  # remotes::install_github("jeffreyevans/rfUtilities", force = TRUE)
library(Boruta)       
library(lidRattR)     # install_github("RCBlackburn/lidRattR")

dir = "C:/Users/mr2988//lidR_4.1_tutorials/"

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

# Read in the shapefile of the inventory plot locations and buffer to the plot radius 
pointList = read_sf(paste0(dir, "input/vector/", "FUWI_Plots_WGS84_4326_Subset.shp"))
pointList = st_transform(pointList, sf::st_crs(project))
plotList = st_buffer(pointList, dist=11.34970)

# Clip the catalog to the plot boundaries and and normalize the point clouds 
new_ctg = lidR::clip_roi(project, plotList)
new_ctg = lapply(new_ctg, function(x) normalize_height(x, knnidw(k = 8, p = 2)))
new_ctg = lapply(new_ctg, function(x) filter_poi(x, Classification != 2))
new_ctg = lapply(new_ctg, function(x) decimate_points(x, homogenize(1,24)))

# Calculate metrics for each plot 
plots_metrics <- lapply(new_ctg, cloud_metrics, func = .stdmetrics)
plots_metrics <- data.table::rbindlist(plots_metrics)
plots = cbind(plotList, plots_metrics)

# Read in our FVS output so we can model biomass across the study area 
FVS_plots = read_xls(paste0(dir, "input/data/FVSOut.xls"), sheet = "FVS_Compute")
biomass = FVS_plots[c("StandID", "STBIOMSS")]
plots = inner_join(biomass, plots, join_by(StandID == Comment))

# Clean up our data so we can plot a correlation matrix to select our predictor variables 
plots_numeric = plots %>% select(where(is.numeric)) %>% select_if(~ !any(is.na(.)))
plots_numeric = plots_numeric[, colSums(plots_numeric != 0) > 0]

# And now we can look at our correlation matrix. 
cor_mat = cor(plots_numeric)
corrplot(cor_mat)
# We'll also list out the 5 strongest correlations (positive and negative)
as.data.frame(cor_mat) %>% arrange(desc(STBIOMSS)) %>% head()
as.data.frame(cor_mat) %>% arrange(STBIOMSS) %>% head()

#Build a few linear models from a few highly correlated variables and pick the best one 
m1 = lm(STBIOMSS ~ zq75, data = plots)
summary(m1)
m2 = lm(STBIOMSS ~ zq75 + p3th, data = plots)
summary(m2)
m3 = lm(STBIOMSS ~ zq75 + ipcumzq70, data = plots)
summary(m3)

#Plot our observed vs predicted biomass for m3
plot(plots$STBIOMSS, predict(m3))
abline(0,1)

#Now we can apply this plot-based model for wall-to-wall coverage across our study area
#first we normalize the catalog 
opt_output_files(project) <- paste0(dir, "output/normalized_ctg")
output <- normalize_height(project, tin())
#then we compute our metrics
opt_output_files(project) <- ""
metrics_w2w <- grid_metrics(project, .stdmetrics, res = 20)

#Plot our predictor variables 
plot(metrics_w2w$zq75)
plot(metrics_w2w$p3th)
plot(metrics_w2w$ipcumzq70)

#Make our predictions and plot!
STBIOMSS_pred = coef(m3)[1] + metrics_w2w$zq75*coef(m3)[2] + metrics_w2w$ipcumzq70*coef(m3)[3]
plot(STBIOMSS_pred)

################################################################################################
#Now let's kick it up a notch and try modeling without having to manually select our variables
#and as a random forest model instead of a linear one. We'll use the package boruta to do our 
#variable selection. We already have our project and plots loaded so we won't do that again

# Define a custom set of las metrics to extract for each plot (this is basically .stdmetrics, 
# repeated to give you an example of a custom metrics function)
source(paste0(dir, "scripts/", "myMetricsScript.R"))

myMetrics = function(x, y, z, i, rn, class, dz, th, minht, above)
{
  S  = stdmetrics(x,y,z,i,rn,class,dz,th)
  ShpM = stdshapemetrics(x,y,z)
  # C  = stdmetrics_ctrl(x, y, z, a)
  # Z  = stdmetrics_z(z, dz)
  # I  = stdmetrics_i(i, z, class, rn)
  # RN = stdmetrics_rn(rn, class)
  mSM = mystdMetrics(z,i,rn,minht,above)
  
  metrics = c(S,mSM,ShpM)
  return(metrics)
}

# calculate our custom myMetrics for each plot and store the results in a dataframe with the associated plot ID
overall_a_metrics = lapply(new_ctg, function(x) cloud_metrics(x, myMetrics(X, Y, Z, Intensity, ReturnNumber, Classification, 
                                                                           dz = 1, th = 2, minht = 1.37, above = 2)))

overall_a_metrics = data.table::rbindlist(overall_a_metrics)
names(overall_a_metrics)=make.names(trimws(names(overall_a_metrics),"b"),unique = TRUE)

new_ctg = lapply(new_ctg, function(x) segment_trees(x, li2012()))
new_ctg = lapply(new_ctg, function(x) ind_tree(x))
overall_t_metrics = lapply(new_ctg, function(x) tree_sum(x))
overall_t_metrics = data.table::rbindlist(overall_t_metrics)

new_ctg = lapply(new_ctg, function(x) ind_vox(x, res = 2))
overall_v_metrics <- lapply(new_ctg, function(x) vox_sum(x, 2))
overall_v_metrics = data.table::rbindlist(overall_v_metrics)

xvars = do.call(cbind,list(overall_a_metrics,overall_t_metrics,overall_v_metrics))
xvars = cbind("StandID"=plotList$Comment, xvars)

# Check your variable for NAs
xvars = xvars %>% select_if(~ !any(is.na(.)))

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
  las = segment_trees(las, li2012())
  las_tree = ind_tree(las)
  tree_metrics <- pixel_metrics(las_tree, tree_sum_raster(Z, npoints, ca), res = 20)
  
  las_vox = ind_vox(las, res = 2)
  voxel_metrics = pixel_metrics(las_vox, 
                                func = vox_sum_raster(SVi_2, FRDi_2, PDi_2, PDi_above_2, vox_res = 2), 
                                res = 20)
  
  metrics = c(area_metrics, tree_metrics, voxel_metrics)
  
  terra::crs(metrics) = "epsg:26912"
  bbox = lidR::extent(cluster)
  metrics = terra::crop(metrics, bbox) # This removes the points in the buffered area
  
  return(metrics)
}
output_myMetrics = catalog_apply(project, myMetrics_Catalog)
output_myMetrics = do.call(terra::merge, output_myMetrics)

# We only need the layers that the model needs, and we should check for NA's (and likely drop zentropy)
output_myMetrics = subset(output_myMetrics, names(modelRF$importance[,1]))
# output_myMetrics = terra::subset(output_myMetrics, names(output_myMetrics)[names(output_myMetrics) != "zentropy"])
# xvars = subset(xvars, select = -c(zentropy))

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

