#Tree Segmentation 


library(lidR)
library(mapview)
library(terra)
library(raster)

plan(multisession, workers = (parallel::detectCores(logical=FALSE)-2))
set_lidr_threads(1L)

dir = "C:/Users/mr2988//lidR_4.1_tutorials/"

las<-readLAS(paste0(projDir, "/input/las/FUWI_retile_437500_3903000.las"), select = "xyzicnr")
las@data$Classification<-0
las<-classify_ground(las, csf(sloop_smooth = FALSE, class_threshold = 0.1, cloth_resolution = 0.33,
                              rigidness = 1L, iterations = 500L, time_step = 0.65))


dem<-rasterize_terrain(las, res=0.5, algorithm = knnidw())
las<-normalize_height(las, knnidw())
las<-classify_noise(las, ivf(res=1, n=3))
las<-filter_poi(las, Classification != LASNOISE)
las<-filter_poi(las, Z>=0, Z<50)
plot(las)
#Where do you expect tree segmentation to work better? 


algo = pitfree(thresholds = c(0,10,20,30,40,50), subcircle = 0.2)
chm  = grid_canopy(las, 0.5, algo)
plot(chm, col = height.colors(50))

#####################################################################################
#Let's start by looking at individual tree detection
#Here we're using the Local Maximum Filter on the point cloud itself 
ttops <- locate_trees(las, lmf(ws = 5))
nrow(ttops) #How many trees were identified? 
#Adjust the window size (ws) and check how that changes the number of trees
#Why does the number of trees change with window size?

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

x <- plot(las, bg = "white", size = 4)
add_treetops3d(x, ttops)

#Now let's try again with a variable window size so that we use a larger ws for large trees 
#and a smaller ws for shorter trees
f <- function(x) {x * 0.1 + 3}
heights <- seq(0,30,5)
ws <- f(heights)
plot(heights, ws, type = "l", ylim = c(0,6))

ttops <- locate_trees(las, lmf(f))
nrow(ttops)

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

p <- plot(las, bg = "white", size = 4)
add_treetops3d(p, ttops)

#We can also add lines to show where the stems should be 

# First we extract the tree's coordinates and then 
# apply the shift to display the lines
# in the rendering coordinate system
x <- sf::st_coordinates(ttops)[,1] - p[1] 
y <- sf::st_coordinates(ttops)[,2] - p[2] 
z <- ttops$Z

# Build a GL_LINES matrix for fast rendering
x <- rep(x, each = 2)
y <- rep(y, each = 2)
tmp <- numeric(2*length(z)) 
tmp[2*1:length(z)] <- z
z <- tmp
M <- cbind(x,y,z)

# Display lines
rgl::segments3d(M, col = "black", lwd = 2)

#We can also find tree tops from our CHM raster instead of from the point cloud directly. 
#Let's explore a few different methods of rendering the CHM 

# Point-to-raster 2 resolutions
chm_p2r_05 <- rasterize_canopy(las, 0.5, p2r(subcircle = 0.2), pkg = "terra")
chm_p2r_1 <- rasterize_canopy(las, 1, p2r(subcircle = 0.2), pkg = "terra")

# Pitfree with and without subcircle tweak
chm_pitfree_05_1 <- rasterize_canopy(las, 0.5, pitfree(), pkg = "terra")
chm_pitfree_05_2 <- rasterize_canopy(las, 0.5, pitfree(subcircle = 0.2), pkg = "terra")

# Post-processing median filter
kernel <- matrix(1,3,3)
chm_p2r_05_smoothed <- terra::focal(chm_p2r_05, w = kernel, fun = median, na.rm = TRUE)
chm_p2r_1_smoothed <- terra::focal(chm_p2r_1, w = kernel, fun = median, na.rm = TRUE)

#Now let's detect trees from each of our CHMs. We'll reduce the redundancy of typing out the locate_trees() 
#function 6 times by using lapply() on our list of CHMs 
chm_list = list(chm_p2r_05, chm_p2r_1, chm_pitfree_05_1, chm_pitfree_05_2, chm_p2r_05_smoothed, chm_p2r_1_smoothed)
ttops_chm = lapply(chm_list, function(X) locate_trees(X, lmf(5)))

#Let's now plot 
par(mfrow=c(1,2))
col <- height.colors(50)
plot(chm_list[[1]], main = "CHM P2R 0.5", col = col); plot(sf::st_geometry(ttops_chm[[1]]), add = T, pch =3)
plot(chm_list[[2]], main = "CHM P2R 1", col = col); plot(sf::st_geometry(ttops_chm[[2]]), add = T, pch = 3)

plot(chm_list[[3]], main = "CHM P2R 0.5 smoothed", col = col); plot(sf::st_geometry(ttops_chm[[3]]), add = T, pch =3)
plot(chm_list[[4]], main = "CHM P2R 1 smoothed", col = col); plot(sf::st_geometry(ttops_chm[[4]]), add = T, pch =3)

plot(chm_list[[5]], main = "CHM PITFREE 1", col = col); plot(sf::st_geometry(ttops_chm[[5]]), add = T, pch =3)
plot(chm_list[[6]], main = "CHM PITFREE 2", col = col); plot(sf::st_geometry(ttops_chm[[6]]), add = T, pch =3)

#With an area this large, it is hard to see the differences in the number of trees, so let's pull that 
#out again to compare each method. Why does each method result in a different number of trees? 
ntrees = lapply(ttops_chm, nrow)
chm_trees = data.frame(CHM =  c("p2r_05", "p2r_1", "pitfree_05_1", "pitfree_05_2", "p2r_05_smoothed", "p2r_1_smoothed"), ntrees = as.array(ntrees))
chm_trees

#####################################################################################################
#Individual tree segmentation

#Just as we can detect trees from our CHM raster, dalponte2016 is a raster-based tree segmentation method
algo_dalponte <- dalponte2016(chm_list[[3]], ttops_chm[[3]])
algo_watershed = lidR::watershed(chm, th_tree = 4)

las <- segment_trees(las, algo_dalponte, attribute = "IDdalponte") # segment point cloud
las  = segment_trees(las, algo_watershed, attribute = "IDwatershed")

plot(las, bg = "white", color = "IDdalponte", colorPalette = pastel.colors(100)) # visualize trees
plot(las, bg = "white", color = "IDwatershed", pal = pastel.colors)

#Now let's look at the methods that work directly on the point cloud
algo_silva = silva2016(chm_list[[3]], ttops_chm[[3]])
algo_li = li2012()

#Watershed
las  = segment_trees(las, algo_silva, attribute = "IDsilva")
las  = segment_trees(las, algo_li, attribute = "IDli")

plot(las, bg = "white", color = "IDsilva", pal = pastel.colors)
plot(las, bg = "white", color = "IDli", pal = pastel.colors)

crowns_dalponte <- crown_metrics(las, func = NULL, attribute = "IDdalponte", geom = "concave")
crowns_watershed <- crown_metrics(las, func = NULL, attribute = "IDwatershed", geom = "concave")
crowns_silva <- crown_metrics(las, func = NULL, attribute = "IDsilva", geom = "concave")
crowns_li <- crown_metrics(las, func = NULL, attribute = "IDli", geom = "concave")

#Let's checkout the crown delineations for each method
par(mfrow=c(2,2),mar=rep(0,4))
plot(sf::st_geometry(crowns_dalponte), reset = TRUE)
plot(sf::st_geometry(crowns_watershed), reset = FALSE)
plot(sf::st_geometry(crowns_silva), reset = FALSE)
plot(sf::st_geometry(crowns_li), reset = FALSE)
#We can immediately see that li2012 did not perform as well as the other methods.
#Adjust the default parameters and see if you can get this method to match the others 


#Let's look at some other ways to visualize our data.

#When our trees are segmented, we can calculate canopy metrics for each segmented tree 
#We'll look at modeling with metrics later on, but for now we're using this to show another
#way to visualize the delineated canopies. 
crowns <- crown_metrics(las, func = .stdtreemetrics, attribute = "IDli", geom = "convex")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")
plot(crowns["Z"], main = "Height (convex hull)")

#We can also add stems to our tree tops 

