# Introduction to lidR workshop 
# Ecological Restoration Institute 
# Day 1 - Introduction to R and lidR
# 8/22/23


##########################################################################################
# For the first portion of this code we're going to practice manipulating dataframes to get 
# comfortable with R syntax and RStudio. We'll then get into reading and visualizing some 
# ALS plots. Play around with any parts of the code you're curious or confused about!
##########################################################################################

# load your libraries 

library(lidR)
library(sf)
library(ggplot2)
library(dplyr)

#set this to the path to where the files are on your machine 
dir = "C:/Users/mr2988/lidR_workshop/"

#We have a csv file with all the trees surveyed across each 0.1 acre plot. We're going to 
#take this and find the BA for each plot to explore some basic data manipulation

#Here are all the trees surveyed on the plots for this study area  
all_trees = read.csv(paste0(dir, "input/data/FUWI_2018_treelocs.csv"), fileEncoding = 'latin1')

#Let's look at what we're working with 
summary(all_trees)
#You'll see that there are a lot of columns of data that we don't really need. You can also 
#click on "all_trees" in your Environment tab in RStudio to look at all the data.

#Let's grab the basics so we can calculate plot BA and keep the locations so we can plot them 
#in the next section. I'm selecting by index, but you could select by column name too if you
#want to type out all the column names 
#all_trees['height.m']
all_trees = all_trees[c(2:5, 7:13)]

#Now let's change how these plots are named so it matches the naming convention from our 
#shapefile that we'll be looking at in the next section. 
all_trees$plotID = paste0("FUWI-", all_trees$block, "-", all_trees$treatment, "-", all_trees$plot)
head(all_trees)

#Now let's get into some of the ways we can transform data with dplyr. We'll start by 
#just rearranging our plots so it groups together trees on the same plot. Run this and
#then check out all_trees again 
all_trees = all_trees %>% arrange(desc(plotID)) 

#The %>% is called a pipe and it is an operator and it feeds whatever is left of the pipe
#into the function to the right of the pipe. The previous line is equivalent to: 
#arrange(all_trees, desc(plotID))

PlotBA = all_trees %>% filter((tree.status == (1 | 2)) & (!is.na(DBH.cm))) %>% #select live trees (including declining) with diameter measurements (saplings have NA dbh)
                       mutate(BA = 0.00107639*pi*(DBH.cm/2)^2) %>% #create a new column called BA and calculate from dbh and convert to sq ft 
                       group_by(plotID, block, treatment) %>% #the following calculation will summarize by plot, keeping the block & treatment 
                       summarise(BA = sum(BA)) #and add up all the tree's BAs for a plot BA 

#Now let's look at our result
head(PlotBA)
hist(PlotBA$BA, breaks = 10)
#We can also look at how this breaks down by treatment type 
PlotBA$treatment = as.factor(PlotBA$treatment)
ggplot(PlotBA, aes(x = treatment, y = BA)) +
  geom_boxplot() +
  labs(title = "Basal area by treatment")

#We can also look at the data aggregated across plots in the same block and treatment type 
treatmentBA = PlotBA %>% group_by(block, treatment) %>% summarise(BA = mean(BA))
treatmentBA
#Now make a boxplot for this summarized data  


#For some extra practice see if you can get TPA for each 0.1 acre plot  
#Don't forget that there are live trees and dead trees and some of those live trees are saplings (look at heights)
#The n() function gets you the tally:
all_trees %>% group_by(plotID, tree.status) %>% summarise(tree_count = n())

###########################################################################################
#Now let's get into visualizing the point clouds for those plots   
#Don't worry too much about this chunk of code, we'll cover this in more detail tomorrow. 
project = readALSLAScatalog(paste0(dir, "input/las/"))
st_crs(project) = 26912
opt_chunk_buffer(project) = 15        
opt_chunk_size(project) = 250         
opt_select(project) = "xyzicra"        
opt_stop_early(project) = FALSE       
opt_wall_to_wall(project) = TRUE
opt_progress(project) = TRUE

#read in the shapefile with your GPS-ed plot locations 
plots = read_sf(paste0(dir, "input/vector/", "FUWI_Plots_WGS84_4326_Subset.shp"))
#set the projection to match the projection set for your ALS data 
plots =st_transform(plots, sf::st_crs(project)) 
#clip out each plot. The plot radius is ~12m, so this has a 7m buffer 
plots = st_buffer(plots, dist=20) 

#select the tiles from our ALS catalog that have plots on them
tile = catalog_intersect(project, plots)
#and then clip out those plots 
las_list = clip_roi(tile, plots)

#Let's look at our first plot 
n = 1
las = las_list[[n]]
PlotID = plots$Comment[[n]]

summary(las)
plot(las)

#You can rotate, zoom in and out, shift the plot around with your mouse to get a full view 
#Some of the major color palettes are: viridis, topo.colors, terrain.colors, rainbow, heat.colors
#Some show up better against a white background instead of the default black 
plot(las, pal = magma, bg = "white")  

#You can also adjust the point size, which can help make sense of the scene in sparse point clouds 
plot(las, point = 3)

#The default is to color by the Z values, but you can color by any of the las attributes
#Looking at the plot by classification, can you guess what the classifications are?
plot(las, color = "Classification")


#You can also add an axis and legend to your plot
#Try this with and without the "breaks" argument and see what difference it makes 
plot(las, color = "Intensity", breaks = "quantile", axis = TRUE, legend = TRUE)

#Instead of plotting points, you can also plot as voxels, which are 3D pixels 
#The "res" argument is the size of the voxel in the same units as the las
vox = voxelize_points(las, res = 0.5) 
plot(vox, voxel = TRUE, bg = "white")

#We can also plot a transect using the clip_transect function and ggplot.
#Test out some different starting and ending points and widths then 
#See if you can figure out how to turn this into a function that takes your 
#las, starting and ending points, and transect width as input and returns 
#your plot
p1 = c(min(las@data$X), min(las@data$Y))
p2 = c(max(las@data$X), max(las@data$Y))       
las_tr = clip_transect(las, p1, p2, width  = 4, xz = TRUE)

p = ggplot(las_tr@data, aes(X,Z, color = Z)) + 
           geom_point(size = 1) + 
           coord_equal() + 
           theme_minimal() +
           scale_color_gradientn(colours = height.colors(50))
p

#Here is a sample function that adds three numbers together so you can see the syntax
function_name = function(input1, input2, input3){
  output = input1 + input2 + input3
  return(output)
}
function_name(1,2,3)

#rast_las = rasterize_canopy(las, res = 0.5)
#plot(rast_las)

