
################################################################################
#
# Journal: Molecular Ecology
#
# Article title: "Corridor-based approach with spatial cross validation reveals 
#                scale-dependent effects of geographic distance, 
#                human footprint, and canopy cover on grizzly bear genetic 
#                connectivity"
#
# Authors: Eric Palm, Erin Landguth, Zachary Holden, Casey Day, Clayton Lamb, 
#          Paul Frame, Andrea Morehouse, Garth Mowat, Michael Proctor, 
#          Michael Sawaya, Gordon Stenhouse, Jesse Whittington, Katherine Zeller
#
################################################################################

# Script for extracting mean covariates along least cost paths between
# pairwise locations of grizzly bear genetic samples. This code uses the real
# spatial prediction from Script #4: "model_straight_40_km.rds" as the input
# resistance surface. The input genetic data
# is an example dataset with randomly generated spatial coordinates.

# Code written by Eric Palm and Erin Landguth

################################################################################

# Load packages

require(raster)
require(foreach)
require(doParallel)
require(sp)
require(sf)
require(gdistance)
require(exactextractr)
require(dplyr)
require(terra)


# Load environmental covariate data
env <- list.files("env/", full.names = T, pattern = ".tif$") %>% 
  raster::stack() 

# Load in example genetic data 
gen_table <- readRDS("data/processed/example/pairwise_data_example.rds") %>% 
  dplyr::mutate(index = row_number())


# Set the number of cores to use for parallel processing
noCLS <- 20

# Read in the resistance raster produced from the 40 km straight-line model 
pred_res <- raster::raster("results/real/pred_straight_40_km.tif")
  
# Create the transition layer to use with least cost path 
# take 1/mean(x) to convert resistance values to conductance as required by the 
# gdistance::shortestPath function
tr_predr <- gdistance::transition(pred_res, function(x) 1/mean(x), directions=8) 
  
# Geocorrect the transition layer
tr_predr <- gdistance::geoCorrection(tr_predr, type = "c", multpl=F)
  
# Set up parallel here
cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)
  
# Use this index in the foreach loop
iseq <- gen_table$index
  
# Begin loop through each start and end point, calculate LCP, extract
# mean covariate value
system.time(  
  env_table <- foreach(i=iseq, .combine=rbind, .verbose = T,
                       .packages = c("gdistance", "exactextractr", "sp", "sf", 
                                     "raster", "dplyr", "terra")) %dopar% 
    {
      index <- gen_table$index[i]
      
      y_1 <- gen_table$y_1[i]
      x_1 <- gen_table$x_1[i]
      y_2 <- gen_table$y_2[i]
      x_2 <- gen_table$x_2[i]
      start_points <- cbind(x_1, y_1) %>% sp::SpatialPoints()
      end_points <- cbind(x_2, y_2) %>% sp::SpatialPoints()
      
      # Some sets of points are within the same raster cell, which will cause 
      # the shortestPath function to fail. For those, we just buffer the start 
      # and end points by 1 km, then union the polygons and extract covariates 
      # within that merged polygon
      if (gen_table$euc_geog[i] < 510)
      {

        lcp_buff <- sf::st_buffer(sf::st_as_sfc(start_points), 1000) %>%
          st_union(., sf::st_buffer(sf::st_as_sfc(end_points), 1000)) %>%
          sf::st_set_crs(26911)
        
      } else {

        lcp <- gdistance::shortestPath(tr_predr, start_points, end_points, 
                                       output = "SpatialLines")
        
        lcp_buff <- sf::st_as_sf(lcp) %>% 
          sf::st_buffer(., 1000) %>% 
          sf::st_set_crs(26911)
      }
      
      # exactextractr prefers terra rasts over rasters
      # You can create terra objects within foreach loops, but you can't import 
      # terra objects created outside the foreach loop because they use C++
      
        line_vals <- exactextractr::exact_extract(terra::rast(env), 
                                                  lcp_buff)[[1]][, -(nlayers(env) + 1)]
        
        #  gather all metrics from line or buffer, updating names along way
        mean_sl <- apply(line_vals, 2, FUN=mean)
      
      # Appending to row vector type
      out <- append(index, mean_sl)
      return(out)
    }
)  
 
# Make sure to close the connection for parallel
parallel::stopCluster(cl)

# Clean up data frame of results
env_table_names <- as.data.frame(env_table) %>%
  `rownames<-`(NULL) %>% 
  `names<-`(c("index", names(env))) %>% 
  dplyr::as_tibble()
  
# Add IDs, genetic and geographic distances to the new data frame.
env_table_w_extract <- dplyr::inner_join(env_table_names, gen_table) 
  
# Save extracted values to file for use in LCP models  
saveRDS(env_table_w_extract,
        "results/example/extracted_lcp_example.rds")
  