
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

# Script for extracting mean covariates along straight lines connecting
# pairwise locations of grizzly bear genetic samples. The input genetic data
# is an example dataset with randomly generated spatial coordinates.

# Code written by Eric Palm and Erin Landguth

################################################################################

# Load packages
require(raster)
require(foreach)
require(doParallel)
require(sf)
require(exactextractr)
require(dplyr)
require(terra)


# Load environmental covariate data
env <- list.files("env/", full.names = T, pattern = ".tif$") %>% 
  raster::stack() 

# Load in example genetic data
gen_table <- readRDS("data/processed/example/pairwise_data_example.rds") 

# Set the number of cores to use for parallel processing
noCLS <- 20

# Set up parallel here
cl <- parallel::makeCluster(noCLS) 
doParallel::registerDoParallel(cl)

# Use this index in the foreach loop
iseq <- gen_table$index

# Begin loop through each start and end point, calculate straight line, extract
# mean covariate value
system.time(
  env_table <- foreach::foreach(i=iseq, .combine=rbind, .verbose = T, 
                       .packages = c("sp", "sf", "exactextractr", "raster", 
                                     "terra")) %dopar% 
    {
      index <- gen_table$index[i]
     
      y_1 <- gen_table$y_1[i]
      x_1 <- gen_table$x_1[i]
      y_2 <- gen_table$y_2[i]
      x_2 <- gen_table$x_2[i]
      ys <- rbind(y_1, y_2)
      xs <- rbind(x_1, x_2)
      
      points <- sp::SpatialPoints(cbind(xs, ys))
      
      # Create straight lines between points and buffer by 1 km
      sl <- as(points, "SpatialLines") %>%
        sf::st_as_sf(.) %>% 
        sf::st_set_crs(26911) %>% 
        sf::st_buffer(., 1000)
      
      # exactextractr prefers terra rasts over rasters
      # You can create terra objects within foreach loops, but you can't import 
      # terra objects created outside
      # the foreach loop because they use C++
      line_vals <- exactextractr::exact_extract(terra::rast(env), 
                                                sl)[[1]][, - (nlayers(env)+1)]
      
      # Calculate means of pixels within buffer polygon
      mean_sl <- apply(line_vals, 2, FUN=mean)
      
      # Append to row vector type
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
  
# Save extracted values to file for use in straight line models
saveRDS(env_table_w_extract,
        "results/example/extracted_straight_example.rds")
