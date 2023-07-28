
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

# Script for running gradient boosting machine models using data from straight-
# line covariate extractions. These models use the caret and CAST packages to
# incorporate spatial cross validation and variable (feature) selection.

# Code written by Eric Palm and Zack Holden

################################################################################

# Load packages
require(foreach)
require(doParallel)
require(caret)
require(raster)
require(dplyr)
require(CAST)
require(gbm)
require(sp)
require(sf)

### Load environmental covariate data
env <- list.files("env/", full.names = T, pattern = ".tif$") %>% 
  raster::stack() 

# Load in genetic distances and straight-line covariate extractions from the 
# full dataset and only use pairs of locations less than 40 km apart
dat <- readRDS("results/real/extracted_straight_full_dataset.rds") %>% 
  dplyr::filter(euc_geog <= 40000) %>% 
  dplyr::mutate(index = row_number())

# Pull out the covariate data as predictors
traindat <- dat %>% 
  dplyr::select(canopy_cover:tri, euc_geog) 

# Pull out the genetic distances as the response data
respdat <- dat %>% 
  dplyr::pull(GD)

# Create a raster with a constant value corresponding to the median
# geographic distance and add it to the raster stack
env_all <- raster::setValues(env$evi, median(dat$euc_geog)) %>%
  `names<-`("euc_geog") %>%
  raster::stack(env, .) %>%
  raster::mask(., env$evi)

# Read in the list containing bear IDs for each of the 10 spatial clusters
cluster_list <- readRDS("data/processed/real/spatial_cluster_list.rds")

# Empty list for storing results of for loop
indices_cv <- list()

# For each of the 10 vectors (training folds) in the spatial clusters list,
# subset only those rows in the pairwise dataframe where either bear 1 or 2
# are not included in that spatial cluster.
# And place the indices (row numbers) of those observations in a list.
# The output here is a list with 10 elements, each with the indices from the 
# pairwise dataframe corresponding to each training dataset (spatial clusters).
for (i in 1:length(cluster_list)){
  indices_cv[[i]] <- dat %>% 
    dplyr::filter(!(id_1 %in% cluster_list[[i]]$animal_id | 
                      id_2 %in% cluster_list[[i]]$animal_id)) %>% 
    dplyr::pull(index)
}

# Specify model cross-validation parameters here
# Add in the indices_cv list for the indices of training folds
ctrl_fit <- caret::trainControl(index = indices_cv,
                                allowParallel = TRUE,
                                method = "cv",
                                classProbs = FALSE,
                                savePredictions = T,
                                returnResamp = "all")

# Choose how many clusters you want to run for parallelization
cl <- parallel::makeCluster(60)
doParallel::registerDoParallel(cl)

# Set a seed for reproducibility
set.seed(1234)

# Define a search grid for n.trees, shrinkage and n.minobsinnode
# This model took an hour with 60 cores on a fast machine
system.time(
  gbm_ffs <- CAST::ffs(traindat, respdat,
                       method = "gbm",
                       trControl = ctrl_fit,
                       metric = "RMSE", 
                       minVar = 2,
                       tuneGrid=expand.grid(interaction.depth = c(2),
                                            n.trees = seq(100, 400, 20),
                                            shrinkage = c(.01),
                                            n.minobsinnode = c(10)))
)

# To run a full model without variable selection, use caret::train, as below
# system.time(
#   gbm_full <- caret::train(traindat, respdat,
#                        method = "gbm",
#                        trControl = ctrl_fit,
#                        metric = "RMSE", 
#                        tuneGrid=expand.grid(interaction.depth = c(2),
#                                             n.trees = seq(100, 400, 20),
#                                             shrinkage = c(.01),
#                                             n.minobsinnode = c(10)))
# )

# Make sure to stop the cluster
parallel::stopCluster(cl)

# Save the model output
saveRDS(gbm_ffs, "results/real/model_straight_40_km.rds")

# Save the selected variables
vars <- gbm_ffs$selectedvars

# View the performance of each model fit in the feature selection process.
gbm_ffs$perf_all 

# View the average RMSE across each of the 10 test folds for each combination
# of hyperparameters in the grid search
gbm_ffs

# Make spatial predictions from your model using the raster stack as new data
pred_rast <- predict(env_all[[vars]], gbm_ffs)

# Save predictions to a resistance raster for use in the LCP extraction step
raster::writeRaster(pred_rast, "results/real/pred_straight_40_km.tif", overwrite=T, 
                    datatype = "FLT4S")

