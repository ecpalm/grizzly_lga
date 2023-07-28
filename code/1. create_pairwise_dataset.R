
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

# Script for taking an example dataset with actual grizzly bear genotypes
# and randomly generated spatial coordinates and creating a pairwise
# data frame for use in gradient boosting machine models with genetic distance
# as the response variable and mean environmental covariates extracted along
# transects (straight-line or LCPs) as the predictor variables.

# Code written by Eric Palm

################################################################################

# Load packages
require(gstudio)
require(sf)
require(dplyr)
require(purrr)

#####  FUNCTIONS ###############################################################

# Function to take a geographic distance matrix and convert it to a data frame
# where each row has a pair of animals and the geographic distance between them
lower_df <-
  function (x) {
    mat <- as.matrix(x)
    ind <- which(lower.tri(mat, diag = F), arr.ind = TRUE)
    nn <- dimnames(mat)
    tibble(id_1 = nn[[1]][ind[, 1]],
           id_2 = nn[[2]][ind[, 2]],
           d = mat[ind])
  }

# Function to normalize the genetic distances to be between 0 and 1.
norm <- function (x) {
  (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}

################################################################################
# Load in genetic data.
# This is a subset (n = 100) of the 1161 grizzly bear genotypes.
# Real spatial coordinates were replaced by random coordinates for this example.
gb_gen <- readRDS("data/raw/example/genetic_data_with_coords_example.rds")

d_geo <-
  gb_gen %>% 
  dplyr::select(x,y) %>% 
  dist(., method = 'euclidean') %>% 
  lower_df() %>% 
  dplyr::rename(euc_geog = d) %>% 
  dplyr::mutate(across(id_1:id_2, ~as.numeric(.)),
                id_1 = gb_gen$animal_id[id_1],
                id_2 = gb_gen$animal_id[id_2]) 

# Set a seed for reproducibility
set.seed(1234)

# Create a list of animal IDs in each of 10 spatial folds that will be used in the
# spatial cross validation during GBM models
cluster_list <-
  gb_gen %>% 
  dplyr::filter(animal_id %in% unique(c(d_geo$id_1, d_geo$id_2))) %>% 
  dplyr::mutate(cluster = kmeans(cbind(.$x, .$y), 10)$cluster) %>% 
  dplyr::select(animal_id, cluster) %>% 
  dplyr::nest_by(cluster) %>% 
  dplyr::pull(data)

# Save the list
saveRDS(cluster_list, "data/processed/example/spatial_cluster_example.rds")

# Format genotype data for the 'gstudio' package to calculate genetic distances
gb_gen_locus <-
  gb_gen %>% 
  dplyr::select(animal_id, CXX110:MU50) %>% 
  dplyr::mutate(across(-1, ~paste0(str_sub(., 1, 3), ":", str_sub(., 5, 7))),
                across(-1, ~locus(., type = "separated"))) %>% 
  as.data.frame()

# Calculate euclidean genetic distance and 1 - proportion of shared alleles
euc <- gstudio::genetic_distance(gb_gen_locus, mode = "Euclidean", 
                                 stratum = "animal_id") %>% 
  lower_df() %>% 
  dplyr::rename(euc_gen = d)

Dps_1 <- gstudio::genetic_distance(gb_gen_locus, mode = "Dps", 
                                   stratum = "animal_id") %>% 
  lower_df() %>% 
  dplyr::rename(Dps_1 = d) %>% 
  dplyr::mutate(Dps_1 = 1 - Dps_1)

# Join these genetic distances to the pairwise data frame with bear IDs
gen_d_all <-
  purrr::reduce(list(euc, Dps_1), dplyr::inner_join)


# Sorting the pairwise data frame so that for each pair of IDs, 
# the two IDs and their associated spatial coordinates are sorted the same way
# in this data frame as in the geographic distances data frame so they can be joined

gb_gd_temp <-
  gen_d_all %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(p1 = pmin(id_1, id_2),
                p2 = pmax(id_1, id_2)) %>% 
  dplyr::select(-c(id_1, id_2)) %>% 
  dplyr::rename(id_1 = p1, id_2 = p2) %>% 
  dplyr::mutate(x_1 = gb_gen$x[gb_gen$animal_id == id_1],
                y_1 = gb_gen$y[gb_gen$animal_id == id_1],
                x_2 = gb_gen$x[gb_gen$animal_id == id_2],
                y_2 = gb_gen$y[gb_gen$animal_id == id_2]) %>% 
  dplyr::ungroup() 

# Taking the first principal component of euclidean genetic distance and 1 - proportion of shared alleles
# and normalizing the result to be between 0 and 1
# Confirm that GD is positively correlated with euc_gen and Dps. If not, change "-1" to "1".

gb_gd_pc <-
  gb_gd_temp %>% 
  dplyr::mutate(GD = (-1*(prcomp(gb_gd_temp %>% 
                                   dplyr::select(euc_gen, Dps_1), scale = T, center = T)$x[,1])) %>%
                  norm(.))

gb_gd_pc %>% 
  dplyr::select(euc_gen, Dps_1, GD) %>% 
  cor()

# Joining the data frames with genetic distances and geographic distances
# Make sure the animal ids are sorted so id_1 is earlier in alphabetical sorting
# than id_2
gb_gd <- 
  gb_gd_pc %>% 
  dplyr::inner_join(., d_geo %>%
                      rowwise() %>%
                      dplyr::mutate(p1 = pmin(id_1, id_2),
                                    p2 = pmax(id_1, id_2)) %>%
                      dplyr::ungroup() %>%
                      dplyr::select(-c(id_1, id_2)) %>%
                      dplyr::rename(id_1 = p1, id_2 = p2))


# Rearrange columns
# Add in index of row number
gb_gd %>% 
  dplyr::select(id_1:y_2, euc_gen, Dps_1, GD, euc_geog) %>% 
  dplyr::mutate(index = row_number()) %>% 
  saveRDS(., "data/processed/example/pairwise_data_example.rds")

# Use this final dataset in the straight line covariate extraction script

