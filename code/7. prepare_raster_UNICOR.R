
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

# Script for taking the final resistance raster from full dataset straight-line
# model and converting it to a format suitable to create resistant kernel
# connectivity models in UNICOR.

# Code written by Eric Palm, Kathy Zeller and Adam Ford

################################################################################

# Load packages
require(terra)
require(dplyr)

# Load in predicted resistance surface from the full dataset straight-line model
# Normalize the values to be between 0 and 1
r <- terra::rast("results/real/pred_straight_full.tif") %>% 
  app(., function(x) 
    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)))

# Write this normalized raster to an .asc file for UNICOR
terra::writeRaster(r, "unicor/asc/resistance_full_dataset.asc", 
                   datatype = "FLT4S", overwrite = T) 

# Pull out the different components of the the .asc text/header file
tx <- readLines("unicor/asc/resistance_full_dataset.asc")
header <- c(tolower(tx[(1:5)]), tx[6]) 
bdy <- (tx)[-(1:6)]

# Then rewrite the .asc file in the preferred format for UNICOR
writeLines(append(header, bdy, after = 6), 
           "unicor/asc/resistance_full_dataset.asc")

# Consult UNICOR user manual for details on how to use the software
# and for input definitions:
# https://www.fs.usda.gov/rm/pubs_other/rmrs_2011_landguth_e002.pdf

# For the UNICOR example in this repository, be sure to install UNICOR in the
# "unicor" folder, which also contains the resistance raster (created in this
# script), a .csv file with starting locations, and a .rip used for setting
# inputs and running the resistant kernels in UNICOR.

# The "resistance_full_dataset.asc_start_pts_coc_random.csv.addedpaths.txt" file
# included in the "unicor" folder is the output resistant kernel raster from
# UNICOR and can be read into R using raster::raster() or terra::rast()

# Download UNICOR here: https://github.com/ComputationalEcologyLab/UNICOR

