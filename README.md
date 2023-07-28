# Landscape genetics and connectivity analysis for grizzly bears in the Rocky Mountains in southeast British Columbia and southwest Alberta, Canada

This repository includes data and code for reproducing analyses in the following article, accepted for publication in *Molecular Ecology*:

**"Corridor-based approach with spatial cross validation reveals scale-dependent effects of geographic distance, human footprint, and canopy cover on grizzly bear genetic connectivity"**
by Eric Palm, Erin Landguth, Zachary Holden, Casey Day, Clayton Lamb, Paul Frame, Andrea Morehouse, Garth Mowat, Michael Proctor, Michael Sawaya, Gordon Stenhouse, Jesse Whittington, Katherine Zeller


Enclosed are data and code to do the steps outlined below. The actual GPS locations of grizzly bear genetic samples are not included at the request of data owners. Therefore, steps that require spatial locations (covariate extractions along straight lines and least cost paths between pairwise locations) use example datasets with randomly generated locations. All other steps use real data from the analysis.

Please open the `grizzly_lga.Rproj` file to start RStudio before opening individual code files to ensure that relative file paths in the code work correctly.

Here is a list of code files with descriptions:
1.	“create_pairwise_dataset.R” – Use a raw example dataset with animal IDs, real genotypes, and randomly generated spatial coordinates, and create a pairwise dataset for analyses.
2.	“extract_straight_example.R” – Use the example pairwise dataset created in Step 1 and extract mean covariate values along straight lines between pairwise locations.
3.	“model_straight_ffs.R” – Run a straight-line model with variable (feature) selection using the real dataset of extracted covariates (real spatial coordinates omitted) along straight lines, filtered to a maximum pairwise geographic distance of 40 km, and make a spatial prediction using only the retained variables.
4.	“plot_ALE.R” – Create an accumulated local effects plot for the retained environmental variables from Step 3 to visualize the relationships between genetic distance and those variables.
5.	“extract_lcp_example.R” – Use the example pairwise dataset created in Step 1 and extract covariate values along least cost paths between pairwise locations using the real spatial prediction raster from Step 3 as the resistance surface.
6.	“model_lcp_ffs_40_km” – Run a least cost path model with variable (feature) selection using the real dataset of extracted covariates (real spatial coordinates omitted) along least cost paths, filtered to a maximum pairwise geographic distance of 40 km, and make a spatial prediction using only the retained variables.
7.	“prepare_raster_UNICOR.R” – Prepare the final prediction (resistance) surface from the real full dataset model (which included all pairwise geographic distances) for running a resistant kernel connectivity model in UNICOR. Links to the UNICOR user manual are provided in this script.


Note: Models in this repository were run using covariate data extracted the same 360-m rasters included in the repository, rather than the 180-m rasters used in the manuscript, to speed up processing times for spatial predictions and UNICOR resistant kernel connectivity models. Therefore, the results from models in this repository differ slightly from those in the manuscript.


