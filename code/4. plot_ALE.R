
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

# Script for creating an accumulated local effects (ALE) plot from the variable-
# selected straight-line model from the "3. model_straight_ffs_40_km.R" script.

# Code written by Eric Palm

################################################################################

# Load packages
require(caret)
require(gbm)
require(iml)
require(dplyr)
require(ggplot2)

# Function to read a saved model object fit using caret::train, and
# use the iml package to generate accumulated locals effects for plotting
# relationships between genetic distance and environmental covariates
get_ale <- function(x) {
 
   gbm <- readRDS(x)
  
   pred_gbm <- 
    iml::Predictor$new(gbm,
                       gbm$trainingData %>% dplyr::select(-.outcome),
                       y = gbm$trainingData$.outcome)
   eff_gbm <- 
     iml::FeatureEffects$new(pred_gbm, grid.size = 100)
   
   ale <- eff_gbm$results %>%
     dplyr::bind_rows() %>% 
     dplyr::rename(., x = ".borders", y = ".value", covariate = ".feature") %>% 
     dplyr::as_tibble() 
   
   return(ale)
}

# Run the get_ale function on the straight-line model from the 40-km dataset 
ale <- "results/real/model_straight_40_km.rds" %>% 
  get_ale()

# Plot the ALE curves faceted by covariate
ale %>% 
  dplyr::mutate(covariate = case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                               TRUE ~ "Paved roads"),
         x = dplyr::if_else(covariate == "Geographic distance (km)", x/1000, x)) %>% 
  ggplot(., aes(x = x, y = y)) +
  facet_wrap(~ covariate, scales = "free", strip.position = "bottom") +
  geom_line(linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
  theme_classic(base_size = 15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        strip.placement = "outside",
        strip.text = element_text(size = 15)) +
  labs(y = "Predicted change from mean genetic distance \n(i.e., landscape resistance)",
       x = NULL) 

# Save the figure to a file
ggsave("figures/ALE_plot_straight_model_40_km.png", dpi = "print",  height = 9, width = 10.5)

