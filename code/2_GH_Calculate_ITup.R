# =============================================================================
# - Calculate logistic models for all combinations of predawn WP, midday WP,
#   min rel. TWD, max rel. TWD and plot results (see Dietrich et al., 2018)
# - Compare modeled rel. TWD with observations with a linear regression
#   (see Dietrich et al., 2018)
# - Calculate 90% confidence intervals with bootstrapping
#   (see Ziegler et al., 2024)
# - Select rel. TWD values for p12 (p12_thres) -> ITup
# - create plots with measurements, logistic models their confidence intervals
#   and derived irrigation threshold ITup
# =============================================================================

library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Load functions
source("./code/Functions.R")

# Psi12 values for the different tree species 
threshold_p12_B <- -2.5 # Fagus sylvatica
threshold_p12_E <- -3.51 # Quercus robur
threshold_p12_D <- -2.7 # Pseudotsuga menziesii

# Load table with min. and max. rel. TWD per Tree and corresponding predawn and 
# midday canopy waterpotential values (WP)
combined <- read_csv("./data/GH/input_data/TWD_WP_combined.csv", 
                     col_select = -1)


# create logistic model for predicting rel. TWD from WP, use all combinations
# of min. rel. TWD, max. rel. TWD, predawn WP and midday WP, 
# calculate confidence intervals with bootstrapping,
# select rel. TWD for p12 (= ITup) 

# create empty lists and dataframes to store results in
model_info <- list()
model_results <- data.frame()
wp_created_df <- data.frame()
model_info_comp_observ <- list()
p12_thres <- data.frame()
R2_df <- data.frame()
CI_df <- data.frame()

tree.species_list <- c("B", "E", "D")
log_pred_results <- lapply(
  tree.species_list, function(species) log_pred_twd_wp(
    data = combined, tree_species = species))

# bind results together
log_model_results <- do.call(rbind, lapply(log_pred_results, 
                                           function(x) x$model_results)) 

CI_results <- do.call(rbind, lapply(log_pred_results, 
                                    function(x) x$CI_df))
wp_created_results <- do.call(rbind, lapply(log_pred_results, 
                                            function(x) x$wp_created_df))
p12_thres_results <- do.call(rbind, lapply(log_pred_results, 
                                           function(x) x$p12_thres))
r2_lin_results <-do.call(rbind, lapply(log_pred_results, 
                                       function(x) x$R2_df))

parameter_combinations <- c("MinTWD_PDWP", "MaxTWD_PDWP", 
                            "MinTWD_MDWP", "MaxTWD_MDWP")
tree.species_list <- c("B", "E", "D")

# p-values
p_values_df <- data.frame(
  Tree.Species = character(),
  model = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(tree.species_list)) {
  for (param in parameter_combinations) {
    # model for Beech "MaxTWD_MDWP" could not be created, therefore write NA
    # for p-value
    if(i == 1 & param == "MaxTWD_MDWP") {
      p_value <- NA
    }
    else {
      coefficients <- 
        log_pred_results[[i]]$model_info_comp_observ[[param]]$coefficients
      # Extract p-value from results table
      p_value <- coefficients[2, "Pr(>|t|)"]
    }
    p_values_df <- rbind(p_values_df, data.frame(
      Tree.Species = tree.species_list[i],
      model = param,
      P_Value = p_value,
      stringsAsFactors = FALSE))
  }}

# change format of p-values
p_values_df <- p_values_df %>%
  mutate(P_Value_rounded = ifelse(P_Value < 0.001, "<0.001", 
                                  format(P_Value, scientific = TRUE)))

p12_thres_results <- p12_thres_results %>%
  pivot_longer(
    cols = starts_with("thres"), 
    names_to = "model",    
    values_to = "p12_thres"  
  ) %>%
  mutate(model = str_remove(model, "threshold_p12_"))

r2_lin_results_long <- r2_lin_results %>%
  pivot_longer(
    cols = starts_with("MinTWD") | starts_with("MaxTWD"), 
    names_to = "model",    
    values_to = "R2_value"  
  ) %>%
  mutate(model = str_remove(model, "_r2"))

model_results_combined <- reduce(list(r2_lin_results_long,
                                      p_values_df,
                                      p12_thres_results), 
                                 ~ left_join(.x, .y, 
                                             by = c("model", "Tree.Species")))

write.csv(model_results_combined, 
          "./data/GH/results/log_model_results/Log_model_results.csv")

# create percentages for rel TWD. approaching the p12 value (ITup)
threshold_perc <- data.frame(Percentage = seq(0, 100, 1))

threshold_perc <- model_results_combined %>%
  crossing(threshold_perc) %>%
  # Calculate p12_thres based on percentages
  mutate(p12_thres_scaled = case_when(
    Percentage == 0   ~ 0,  
    Percentage == 100 ~ p12_thres,  # use p12 thres when Percentage = 100 
    # fill values via linear interpolation
    TRUE              ~ p12_thres * (Percentage / 100)  
  )) %>%
  select(Tree.Species, model, Percentage, p12_thres_scaled)

# bind R² values
threshold_perc <- left_join(threshold_perc, model_results_combined[, 1:3],
                            by = c("Tree.Species", "model"))

write.csv(threshold_perc, 
          "./data/GH/results/log_model_results/p12thres_percentages.csv")

# combine all possible model combinations
model_combinations <- list(c("Min_rel_TWD", "WP_predawn"),
                           c("Max_rel_TWD", "WP_predawn"),
                           c("Min_rel_TWD", "WP_midday"),
                           c("Max_rel_TWD", "WP_midday"))

# plot logistic model results and derived irrigation thresholds
for (tree_species in tree.species_list) {
  for (i in 1:length(model_combinations)) {
    plot_log_models_p12_thres(model_data_df = log_model_results, 
                              CI_data = CI_results,
                              p12_data = p12_thres_results,
                              tree_species = tree_species,
                              w = wp_created_results, 
                              i = i,
                              logreg_path = 
                          "./graphics/GH/Irrigation_Thresholds/Log_Models_ITs/")
  }
}
