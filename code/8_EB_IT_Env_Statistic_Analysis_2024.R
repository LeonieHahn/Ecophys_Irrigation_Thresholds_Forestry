# =============================================================================
# Analyse how the irrigation demand of F. sylvatica and Q. robur on EB in 2024
# relates to environmental drought conditions
# - Calculate rolling statistics for the VPD over different timesteps
# - Calculate Correlations between the irrigation demand and environemntal 
#   drought conditions for all tree species and treatments 
# - Perform linear regression between irrigation demand and environmental
#   drought conditions
# - Perform multiple regression and partial regression between irrigation demand 
#   and environmental drought conditions to investigate the influence of the 
#   single variable
# =============================================================================

library(readr)
library(ggplot2)
library(tidyr)
library(zoo)
library(purrr)
library(corrplot)
library(tidyverse)
library(Hmisc)
library(scales)
library(dplyr)
library(car)

source("./code/Functions.R")

irrig_thres_EB <- read_csv(
"./data/EB/results/Irrigation_thresholds/Irrigation_thres_EB_2024.csv", 
col_select = -1)

irrig_thres_EB$Date <- as.POSIXct(format(as.POSIXct(irrig_thres_EB$Date, 
                                                format="%Y-%m-%d %H:%M:%S", 
                                                tz = "Europe/Berlin"),
                                     format="%Y-%m-%d"))

# calculate rolling statistics for min, max and mean daily VPD,  
# for different timesteps, rename so the applied function (eg. min, mean) is 
# written at the end
columns <- c("VPD_kPA_mean", "VPD_kPA_min", "VPD_kPA_max")
windows <- c(seq(2, 7, 1)) 
funcs <- list(mean = mean, max = max, min = min) 

irrig_thres_EB <- calculate_rolling_stats(irrig_thres_EB, columns, windows, 
                                          funcs)
irrig_thres_EB <- calculate_rolling_stats(irrig_thres_EB, 
                                          c("Precip_mm_10min_sum", "SWC_mean"), 
                                          windows, 
                                          list(mean = mean))

# use only complete SWC observations and remove rain sheltered trees
irrig_thres_EB <- irrig_thres_EB %>% 
  filter(!is.na(SWC_mean) & RainShelter == "ns")

# Create empty dataframes for storing the results
final_cor_matrix <- data.frame()
final_p_matrix <- data.frame()

# Define different combinations of treatments and tree species
tree_species_list <- c("Quercus robur", "Fagus sylvatica", "Overall")
treatment_list <- c("Drip irrigated", "Unirrigated", "Overall")

# Calculate Correlations for all tree species and treatments
for (species in tree_species_list) {
  for (treatment in treatment_list) {
    
    # Special case: if tree species and treatment is "Overall" 
    if (species == "Overall" & treatment == "Overall") {
      correlation_data <- irrig_thres_EB %>%
        filter(Tree.Species %in% c("Quercus robur", "Fagus sylvatica"),
               Treatment %in% c("Drip irrigated", "Unirrigated")) %>%
        calculate_all_correlations_EB_2024(.)
      
      # Special case: if tree species is "Overall" 
    } else if (species == "Overall") {
      correlation_data <- irrig_thres_EB %>%
        filter(Tree.Species %in% c("Quercus robur", "Fagus sylvatica"),
               Treatment == treatment) %>%
        calculate_all_correlations_EB_2024(.)
      
      # Special case: if treatment is "Overall" 
    } else if (treatment == "Overall") {
      correlation_data <- irrig_thres_EB %>%
        filter(Tree.Species == species,
               Treatment %in% c("Drip irrigated", "Unirrigated")) %>%
        calculate_all_correlations_EB_2024(.)
      
      # Standard case: specific tree species and treatment
    } else {
      correlation_data <- irrig_thres_EB %>%
        filter(Tree.Species == species, Treatment == treatment) %>%
        calculate_all_correlations_EB_2024(.)
    }
    
    # Check, if enough data is available
    if (nrow(correlation_data) > 4) {
      cor_results <- calculate_correlations(correlation_data, species, 
                                            treatment)
      final_cor_matrix <- rbind(final_cor_matrix, cor_results$cor_matrix)
      final_p_matrix <- rbind(final_p_matrix, cor_results$p_matrix)
    } else {
      message(paste("Not enough observations for", species, "with treatment", 
                    treatment))
    }
  }
}

write.csv(final_cor_matrix, 
          "./data/EB/results/Irrigation_thresholds/Cor_EB_2024_all.csv")
write.csv(final_p_matrix, 
  "./data/EB/results/Irrigation_thresholds/Cor_EB_2024_all_p_value.csv")

final_cor_matrix$Variable <- rownames(final_cor_matrix)
final_p_matrix$Variable <- rownames(final_p_matrix)

final_cor_matrix$Variable <- gsub("(\\d{1,2})$", "", final_cor_matrix$Variable)
final_p_matrix$Variable <- gsub("(\\d{1,2})$", "", final_p_matrix$Variable)

cor_matrix_perc <- final_cor_matrix %>%
  select(Tree.Species, Treatment, Variable, Percentage_w_cros) %>%
  mutate(Percentage_w_cros_abs_cor = abs(Percentage_w_cros)) %>%
  filter(Variable != "Percentage_w_cros")

colnames(cor_matrix_perc)[4] <- "Percentage_w_cros_cor"

p_matrix_perc <- final_p_matrix %>%
  select(Tree.Species, Treatment, 
         Variable, Percentage_w_cros) %>%
  filter(Variable != "Percentage_w_cros") %>%
  mutate(p_value_cor_formatted = case_when(
    Percentage_w_cros < 0.001 ~ "< 0.001",
    Percentage_w_cros < 0.01  ~ "< 0.01",
    TRUE ~ formatC(Percentage_w_cros, format = "f", digits = 2)))
colnames(p_matrix_perc)[4] <- "p_value_cor"

cor_matrix_perc_sel <- cor_matrix_perc %>%
  mutate(variable_cat = case_when(
    str_detect(Variable, "Precip") ~ "Precip",
    str_detect(Variable, "VPD") ~ "VPD",
    str_detect(Variable, "SWC") ~ "SWC",
    TRUE ~ "Other")  
  ) %>%
  filter((variable_cat == "VPD" & Percentage_w_cros_cor > 0) | 
           (variable_cat %in% c("Precip", "SWC") & 
              Percentage_w_cros_cor < 0)) %>%  
  group_by(Treatment, Tree.Species, variable_cat) %>%
  slice_max(order_by = Percentage_w_cros_abs_cor, n = 1)

cor_matrix_perc_sel <- left_join(cor_matrix_perc_sel, p_matrix_perc, 
                                by = c("Tree.Species", "Treatment", "Variable"))

write.csv(cor_matrix_perc_sel, 
          "./data/EB/results/Irrigation_thresholds/Cor_EB_2024_best.csv")

cor_matrix_perc_sel_2023 <- read_csv(
  "./data/EB/results/Irrigation_thresholds/Cor_EB_2023_best.csv",
  col_select = -1)
cor_matrix_perc_sel_2023 <- cor_matrix_perc_sel_2023 %>%
  rename_with(~ paste0(.x, "_2023"), 
              .cols = !c("Tree.Species", "Treatment", 
                         "Variable", "variable_cat"))

cor_matrix_perc_sel <- cor_matrix_perc_sel %>%
  rename_with(~ paste0(.x, "_2024"), 
              .cols = !c("Tree.Species", "Treatment", "Variable", 
                         "variable_cat"))
cor_matrix_perc_sel_allyears <- full_join(cor_matrix_perc_sel_2023, 
                                          cor_matrix_perc_sel,
                                          by = c("Tree.Species", "Treatment", 
                                                 "Variable", "variable_cat"))

write.csv(cor_matrix_perc_sel_allyears, 
      "./data/EB/results/Irrigation_thresholds/Cor_EB_2023_2024_best.csv")

# (Simple) Regression analysis
results_reg <- data.frame(Tree.Species = character(),
                          Treatment = character(),
                          Variable = character(),
                          Adj_R2_reg = numeric(),
                          P_Value_reg = numeric(),
                          stringsAsFactors = FALSE)

apply(cor_matrix_perc_sel, 1, function(row) {
  species <- row["Tree.Species"]
  treatment <- row["Treatment"]
  variable <- as.character(row["Variable"])  
  percentage_column <- "Percentage_w_cros" 
  
  # if tree species or treatment are "Overall", 
  # replace with the corresponding values 
  if (species == "Overall") {
    species2 <- c("Fagus sylvatica", "Quercus robur")  # use both tree species
  } else {
    species2 <- as.character(species)  # otherwise just one tree species
  }
  
  if (treatment == "Overall") {
    treatment2 <- c("Unirrigated", "Drip Irrigated")   # use both treatments
  } else {
    treatment2 <- as.character(treatment)  # otherwise just one treatment
  }

  # perform regression
  df_filtered <- irrig_thres_EB %>%
    filter(Tree.Species %in% species2, Treatment %in% treatment2) %>%
    select(all_of(c(percentage_column, variable))) %>%
    na.omit() 
  
  # check if enough data is available for regression
  if (nrow(df_filtered) > 2) {
    model <- lm(as.formula(paste(percentage_column, "~", variable)), 
                data = df_filtered)
    adj_r2 <- summary(model)$adj.r.squared
    p_value <- summary(model)$coefficients[2, 4] 
    
    # combine results to dataframe
    results_reg <<- rbind(results_reg, data.frame(Tree.Species = species,
                                                  Treatment = treatment,
                                                  Variable = variable,
                                                  Adj_R2_reg = adj_r2,
                                                  P_Value_reg = p_value))
  } else {
    print("not enough data for regression")
  }
})

results_cor_reg <- left_join(cor_matrix_perc_sel, results_reg, 
                             by = c("Treatment", "Tree.Species", "Variable"))
results_cor_reg <- results_cor_reg %>%
  mutate(P_Value_reg_formatted = case_when(
    P_Value_reg < 0.001 ~ "< 0.001",
    P_Value_reg < 0.01  ~ "< 0.01",
    P_Value_reg < 0.05  ~ "< 0.05",
    TRUE ~ formatC(P_Value_reg, format = "f", digits = 3)),
    Percentage_w_cros_cor_2024 = formatC(Percentage_w_cros_cor_2024, 
                                    format = "f", digits = 2),
    Adj_R2_reg = formatC(Adj_R2_reg, format = "f", digits = 2),
    variable_cat = paste0(variable_cat, "_2024"))

write.csv(results_cor_reg,
      "./data/EB/results/Irrigation_thresholds/Cor+RegModels_EB_2024.csv",
          row.names = FALSE)

# add correlations results from EB 2023
results_cor_reg_2023 <- read_csv(
  "./data/EB/results/Irrigation_thresholds/Cor+RegModels_EB_2023.csv")
results_cor_reg_2023  <- results_cor_reg_2023  %>%
  mutate(variable_cat = paste0(variable_cat, "_2023")) %>%
  rename_with(~ paste0(.x, "_2023"), 
              .cols = !c("Tree.Species", "Treatment", "Variable", 
                         "variable_cat"))

results_cor_reg  <- results_cor_reg %>%
  rename_with(~ paste0(.x, "_2024"), 
              .cols = !c("Tree.Species", "Treatment", "Variable", 
                         "variable_cat"))
results_cor_reg_allyears <- full_join(results_cor_reg_2023, 
                                      results_cor_reg,
                                          by = c("Tree.Species", "Treatment", 
                                                 "Variable", "variable_cat"))

write.csv(results_cor_reg_allyears, 
  "./data/EB/results/Irrigation_thresholds/Cor+RegModels_EB_allyears.csv")

# multiple regression analysis
variable_combinations_df <- results_cor_reg %>%
  group_by(Tree.Species, Treatment) %>%
  summarise(Variables = list(unique(Variable)), .groups = "drop") %>%
  mutate(
    Variable_1 = map_chr(Variables, ~ .x[1] %||% NA),
    Variable_2 = map_chr(Variables, ~ .x[2] %||% NA),
    Variable_3 = map_chr(Variables, ~ .x[3] %||% NA)
  ) %>%
  select(Tree.Species, Treatment, Variable_1, Variable_2, Variable_3)

# Multiple regressions only for all treatments together, divided in tree species
# and all species combined

# Part 1: Fagus: Precip_mm_10min_sum_3_mean, SWC_mean, VPD_kPA_min_2_min
# Part 2: Quercus + Overall: Precip_mm_10min_sum_2_mean, SWC_mean, 
#         VPD_kPA_min_2_min

# calculate effect of drought variables on the irrigation demand 
# for F. sylvatica
part_1 <- irrig_thres_EB %>% 
  filter(!is.na(Precip_mm_10min_sum_3_mean) & !is.na(VPD_kPA_min_2_min) &
           !is.na(SWC_mean))

species_results_part_1 <- 
  evaluate_irrigation_models_by_species_part_1_2024(part_1)
species_results_part_1_sel <- filter(species_results_part_1, 
                                     Tree.Species == "Fagus sylvatica")

# Part 2: Quercus + Overall: Precip_mm_10min_sum_2_mean, SWC_mean, VPD_kPA_min_2_min
part_2 <- irrig_thres_EB %>% 
  filter(!is.na(Precip_mm_10min_sum_2_mean) & !is.na(VPD_kPA_min_2_min) &
           !is.na(SWC_mean))

# calculate effect of drought variables on the irrigation demand 
# for Q. robur
species_results_part_2 <- 
  evaluate_irrigation_models_by_species_part_2_2024(part_2)

overall_results_part_2 <- evaluate_irrigation_models_overall_part_2_2024(part_2)

results_part_2 <- rbind(species_results_part_2, overall_results_part_2)

species_results_part_2_sel <- filter(results_part_2, 
                                     Tree.Species == "Quercus robur" | 
                                       Tree.Species == "Overall")


mr_results_sel <- rbind(species_results_part_1_sel, species_results_part_2_sel)
mr_results_sel <- left_join(mr_results_sel, variable_combinations_df, 
                            by = c("Tree.Species", "Treatment"))

mr_results_sel <- mr_results_sel %>%
  mutate(across(contains("AR2"), ~ round(.x, 2))) %>%
  mutate(across(
    .cols = all_of(names(mr_results_sel)[startsWith(names(mr_results_sel), "p_")]),  
    .fns = ~ case_when(
      .x < 0.001 ~ "< 0.001",
      .x < 0.01  ~ "< 0.01",
      .x < 0.05  ~ "< 0.05",
      TRUE       ~ sprintf("%.3f", .x)
    )
  ))

write.csv(mr_results_sel,
          "./data/EB/results/Irrigation_thresholds/Part_Reg_EB_2024.csv",
          row.names = FALSE)
