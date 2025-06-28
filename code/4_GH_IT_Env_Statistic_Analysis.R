# =============================================================================
# - calculate rolling mean, min. and max. values for different environmental
#   drought stress variables
# - perform correlation analysis between irrigation demand and environmental
#   drought stress variables
# - create scatterplots for the highest correlating variables
# - perform linear, multiple and partial regression between irrigation demand 
#   and environmental drought stress variables
# =============================================================================

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(zoo)
library(purrr)
library(corrplot)
library(tidyverse)
library(Hmisc)
library(readxl)
library(car)

source("./code/Functions.R")

irrig_thres_GH <- read_csv(
  "./data/GH/results/Irrigation_thresholds/Irrigation_thres_GH_2023.csv",
  col_select = -1)

colnames(irrig_thres_GH)[8] <- "SWC"

irrig_thres_GH <- filter(irrig_thres_GH, Date < as.Date("2023-09-07")) 

# calculate rolling statistics for min, max and mean daily VPD 
# for different timesteps
columns <- c("VPD_adj_dailymin_kPa", "VPD_adj_dailymax_kPa",
             "VPD_adj_dailyavg_kPa")  
windows <- c(seq(2, 7, 1)) 
funcs <- list(mean = mean, max = max, min = min) 

irrig_thres_GH <- calculate_rolling_stats(
  irrig_thres_GH, columns, windows, funcs)
irrig_thres_GH <- calculate_rolling_stats(irrig_thres_GH, 
                                          "SWC", 
                                          windows, 
                                          list(mean = mean))

# calculate correlations between irrigation demand and environmental drought 
# stress variables
final_cor_matrix <- data.frame()
final_p_matrix <- data.frame()

# Create combinations between tree species and treatments
tree_species_list <- c("Quercus robur", "Fagus sylvatica", 
                       "Pseudotsuga menziesii", "Overall")
treatment_list <- c("Droughted", "Watered", "Overall")

# Calculate correlations for all tree species and treatments
for (species in tree_species_list) {
  for (treatment in treatment_list) {
    
    # Exception: tree species and treatment are both "Overall" 
    if (species == "Overall" & treatment == "Overall") {
      correlation_data <- irrig_thres_GH %>%
        filter(Tree.Species %in% c("Quercus robur", "Fagus sylvatica",
                                   "Pseudotsuga menziesii"),
               Treatment %in% c("Droughted", "Watered")) %>%
        calculate_all_correlations(.)
      
      # if tree species is "Overall", but the treamtent is specific
    } else if (species == "Overall") {
      correlation_data <- irrig_thres_GH %>%
        filter(Tree.Species %in% c("Quercus robur", "Fagus sylvatica",
                                   "Pseudotsuga menziesii"),
               Treatment == treatment) %>%
        calculate_all_correlations(.)
      
      # if treatment is "Overall", but the tree species is specific
    } else if (treatment == "Overall") {
      correlation_data <- irrig_thres_GH %>%
        filter(Tree.Species == species,
               Treatment %in% c("Droughted", "Watered")) %>%
        calculate_all_correlations(.)
      
      # Standard case: tree species and treatment are given
    } else {
      correlation_data <- irrig_thres_GH %>%
        filter(Tree.Species == species, Treatment == treatment) %>%
        calculate_all_correlations(.)
    }
    
    # check if enough observations exist
    if (nrow(correlation_data) > 4) {
      
      # calculate correlations and p-values
      cor_results <- calculate_correlations(correlation_data, species, 
                                            treatment)
      
      # save results
      final_cor_matrix <- rbind(final_cor_matrix, cor_results$cor_matrix)
      final_p_matrix <- rbind(final_p_matrix, cor_results$p_matrix)
      
    } else {
      message(paste("Not enough observations for", species, 
                    "with treatment", treatment))
    }
  }
}

write.csv(final_cor_matrix, 
          "./data/GH/results/Irrigation_thresholds/Cor_GH_2023.csv")
write.csv(final_p_matrix, 
          "./data/GH/results/Irrigation_thresholds/Cor_GH_2023_p_value.csv")

final_cor_matrix$Variable <- rownames(final_cor_matrix)
final_p_matrix$Variable <- rownames(final_p_matrix)

final_cor_matrix$Variable <- gsub("(\\d{1,2})$", "", final_cor_matrix$Variable)
final_p_matrix$Variable <- gsub("(\\d{1,2})$", "", final_p_matrix$Variable)

# select only correlations with the irrigation demand
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

cor_matrix_perc <- left_join(cor_matrix_perc, p_matrix_perc, 
                              by = c("Tree.Species", "Treatment", "Variable"))

# select only significant negative correlations between irrigation demand and
# SWC and significant positive correlations with VPD
cor_matrix_perc_sel <- cor_matrix_perc %>%
  mutate(variable_cat = case_when(
    str_detect(Variable, "VPD") ~ "VPD",
    str_detect(Variable, "SWC") ~ "SWC",
    TRUE ~ "Other")  
  ) %>%
  filter(p_value_cor < 0.01) %>%
  filter((variable_cat == "VPD" & Percentage_w_cros_cor > 0) | 
           (variable_cat == "SWC" & Percentage_w_cros_cor < 0)) %>%  
  group_by(Treatment, Tree.Species, variable_cat) %>%
  slice_max(order_by = Percentage_w_cros_abs_cor, n = 1)

write.csv(cor_matrix_perc_sel, 
          "./data/GH/results/Irrigation_thresholds/Cor_GH_2023_best.csv")

# Scatterplot with highest correlation variables
ggplot(data = irrig_thres_GH,
       aes(x = SWC_5_mean,
           y = Percentage_w_cros)) +
  geom_point(aes(fill = Percentage_color), shape = 21, color = "black",
             size = 2.5) +
  scale_fill_identity() +
  facet_grid(Treatment ~ Tree.Species) +
  theme_bw()+
  ylab("Irrigation demand [%]") + xlab("SWC [%]")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16))
ggsave(path =
       "./graphics/GH/Irrigation_Thresholds/Correlation_Analysis/",
       filename = "Scatterplot_SWC_5mean_Perc.png",
       width = 2560, height = 1489, units = "px")

ggplot(data = irrig_thres_GH,
       aes(x = VPD_adj_dailyavg_kPa_7_max,
           y = Percentage_w_cros )) +
  geom_point(aes(fill = Percentage_color), shape = 21, color = "black",
             size = 2.5) +
  scale_fill_identity() +
  facet_grid(Treatment ~ Tree.Species) +
  theme_bw()+
  xlab("Rolling max (7 days) over daily mean. VPD [kPa]") + 
  ylab("Irrigation demand [%]")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom",
        axis.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16))
ggsave(path =
       "./graphics/GH/Irrigation_Thresholds/Correlation_Analysis/",
       filename = "Scatterplot_VPDmean_7_max_Perc.png",
       width = 2560, height = 1489, units = "px")

# calculate linear regression between irrigation demand and environmental 
# drought stress variables
results_reg <- data.frame(Tree.Species = character(),
                      Treatment = character(),
                      Variable = character(),
                      Adj_R2_reg = numeric(),
                      P_Value_reg = numeric(),
                      stringsAsFactors = FALSE)

apply(cor_matrix_perc_sel, 1, function(row) {
  species <- row["Tree.Species"]
  treatment <- row["Treatment"]
  variable <- as.character(row["Variable"])  # independent regression variable 
  percentage_column <- "Percentage_w_cros"  # dependent regression variable 

  # If Tree.Species or treatment is "Overall", replace with corresponding values
  # of all tree species
  if (species == "Overall") {
    species2 <- c("Fagus sylvatica", "Quercus robur", "Pseudotsuga menziesii")  
  } else {
    species2 <- as.character(species)  
  }

  if (treatment == "Overall") {
    treatment2 <- c("Droughted", "Watered")  
  } else {
    treatment2 <- as.character(treatment) 
  }

        # If the combination is unique, perform regression
        df_filtered <- irrig_thres_GH %>%
          filter(Tree.Species %in% species2, Treatment %in% treatment2) %>%
          select(all_of(c(percentage_column, variable))) %>%
          na.omit() 

        # check if enough data is available
        if (nrow(df_filtered) > 2) {
          model <- lm(as.formula(paste(percentage_column, "~", variable)), 
                      data = df_filtered)
          adj_r2 <- summary(model)$adj.r.squared
          p_value <- summary(model)$coefficients[2, 4]  

          # create dataframe with results
          results_reg <<- rbind(results_reg, data.frame(Tree.Species = species,
                                                Treatment = treatment,
                                                Variable = variable,
                                                Adj_R2_reg = adj_r2,
                                                P_Value_reg = p_value))
        } else {
          print("Not enough data available for regression")
        }
})

# combine results from correlation and regression analysis
results_cor_reg <- left_join(cor_matrix_perc_sel, results_reg, 
                             by = c("Treatment", "Tree.Species", "Variable"))
results_cor_reg <- results_cor_reg %>%
mutate(P_Value_reg_formatted = case_when(
  P_Value_reg < 0.001 ~ "< 0.001",
  P_Value_reg < 0.01  ~ "< 0.01",
  TRUE ~ formatC(P_Value_reg, format = "f", digits = 3)),
  Percentage_w_cros_cor = formatC(Percentage_w_cros_cor, 
                                  format = "f", digits = 2),
  Adj_R2_reg = formatC(Adj_R2_reg, format = "f", digits = 2))

write.csv(results_cor_reg,
          "./data/GH/results/Irrigation_thresholds/Cor+RegModels_GH_2023.csv",
          row.names = FALSE)

# only Quercus robur watered had 2 variables correlating, therefore perform a
# partial regression only for this case

# use only complete VPD_adj_dailymax_kPa_7_max observations
irrig_thres_GH <- irrig_thres_GH %>% 
  filter(Tree.Species == "Quercus robur" & Treatment == "Watered" &
           !is.na(SWC_5_mean) & !is.na(VPD_adj_dailyavg_kPa_7_max))

# calculate correlation and regression between irrigation demand and 
# environmental drought stress variables
species_treatment_results <- 
  evaluate_irrigation_models_by_species_and_treatment(irrig_thres_GH)

species_treatment_results <-species_treatment_results %>%
  mutate(across(contains("AR2"), ~ round(.x, 2))) %>% 
  mutate(across(contains("Extra"), ~ round(.x, 2))) %>% 
  mutate(across(contains("p_"), ~ case_when(
    .x < 0.001 ~ "< 0.001",
    .x < 0.01  ~ "< 0.01",
    .x < 0.05  ~ "< 0.05",
    TRUE       ~ sprintf("%.3f", .x)
  )))

write.csv(species_treatment_results,
"./data/GH/results/Irrigation_thresholds/RegModels_GH_2023_SWC_5_mean_VPD_avg7max_EW.csv",
           row.names = FALSE)
