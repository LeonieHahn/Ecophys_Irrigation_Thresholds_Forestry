# =============================================================================
# - Calculate irrigation thresholds ITlow and ITup per tree in the greenhouse 
# - Visualize the thresholds in combination with environmental measurements 
#   per tree and as overview heatmaps for all trees
# - Select max. irrigation demand for when the trees still recovered back to 
#   no irrigation demand (for the droughted trees only)
# =============================================================================

library(readr)
library(readxl)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(hrbrthemes)
library(tidyr)
library(ggtext)
library(zoo)
library(purrr)

# Load functions
source("./code/Functions.R")

model_results <- read_csv(
  "./data/GH/results/log_model_results/Log_model_results.csv", 
  col_select = -1)

threshold_perc <- read_csv(
  "./data/GH/results/log_model_results/p12thres_percentages.csv", 
  col_select = -1)

# table with output from DendRoAnalyst (phase.zg) and rel. log TwD
rel_TWD_df<-read_csv(
  "./data/GH/input_data/ZG_phase_rel_logTWD.csv",
  col_select = c(2:12))

rel_TWD_df$TIME <- sub(" UTC$", " CEST", rel_TWD_df$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
rel_TWD_df$TIME <- ifelse(nchar(rel_TWD_df$TIME) == 10,
                          paste(rel_TWD_df$TIME, "00:00:00"),
                          rel_TWD_df$TIME)
rel_TWD_df$TIME <- as.POSIXct(rel_TWD_df$TIME, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")
rel_TWD_df$Date <- as.Date(rel_TWD_df$TIME, tz = "Europe/Berlin")

# table with output from DendRoAnalyst (daily.data)
daily_df <- read_csv("./data/GH/input_data/daily_data_DF.csv", 
                     col_select = -1)
daily_df$DATE <- as.Date(daily_df$DATE, tz = "Europe/Berlin")
colnames(daily_df)[1] <- "Date"

# climate data
climate <- read_excel(
  "./data/GH/input_data/GH2023_T-Sensor-comparison.xlsx",
  sheet = 2)
colnames(climate) <- c("Date", "VPD_adj_dailyavg_kPa", "VPD_adj_dailymax_kPa",
                       "VPD_adj_dailymin_kPa", "n_dailyavg", 
                       "VPD_meas_dailyavg_kPa", "VPD_meas_dailymax_kPa",
                       "VPD_meas_dailymin_kPa")
climate$Date <- as.Date(climate$Date, tz = "Europe/Berlin")

# vol. soil water content 
vwc <- read_csv(
  "./data/GH/input_data/GH_VWC_dendrotrees.csv", 
  col_select = -1)
vwc$Tree.ID <- paste0(vwc$Tree.ID, "_", vwc$Treatment)
vwc$Date <- as.Date(vwc$Date, tz = "Europe/Berlin")

# table with TWD and canopy water potentials
combined <- read_csv("./data/GH/input_data/TWD_WP_combined.csv", 
                     col_select = -1)
combined$Date <- as.Date(combined$Date, tz = "Europe/Berlin")

# select best model for predicting TWD according to R²
best_model <- model_results %>%
  group_by(Tree.Species) %>%
  slice_max(order_by = R2_value, n = 1) %>%
  ungroup()

# select rel. TWD percentages (= irrigation demand) according to best R²
threshold_perc_sel <- threshold_perc %>%
  group_by(Tree.Species) %>%
  slice_max(order_by = R2_value, n = 1) %>%
  ungroup()

# select p12 thresholds according to tree species and best R²
thres_B <- filter(best_model, Tree.Species == "B")$p12_thres
thres_E <- filter(best_model, Tree.Species == "E")$p12_thres
thres_D <- filter(best_model, Tree.Species == "D")$p12_thres

# combine dataframes
combined <- left_join(rel_TWD_df, climate[, 1:4], by = "Date")
combined <- left_join(combined, vwc[, c(1:3)], by = c("Date", "Tree.ID")) 

# add min. rel TWD per day and Tree.ID 
combined  <- combined  %>%
  group_by(Tree.ID, Date) %>%
  mutate(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) %>%
  ungroup()

## calculate irrigation thresholds
# ITlow - calculate when min TWD > 0
MinTWD_not0 <- combined %>%
  group_by(Date = Date, Tree.ID, Treatment, Tree.Species) %>%
  summarise(Min_TWD = min(TWD)) %>%
  filter(Min_TWD != 0) 
MinTWD_not0$MinTWDnot0_Date <- as.Date(MinTWD_not0$Date)

# ITup - calculate when 12% Xylem conductivity loss corresponding min. TWD value 
# threshold is crossed
p12_thres_crossed <- combined %>%
  mutate(Date_rel_minTWD_thres_crossed = case_when(
    Tree.Species == "Fagus sylvatica" & min_rel_TWD >= thres_B ~ as.Date(Date),
    Tree.Species == "Quercus robur" & min_rel_TWD >= thres_E ~ as.Date(Date),
    Tree.Species == "Pseudotsuga menziesii" & 
      min_rel_TWD >= thres_D ~ as.Date(Date),
    TRUE ~ as.Date(NA)
  )) %>%
  select(Tree.ID, Date, Date_rel_minTWD_thres_crossed) %>%
  unique()

#  calculate last day of GRO and the day after last GRO
combined_sorted <- combined %>%
  arrange(Tree.ID, Date)
lastGRO <- combined_sorted %>%
  group_by(Tree.ID) %>%
  filter(Phases == "2") %>%
  slice_tail(n = 1)
lastGRO$FirstDayAfterLastGRO <- as.Date(lastGRO$Date) + 1
lastGRO$LastGRODate <- as.Date(lastGRO$Date)

FirstDayAfterLastGRO_df <- lastGRO[, c(6, 18, 18)]
colnames(FirstDayAfterLastGRO_df)[3] <- "Date"

# set "2023-09-07 and following dates to NA since this Date does not exist in 
# dataframe (data only exists until "2023-09-06")
FirstDayAfterLastGRO_df <- FirstDayAfterLastGRO_df %>%
  mutate(FirstDayAfterLastGRO = if_else(
    FirstDayAfterLastGRO >= as.Date("2023-09-07"), 
    as.Date(NA), 
    FirstDayAfterLastGRO
  ))

# merge results in one dataframe
results <- reduce(list(combined, MinTWD_not0[, c(1, 2, 6)], 
                       FirstDayAfterLastGRO_df,
                       p12_thres_crossed),
                  left_join, by = c("Date", "Tree.ID")) 

# fill in FirstDayAfterLastGRO all FirstDayAfterLastGRO-dates after 
# that occurred
results <- results %>%
  group_by(Tree.ID) %>%
  mutate(FirstDayAfterLastGRO = if_else(is.na(FirstDayAfterLastGRO), 
                                        lag(FirstDayAfterLastGRO), 
                                        FirstDayAfterLastGRO)) %>%
  fill(FirstDayAfterLastGRO, .direction = "down") %>%
  ungroup()

# Round threshold_perc_sel$rounded_rel_TWD and change
# add color to percentage value
threshold_perc_sel <- threshold_perc_sel %>%
  group_by(Tree.Species) %>%
  mutate(
    rounded_rel_TWD = round(p12_thres_scaled, 8) 
    )

# tree species and corresponding tree species abbreviation in df_base_perc
species_mapping <- tibble(
  results_species = c("Fagus sylvatica", "Quercus robur", 
                      "Pseudotsuga menziesii"),
  threshold_species = c("B", "E", "D")
)

# apply calc_thres_perc on all tree species and combine results
results <- species_mapping %>%
  pmap_dfr(function(results_species, threshold_species) {
    calc_thres_perc(
      df_results = results, 
      df_base_perc = threshold_perc_sel, 
      treespecies_l = results_species, 
      treespecies_s = threshold_species
    )
  })

results$Tree.Species.y <- NULL
colnames(results)[8] <- "Tree.Species"

# create new percentage column, where 0 = no threshold crossed,
# 1-99 threshold ITlow is crossed, 100 = threshold ITup is reached and 
# 101 = threshold ITup is crossed
results <- results %>%
  mutate(
    Percentage_w_cros = case_when(
      Percentage == 100 &
      !is.na("Date_rel_minTWD_thres_crossed") ~ 101, 
      TRUE ~ Percentage                               
    )
  )

results <- results %>%
  mutate(
    Percentage_color = case_when(
      Percentage_w_cros == 0 ~ "lightblue",      
      Percentage_w_cros == 101 ~ "darkmagenta",
      Percentage_w_cros > 0 & Percentage_w_cros <= 100 ~ scales::col_numeric(
        palette = c("goldenrod2", "mediumorchid3"),
        domain = c(1, 100)
      )(Percentage_w_cros),  # Continous color scheme for values from 1-99
      TRUE ~ "black"                  
    )
  )

# Create Demo color palette
palette_df <- data.frame(Percentage_w_cros = -10:110) %>%
  mutate(
  Percentage_color = case_when(
    Percentage_w_cros <= 0 ~ "lightblue",      
    Percentage_w_cros >= 101 ~ "darkmagenta",
    Percentage_w_cros > 0 & Percentage_w_cros <= 100 ~ scales::col_numeric(
      palette = c("goldenrod2", "mediumorchid3"),
      domain = c(1, 100)
    )(Percentage_w_cros),  # Continous color scheme for values from 1-99
    TRUE ~ "black"                  
  )
)

ggplot(palette_df, aes(x = Percentage_w_cros, fill = Percentage_color)) +
    geom_tile(aes(y = 1)) +
    scale_fill_identity() +
    scale_x_continuous(name = "%", 
                       breaks = c(-5, 1, 25, 50, 75, 100, 110),
                       labels = c( "0", "1", "25", "50", 
                                   "75", "100", ">100")) +
    theme_minimal() +
  theme(
      axis.ticks.y = element_blank(),       
    axis.text.y = element_blank(),       
    axis.title.y = element_blank(),       
    axis.text.x = element_text(size = 22), 
    axis.title.x = element_text(size = 24)  
  )

ggsave(path = "./graphics/GH/Irrigation_Thresholds/",
                filename = "Color_Pallete_Irrigation_Threshold_Percentages.png",
                width = 2560, height = 1489, units = "px")


# reverse plot
ggplot(palette_df, aes(x = Percentage_w_cros, fill = Percentage_color)) +
  geom_tile(aes(y = 1)) +
  scale_fill_identity() +
  scale_x_reverse( 
    name = "%", 
    breaks = c(110, 100, 75, 50, 25, 1, -5), 
    labels = c(">100", "100", "75", "50", "25", "1", "0")  
  ) +
  theme_minimal() +
  theme(
    axis.ticks.y = element_blank(),       
    axis.text.y = element_blank(),       
    axis.title.y = element_blank(),       
    axis.text.x = element_text(size = 22),  
    axis.title.x = element_text(size = 24) 
  )

ggsave(path =
         "./graphics/GH/Irrigation_Thresholds/",
       filename = "Color_Pallete_Irrigation_Threshold_Percentages_rev.png",
       width = 2560, height = 1489, units = "px")

# plot irrigation thresholds per tree
trees <- as.list(unique(results$Tree.ID))

for (i in trees){
  Tree.ID_Nr <- i
  print(Tree.ID_Nr)
  selected_tree_species <- substr(Tree.ID_Nr, 1, 1)
  rel_TWD_max_thres <- ifelse(selected_tree_species == "E", thres_E,
                              ifelse(selected_tree_species == "B", thres_B,
                                            thres_D))
  filtered_data <- filter(results, Tree.ID == Tree.ID_Nr)
  filtered_data$Date <- as.POSIXct(filtered_data$Date)
  plot_irrigation_thresholds(filtered_data, 
  "./graphics/GH/Irrigation_Thresholds/Irrigation_Thresholds_Singletree/")
}

# create stacked area charts with ggplot to visualize the dates when the 
# irrigation thresholds occurred
thresholds_daily <- results[, c(6:8, 12, 13:21, 23:24)] %>%
  group_by(Date, Tree.Species, Treatment)
thresholds_daily_unique <- thresholds_daily %>% distinct()

# calculate mean SWC, standard deviation and rel. SWC per treatment 
# and tree species for plotting later on
vwc_stats <- results[, c(7, 8, 12, 16)] %>%
  group_by(Date, Treatment, Tree.Species) %>%
  summarise(SWC_mean_treat = mean(VWC_mean, na.rm = TRUE),
            SWC_sd_treat = sd(VWC_mean, na.rm = TRUE))

# create new column with information that no threshold was reached 
# (but we would have had data to differentiate NAs from thresholds not reached)
thresholds_daily_unique <- thresholds_daily_unique %>%
  mutate(Date_nothres = case_when(
    is.na(MinTWDnot0_Date) & 
      is.na(FirstDayAfterLastGRO) & 
      is.na(Date_rel_minTWD_thres_crossed) ~  as.Date(Date),
    FALSE ~ as.Date(NA)))

# create heatmaps with irrigation thresholds per tree
heatmap_data <- thresholds_daily_unique %>%
  pivot_longer(
    cols = c(MinTWDnot0_Date,  
             Date_rel_minTWD_thres_crossed, 
             Date_nothres),
    names_to = "Threshold",
    values_to = "ThresholdDate"
  ) %>%
  filter(!is.na(ThresholdDate)) %>%  # Filter dates where threshold is crossed
  mutate(ThresholdDate = as.Date(ThresholdDate),  
         Date = as.Date(Date))                    
  
# Prioritize irrigation thresholds and keep only threshold with highest priority
heatmap_data_filtered <- heatmap_data %>%
  mutate(Threshold_Priority = case_when(
    Threshold == "Date_rel_minTWD_thres_crossed" ~ 2,
    Threshold == "MinTWDnot0_Date" ~ 1,
    Threshold == "Date_nothres" ~ 0
  )) %>%
  group_by(Tree.ID, Date) %>%
  slice_max(order_by = Threshold_Priority, with_ties = FALSE) %>%  
  ungroup()

heatmap_data_filtered$Date <- as.POSIXct(heatmap_data_filtered$Date, 
                                         format = "%Y-%m-%d")

# sort Tree.IDs descending after irrigation threshold percentages
threshold_summary <- heatmap_data_filtered %>%
  group_by(Tree.ID) %>%
  summarise(
    High_Percentage_Count = sum(Percentage_w_cros == "101", na.rm = TRUE),
    Mean_Percentage = mean(as.numeric(Percentage_w_cros), na.rm = TRUE)
  ) %>%
  arrange(desc(High_Percentage_Count), desc(Mean_Percentage))

heatmap_data_sorted <- heatmap_data_filtered %>%
  mutate(Tree.ID = factor(Tree.ID, levels = threshold_summary$Tree.ID)) %>%
  arrange(Tree.ID, Date)

FirstDayAfterLastGRO_df2 <- FirstDayAfterLastGRO_df %>%
  mutate(
    Tree.Species = case_when(
      grepl("B", Tree.ID) ~ "Fagus sylvatica",
      grepl("E", Tree.ID) ~ "Quercus robur",
      TRUE ~ "Pseudotsuga menziesii"),
    Treatment = ifelse(grepl("W", Tree.ID), "Watered", "Droughted")
  )

x_axis_limits <- c(as.POSIXct("2023-06-15"), as.POSIXct("2023-09-06"))
x_axis_settings <- scale_x_datetime(
  date_labels = "%d %b" ,
  limits = x_axis_limits
)

# Scaling factor for distributing SWC evenly on the 5 trees in the plot
swc_scaling_factor <- 0.1

heatmap_data_sorted <- heatmap_data_sorted %>%
  mutate(Tree.ID_num = as.numeric(Tree.ID)) %>%
  group_by(Tree.Species, Treatment) %>%
  mutate(Tree.ID_num_cat = as.numeric(factor(Tree.ID))) %>%
  ungroup()  

FirstDayAfterLastGRO_df2 <- left_join(FirstDayAfterLastGRO_df2,
                                      heatmap_data_sorted[, c(1, 13:18)], 
                                      by = "Tree.ID")

FirstDayAfterLastGRO_df2 <- unique(FirstDayAfterLastGRO_df2)


# Create Heatmaps with irrigation thresholds per tree and additional plots with
# the environmental conditions
VPD_plot_1 <- ggplot(filter(combined, Date >= "2023-06-19")) +
  geom_ribbon(aes(x = as.POSIXct(Date), 
                  ymin = VPD_adj_dailymin_kPa,
                  ymax = VPD_adj_dailymax_kPa), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = as.POSIXct(Date), y = VPD_adj_dailyavg_kPa), 
            color = "grey35", linewidth = 1, linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name = "VPD [kPa]",
    breaks = c(1, 3, 5),
    labels = c("1", "3", "5"),
    sec.axis = sec_axis(
      ~ .,  
      name = "")) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.y.left = element_text(size = 16),  
    axis.title.x = element_blank(),
    axis.text.y.left = element_text(size = 14),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  )


VPD_plot_2 <- ggplot(filter(combined, Date >= "2023-06-19")) +
  geom_ribbon(aes(x = as.POSIXct(Date), 
                  ymin = VPD_adj_dailymin_kPa,
                  ymax = VPD_adj_dailymax_kPa), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = as.POSIXct(Date), y = VPD_adj_dailyavg_kPa), 
            color = "grey35", linewidth = 1, linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name = "VPD [kPa]",
    breaks = c(1, 3, 5),
    labels = c("1", "3", "5"),
    sec.axis = sec_axis(
      ~ ., 
      name = "")) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.y.left = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y = element_blank()
    )

heatmap_plot_D <- ggplot() +
  geom_tile(
    data = filter(heatmap_data_sorted, Treatment == "Droughted"), 
    aes(x = as.POSIXct(Date), y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white") +
  scale_fill_identity() +  
  # add scaled SWC
  geom_ribbon(data = filter(vwc_stats, Treatment == "Droughted"),
              aes(x = as.POSIXct(Date),
                  ymin =(SWC_mean_treat - SWC_sd_treat) * swc_scaling_factor +1,
                  ymax = (SWC_mean_treat + SWC_sd_treat) * swc_scaling_factor +1),
              alpha = 0.5) +
  geom_line(
    data = filter(vwc_stats, Treatment == "Droughted"), 
    aes(
      x = as.POSIXct(Date), 
      y = SWC_mean_treat * swc_scaling_factor + 1,  
      group = 1
    ),  
    color = "black",
    linewidth = 1
  ) +
  # add FirstDayAfterLastGRO
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2, Treatment == "Droughted"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted, 
                            Treatment == "Droughted")$Tree.ID_num_cat),
    labels = 
      unique(filter(heatmap_data_sorted, 
                    Treatment == "Droughted")$Tree.ID_num_cat),
    sec.axis = sec_axis(
      ~ (. - 1) / swc_scaling_factor,  
      name = "SWC [%]",
      breaks = c(0, 20, 40, 60),
      labels = c(0,20, 40, 60)
    )
  ) +
  geom_label(
    data = filter(heatmap_data_sorted, 
                  Treatment == "Droughted"),
    aes(x = min(as.POSIXct(Date)) - (86400 * 5), y = Tree.ID_num_cat,
        label = gsub("_", "", Tree.ID)),
    color = "black",
    fill = "white",  
    size = 5,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),  
    label.r = unit(0.2, "lines") 
  ) +
  facet_grid(Tree.Species ~ Treatment, scales = "free", space = "free",
             switch = 'y',
             labeller = labeller(Tree.Species = as_labeller(c(
               "Fagus sylvatica" = "F. sylvatica",
               "Pseudotsuga menziesii" = "P. menziesii",
               "Quercus robur" = "Q. robur"
             )))) +
  theme_bw() +
  x_axis_settings +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.y.left = element_blank(), 
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.text.y.left = element_blank(),  
    axis.title.y.right = element_blank()
  )

print(heatmap_plot_D)

heatmap_plot_W <- ggplot() +
   geom_tile(
    data = filter(heatmap_data_sorted, Treatment == "Watered"), 
    aes(x = as.POSIXct(Date), y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white"
  ) +
  scale_fill_identity() +  
  # add scaled SWC
  geom_ribbon(data = filter(vwc_stats, Treatment == "Watered"),
              aes(x = as.POSIXct(Date),
               ymin =(SWC_mean_treat - SWC_sd_treat) * swc_scaling_factor + 1,
              ymax = (SWC_mean_treat + SWC_sd_treat) * swc_scaling_factor + 1),
              alpha = 0.5) +
  geom_line(
    data = filter(vwc_stats, Treatment == "Watered"), 
    aes(
      x = as.POSIXct(Date), 
      y = SWC_mean_treat * swc_scaling_factor + 1,  
      group = 1
    ),  
    color = "black",
    linewidth = 1
  ) +
  # add FirstDayAfterLastGRO
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2, Treatment == "Watered"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted, 
                            Treatment == "Watered")$Tree.ID_num_cat),
    labels = 
      unique(filter(heatmap_data_sorted, 
                    Treatment == "Watered")$Tree.ID_num_cat),
    sec.axis = sec_axis(
      ~ (. - 1) / swc_scaling_factor,  
      name = "SWC [%]",
      breaks = c(0, 20, 40, 60),
      labels = c(0, 20, 40, 60)
    )
  ) +
  geom_label(
    data = filter(heatmap_data_sorted, Treatment == "Watered"),
    aes(x = min(Date, na.rm = TRUE) - (86400 * 5), y = Tree.ID_num_cat,
        label = gsub("_", "", Tree.ID)),
    color = "black",
    fill = "white",  
    size = 5,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),  
    label.r = unit(0.2, "lines")  
    
  ) +
  facet_grid(Tree.Species ~ Treatment, scales = "free", space = "free",
             switch = 'y',
             labeller = labeller(Tree.Species = as_labeller(c(
               "Fagus sylvatica" = "F. sylvatica",
               "Pseudotsuga menziesii" = "P. menziesii",
               "Quercus robur" = "Q. robur"
             )))) +
  theme_bw() +
  x_axis_settings +
  theme(
    legend.position = "none",
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.y.right = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y.left = element_blank(),  
    axis.ticks.y.left = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y.left = element_blank(),  
    axis.title.y.right = element_text(size = 16, angle = 270, color = "black")
  )

print(heatmap_plot_W)

combined_plot_all_perc <- heatmap_plot_D /  heatmap_plot_W / 
  VPD_plot_1 / VPD_plot_2 +
  plot_layout(ncol = 2, nrow = 2, heights = c(1, 0.15), widths = c(1, 1))
combined_plot_all_perc

write.csv(heatmap_data_sorted, 
         "./data/GH/results/Irrigation_thresholds/Irrigation_thres_GH_2023.csv")

# find max. percentage_w_cros where each tree did still recover to an irrigation
# demand of 0
result_max_still_recovery <- heatmap_data_sorted %>%
  filter(Date >= "2023-06-15" & Date <= "2023-09-06",
         Treatment == "Droughted") %>%
  group_by(Tree.ID) %>%
  group_split() %>%
  lapply(find_threshold_change) %>%
  bind_rows()
