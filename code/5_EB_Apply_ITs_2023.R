# =============================================================================
# For the growing season 2023 on EB:
# - Combine dendrometer, climate and swc data
# - Calculate ITlow and ITup as well as the irrigation demand [%]
# - Create heatmaps for visualizing the irrigation demand in combination with 
#   environmental data (for fig. 5)
# =============================================================================

library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(hrbrthemes)
library(tidyr)
library(stringr)
library(lubridate)
library(forcats)
library(ggtext)

# Load functions
source("./code/Functions.R")

model_results <- read_csv(
  "./data/GH/results/log_model_results/Log_model_results.csv", 
  col_select = -1)

threshold_perc <- read_csv(
  "./data/GH/results/log_model_results/p12thres_percentages.csv", 
  col_select = -1)

# dendrometer data
rel_TWD_df <- read_csv("./data/EB/input_data/ZG_phase_rel_logTWD_2023.csv",
                       col_select = -1)
rel_TWD_df$TIME <- sub(" UTC$", " CEST", rel_TWD_df$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
rel_TWD_df$TIME <- ifelse(nchar(rel_TWD_df$TIME) == 10,
                       paste(rel_TWD_df$TIME, "00:00:00"),
                       rel_TWD_df$TIME)
rel_TWD_df$TIME <- as.POSIXct(rel_TWD_df$TIME, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")
rel_TWD_df$Date <- as.Date(rel_TWD_df$TIME, tz = "Europe/Berlin")
rel_TWD_df <- rel_TWD_df %>%
  filter(Treatment != "hose" & Treatment != "chip")

daily_df <- 
read_csv("./data/EB/input_data/daily_data_DF_2023.csv", 
                     col_select = -1)
daily_df$DATE <- as.Date(daily_df$DATE)
colnames(daily_df)[1] <- "Date"

daily_df <- daily_df %>%
  filter(Treatment != "hose" & Treatment != "chip")

# climate data
climate <- read_csv(
  "./data/EB/input_data/Meteo_2023_daily.csv",
  col_select = -1)

ggplot(climate)+
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/10), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, ymin = VPD_kPA_min, ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), color = "black", 
            linetype = "dashed") +   
  scale_y_continuous(
    name="Daily mean VPD [kPa]",  
    sec.axis = sec_axis(~.*10, name= "Precipitation [mm/day]") 
  ) +
  labs(x = "Date") +
  xlim(as.Date("2023-06-21"), as.Date("2023-10-01")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    axis.title.y.right = element_text(color = "dodgerblue3", size = 18),  
    axis.text.y.right = element_text(color = "dodgerblue3", size = 16),   
    axis.title.y.left = element_text(color = "black", size = 18),         
    axis.text.y.left = element_text(color = "black", size = 16)
  ) 

ggsave(path =
         "./graphics/EB/climate/2023/",
       filename = "Mean_VPD+Precipitation_2023.png",
       width = 2560, height = 1489, units = "px")

ggplot(climate)+
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/20), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, ymin = VPD_kPA_min, ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), color = "black", 
            linetype = "dashed") + 
  # add irrigation dates
  geom_vline(data = climate,
             aes(xintercept = as.Date("2023-07-11")),
             color = "chartreuse4", size = 1) +
  # 2023-07-12 only Eiche drip was irrigated for a second time, the other 
  # treatment did not receive irrigation
  geom_vline(data = climate,
             aes(xintercept = as.Date("2023-07-12")),
             color = "chartreuse4", size = 1) +
  geom_vline(data = climate,
             aes(xintercept = as.Date("2023-09-08")),
             color = "chartreuse4", size = 1) +
  geom_vline(data = climate,
             aes(xintercept = as.Date("2023-09-28")),
             color = "chartreuse4", size = 1) +
  scale_y_continuous(
    name="Daily mean VPD [kPa]",  
    sec.axis = sec_axis(~.*20, name= "Precipitation [mm/day]") 
  ) +
  labs(x = "Date") +
  xlim(as.Date("2023-06-21"), as.Date("2023-10-01")) +
  theme_bw() + 
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16),
    axis.title.y.right = element_text(color = "dodgerblue3", size = 18),  
    axis.text.y.right = element_text(color = "dodgerblue3", size = 16),   
    axis.title.y.left = element_text(color = "black", size = 18),         
    axis.text.y.left = element_text(color = "black", size = 16),
  ) 

ggsave(path =
         "./graphics/EB/climate/2023/",
       filename = "Mean_VPD+Precipitation+Irrig.dates.png",
       width = 2560, height = 1489, units = "px")


# vol. soil water content 
vwc_orig <- read_csv(
  "./data/EB/input_data/vSWC_EB_2023-2025_clean.csv",
  col_select = -1)

vwc_orig$TIME <- sub(" UTC$", " CEST", vwc_orig$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
vwc_orig$TIME <- ifelse(nchar(vwc_orig$TIME) == 10,
                          paste(vwc_orig$TIME, "00:00:00"),
                          vwc_orig$TIME)
vwc_orig$TIME <- as.POSIXct(vwc_orig$TIME, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")

vwc_data <- vwc_orig %>%
  pivot_longer(cols = 2: ncol(vwc_orig), values_to = "SWC", 
               names_to = "Sensor.ID") %>%
mutate(
    Treatment = case_when(
      str_detect(Sensor.ID, "drip") ~ "drip",
      str_detect(Sensor.ID, "control") ~ "control",
      str_detect(Sensor.ID, "hose") ~ "hose",
      str_detect(Sensor.ID, "chip") ~ "chip",
      TRUE ~ "unknown"
    ),
    Tree.Species = case_when(
      str_detect(Sensor.ID, "Buche") ~ "Buche",
      str_detect(Sensor.ID, "Eiche") ~ "Eiche",
      TRUE ~ "unknown"
    ),
    # RainShelter = if_else(str_detect(Sensor.ID, "_rs_"), "rs", "ns")
    Tree.ID = str_remove(Sensor.ID, "_vSWC"),
    Tree.ID = str_replace(Tree.ID, "rs", "ns")
  )

vwc_data <- vwc_data %>%
  filter(Treatment != "hose" & Treatment != "chip")

vwc_daily <- vwc_data %>%
  group_by(as.Date(TIME, tz = "Europe/Berlin"), Tree.ID) %>%
  summarise(SWC_mean = mean(SWC, na.rm = TRUE),
            SWC_min = min(SWC, na.rm = TRUE),
            SWC_max = max(SWC, na.rm = TRUE)) %>%
  mutate(
    Treatment = case_when(
      str_detect(Tree.ID, "drip") ~ "drip",
      str_detect(Tree.ID, "control") ~ "control",
      str_detect(Tree.ID, "hose") ~ "hose",
      str_detect(Tree.ID, "chip") ~ "chip",
      TRUE ~ "unknown"
    ),
    Tree.Species = case_when(
      str_detect(Tree.ID, "Buche") ~ "Buche",
      str_detect(Tree.ID, "Eiche") ~ "Eiche",
      TRUE ~ "unknown"
    ),
   RainShelter = if_else(str_detect(Tree.ID, "_rs_"), "rs", "ns")
  )
colnames(vwc_daily)[1] <- "Date"

vwc_daily <- vwc_daily %>%
  filter(Treatment != "hose" & Treatment != "chip")

# calculate mean VWC per treatment and tree species from daily sensor means
vwc_daily_treat <- vwc_daily %>%
  group_by(Date, Treatment, Tree.Species) %>%
  mutate(SWC_mean_treat = mean(SWC_mean, na.rm = TRUE),
         SWC_min_treat = min(SWC_mean, na.rm = TRUE),
         SWC_max_treat = max(SWC_mean, na.rm = TRUE),
         SWC_sd_treat = sd(SWC_mean, na.rm = TRUE))

ggplot(vwc_daily_treat) +
  geom_ribbon(aes(x = Date, ymin = SWC_min_treat, ymax = SWC_max_treat)) +
  geom_line(aes(x = Date, y = SWC_mean_treat)) +
  facet_grid(Tree.Species ~ RainShelter ~ Treatment)

ggplot(vwc_daily_treat %>%
         filter(Treatment != "hose")) +
  geom_ribbon(aes(x = Date, ymin = SWC_mean_treat - SWC_sd_treat, 
                  ymax = SWC_mean_treat + SWC_sd_treat), alpha = 0.4) +
  geom_line(aes(x = Date, y = SWC_mean_treat)) +
  facet_grid(Tree.Species ~ RainShelter ~ Treatment)

# select best model according to R²
best_model <- model_results %>%
  group_by(Tree.Species) %>%
  slice_max(order_by = R2_value, n = 1) %>%
  ungroup()

# select threshold percentages according to best model according to R²
threshold_perc_sel <- threshold_perc %>% 
  filter(Tree.Species == "B" | Tree.Species == "E") %>%
  group_by(Tree.Species) %>%
  slice_max(order_by = R2_value, n = 1) %>%
  ungroup()

# select p12 thresholds according to tree species and best R²
thres_B <- filter(best_model, Tree.Species == "B")$p12_thres
thres_E <- filter(best_model, Tree.Species == "E")$p12_thres

# combine dataframes
climate_daily <- filter(climate, 
                        Date >= as.Date("2023-06-03", tz = "Europe/Berlin"))

combined <- left_join(rel_TWD_df, climate[, c(1, 4:7)], 
                      by = "Date")
combined <- left_join(combined, vwc_daily_treat[, c(1:5, 9, 12)], 
                      by = c("Date", "Tree.ID"))

# subset dataframe until 01.10. since this will be the "critical" irrigation
# phase
combined <- filter(combined, Date <= as.Date("2023-10-01", 
                                             tz = "Europe/Berlin"))

# add min. rel TWD per day and Tree.ID 
combined  <- combined  %>%
  group_by(Tree.ID, Date) %>%
  mutate(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) %>%
  ungroup()

## calculate irrigation indicators

# ITlow - calculate when min TWD != 0
MinTWD_not0 <- combined %>%
  group_by(Date = Date, Tree.ID, Treatment, Tree.Species) %>%
  summarise(Min_TWD = min(TWD)) %>%
  filter(Min_TWD != 0) 
MinTWD_not0$MinTWDnot0_Date <- as.Date(MinTWD_not0$Date)

# ITup - calculate when 12% Xylem conductivity loss corresponding min. TWD value 
# threshold is crossed

p12_thres_crossed <- combined %>%
  group_by(Tree.Species, Date, Tree.ID) %>% 
  summarise(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) %>%
  ungroup()  %>%
  mutate(Date_rel_minTWD_thres_crossed = case_when(
    Tree.Species == "Fagus sylvatica" & min_rel_TWD >= thres_B ~ as.Date(Date),
    Tree.Species == "Quercus robur" & min_rel_TWD >= thres_E ~ as.Date(Date),
    TRUE ~ as.Date(NA)
  ))

# calculate last day of GRO and the day after last GRO
combined_sorted <- combined %>%
  arrange(Tree.ID, Date)
lastGRO <- combined_sorted %>%
  group_by(Tree.ID) %>%
  filter(Phases == "2") %>%
  slice_tail(n = 1)
lastGRO$FirstDayAfterLastGRO <- as.Date(lastGRO$Date) + 1
lastGRO$LastGRODate <- as.Date(lastGRO$Date)

FirstDayAfterLastGRO_df <- lastGRO[, c(6, 24, 24)]
colnames(FirstDayAfterLastGRO_df)[3] <- "Date"

# merge results in one dataframe
results <- reduce(list(combined, MinTWD_not0[, c(1, 2, 6)],  
                       FirstDayAfterLastGRO_df,
                       p12_thres_crossed[, c(2, 3, 5)]),
                  left_join, by = c("Date", "Tree.ID"))

results$Site <- "Eilaberg"

# Round threshold_perc_sel$rounded_rel_TWD and change
# add color to percentage value
threshold_perc_sel <- threshold_perc_sel %>%
  group_by(Tree.Species) %>%
  mutate(
    rounded_rel_TWD = round(p12_thres_scaled, 8) 
  )

# tree species and corresponding tree species abbreviation in df_base_perc
species_mapping <- tibble(
  results_species = c("Fagus sylvatica", "Quercus robur"),
  threshold_species = c("B", "E")
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
# 1-99 threshold A is crossed, 100 = threshold B is reached and 101 = threshold
# B is crosse
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

# list all Tree.IDs
trees <- as.list(unique(results$Tree.ID))

x_limits <- range(results$TIME)

for (i in trees){
  Site = "Eilaberg"
  Tree.ID_Nr <- i
  print(Tree.ID_Nr)
  filtered_data <- filter(results, Tree.ID == Tree.ID_Nr)
  selected_tree_species <- filtered_data$Tree.Species[[1]]
  rel_TWD_max_thres <- ifelse(selected_tree_species == 
                                "Quercus robur", thres_E,
                              ifelse(selected_tree_species == 
                                       "Fagus sylvatica", thres_B,
                                     ifelse(selected_tree_species == 
                                              "Abies alba", thres_T, 
                                            thres_D)))
  plot_irrigation_warnings(filtered_data, 
 "./graphics/EB/Irrigation_Thresholds/Irrigation_Thresholds_Singletree/2023/")
}
 
# select daily irrigation warnings
warnings_daily <- results[, c(6:8, 13:31)] %>%
  group_by(Date, Tree.Species, Treatment)
warnings_daily_unique <- warnings_daily %>% distinct()

climate_daily_EB <- filter(climate_daily, 
                             Date >= as.Date("2023-06-01", 
                                             tz = "Europe/Berlin"))

climate_daily_EB$Date <- as.POSIXct(climate_daily_EB$Date, tz = "Europe/Berlin")

# Create heatmaps
# create new column with information that no threshold was reached 
# (but we would have had data to differentiate NAs from thresholds not reached)
warnings_daily_unique <- warnings_daily_unique %>%
  mutate(Date_nothres = case_when(
    is.na(MinTWDnot0_Date) & 
      is.na(FirstDayAfterLastGRO) & 
      is.na(Date_rel_minTWD_thres_crossed)~as.Date(Date, tz = "Europe/Berlin"),
      FALSE ~ as.Date(NA)))

heatmap_data <- warnings_daily_unique %>%
  pivot_longer(
    cols = c(MinTWDnot0_Date,
             Date_rel_minTWD_thres_crossed, 
             Date_nothres),
    names_to = "Threshold",
    values_to = "ThresholdDate"
  ) %>%
  filter(!is.na(ThresholdDate)) %>%  
  mutate(ThresholdDate = as.Date(ThresholdDate, tz = "Europe/Berlin"),  
         Date = as.Date(Date, tz = "Europe/Berlin"))                   

# order thresholds after severity
heatmap_data_filtered <- heatmap_data %>%
  mutate(Threshold_Priority = case_when(
    Threshold == "Date_rel_minTWD_thres_crossed" ~ 3,
    Threshold == "MinTWDnot0_Date" ~ 2,
    Threshold == "Date_nothres" ~ 1
  )) %>%
  group_by(Tree.ID, Date) %>%
  # keep only the most severe threshold
  slice_max(order_by = Threshold_Priority, with_ties = FALSE) %>%  
  ungroup()

heatmap_data_filtered$Date <- as.POSIXct(heatmap_data_filtered$Date, 
                                         format = "%Y-%m-%d", 
                                         tz = "Europe/Berlin")

# sort descending after irrigation threshold percentages
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
      grepl("Buche", Tree.ID) ~ "Fagus sylvatica",
      TRUE ~ "Quercus robur"  
    ),
    Treatment = ifelse(grepl("control", Tree.ID), "control", 
                       ifelse(grepl("drip", Tree.ID), "drip",
                              ifelse(grepl("hose", Tree.ID), "hose", "chip"))),
    RainShelter = ifelse(grepl("rs", Tree.ID), "rs", "ns")
    
  )

x_axis_limits <- c(as.POSIXct("2023-05-25", tz = "Europe/Berlin"), 
                   as.POSIXct("2023-10-01", tz = "Europe/Berlin"))

# Define x-axis settings
x_axis_settings <- scale_x_datetime(
  date_labels = "%b",
  date_breaks = "1 month", 
  limits = x_axis_limits
)

combined$Date <- as.POSIXct(combined$Date)
climate$Date <- as.POSIXct(climate$Date)

VPD_plot_1 <- ggplot(filter(climate, 
                            Date >= as.POSIXct("2023-05-20", 
                                               tz = "Europe/Berlin"))) +
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/10), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, 
                  ymin = VPD_kPA_min,
                  ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), 
            color = "grey35", linewidth = 1, , linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name="VPD [kPa]",  
    breaks = c(0, 2, 4),
    sec.axis = sec_axis(~.*10, 
                        name= "P [mm/day]",
                        breaks = c(0, 20, 40))
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.y.left = element_line(color = "grey35"),
    axis.title.y.left = element_text(color = "grey35", size = 18),  
    axis.text.y.left = element_text(color = "grey35", size = 14)
  )


VPD_plot_2 <- ggplot(filter(climate, 
                            Date >= as.POSIXct("2023-05-20", 
                                               tz = "Europe/Berlin"))) +
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/10), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, 
                  ymin = VPD_kPA_min,
                  ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), 
            color = "grey35", linewidth = 1, linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name="VPD [kPa]",  
    breaks = c(0, 2, 4),
    sec.axis = sec_axis(~.*10, 
                        name= "P [mm/day]",
                        breaks = c(0, 20, 40)) 
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.y.left = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_line(color = "dodgerblue3"),
    axis.title.y.right = element_text(color = "dodgerblue3", size = 18),  
    axis.text.y.right = element_text(color = "dodgerblue3", size = 14)
  )


# Define scaling factor to distribute SWC over the 5 trees per treatment and
# species
vwc_scaling_factor <- 0.1

heatmap_data_sorted <- heatmap_data_sorted %>%
  mutate(Tree.ID_num = as.numeric(Tree.ID)) %>%
  group_by(Tree.Species, Treatment) %>%
  mutate(Tree.ID_num_cat = as.numeric(factor(Tree.ID))) %>%
  ungroup()  

FirstDayAfterLastGRO_df2 <- left_join(FirstDayAfterLastGRO_df2,
                                      heatmap_data_sorted[, c(1, 14, 15, 19, 25)], 
                                      by = c("Tree.ID", "FirstDayAfterLastGRO"))

FirstDayAfterLastGRO_df2 <- unique(FirstDayAfterLastGRO_df2)

# change order and labels of treatments
FirstDayAfterLastGRO_df2$Treatment <- FirstDayAfterLastGRO_df2$Treatment %>%
  fct_relevel("control", "drip") %>%  
  fct_recode(
    "Unirrigated" = "control",
    "Drip irrigated" = "drip"
  )

heatmap_data_sorted$Treatment <- heatmap_data_sorted$Treatment %>%
  fct_relevel("control", "drip") %>%  
  fct_recode(
    "Unirrigated" = "control",
    "Drip irrigated" = "drip")

vwc_daily_treat$Treatment <- vwc_daily_treat$Treatment %>%
  fct_relevel("control", "drip") %>%  
  fct_recode(
    "Unirrigated" = "control",
    "Drip irrigated" = "drip")

vwc_daily_treat$Tree.Species <- vwc_daily_treat$Tree.Species %>%
  fct_recode(
    "Fagus sylvatica" = "Buche",
    "Quercus robur" = "Eiche"
  )


heatmap_plot_B <- ggplot() +
  # Heatmap in the background
  geom_tile(
    data = filter(heatmap_data_sorted, Tree.Species == "Fagus sylvatica" &
                    Treatment != "hose"), 
    aes(x = Date, y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white"
  ) +
  scale_fill_identity() + 
  # add scaled SWC
  geom_ribbon(data = filter(heatmap_data_sorted, Treatment != "hose" &
                              Tree.Species == "Fagus sylvatica"),
             aes(x = as.POSIXct(Date),
             ymin =(SWC_mean_treat - SWC_sd_treat) * vwc_scaling_factor + 1,
             ymax = (SWC_mean_treat + SWC_sd_treat) * vwc_scaling_factor + 1),
             alpha = 0.5) +
  geom_line(
    data = filter(heatmap_data_sorted, Treatment != "hose" &
                    Tree.Species == "Fagus sylvatica"),
    aes(
      x = as.POSIXct(Date),
      y = SWC_mean_treat * vwc_scaling_factor + 1,  # adjust scaling
      group = 1
    ),
    color = "black",
    linewidth = 1
  ) +
  # add FirstDayAfterLastGRO 
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2, Tree.Species == "Fagus sylvatica" &
                    Treatment != "hose"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted, 
                            Tree.Species == "Fagus sylvatica" &
                              Treatment != "hose")$Tree.ID_num_cat),
    labels = 
      # Label after Tree.ID_num and Tree.ID 
      unique(filter(heatmap_data_sorted, 
                    Tree.Species == "Fagus sylvatica")$Tree.ID_num_cat)
   ) +
  geom_label(
    data = filter(heatmap_data_sorted, 
                  Tree.Species == "Fagus sylvatica" & Treatment != "hose") %>%
      mutate(
        label_short = paste0(
          case_when(
            grepl("Buche", Tree.ID) ~ "B",
            grepl("Eiche", Tree.ID) ~ "E",
            TRUE ~ "X"  
          ),
          gsub("^D(\\d+)_.*", "\\1", Tree.ID),
          case_when(
            grepl("control", Tree.ID) ~ "U",
            grepl("drip", Tree.ID) ~ "I",
            TRUE ~ "X"
          )
        )
      ),
    aes(
      x = as.POSIXct("2023-05-25"),
      y = Tree.ID_num_cat,
      label = label_short
    ),
    color = "black",
    fill = "white",
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.2, "lines")
  ) +
  facet_grid(Treatment ~ Tree.Species, scales = "free", space = "free",
             switch = 'y',
             labeller = labeller(Tree.Species = as_labeller(c(
               "Fagus sylvatica" = "F. sylvatica",
               "Quercus robur" = "Q. robur"
             )))) +
  theme_bw() +
  x_axis_settings +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 16),
    strip.text.x = element_text(size = 16),
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

print(heatmap_plot_B)

heatmap_plot_E <- ggplot() +
  geom_tile(
    data = filter(heatmap_data_sorted, Tree.Species == "Quercus robur" &
                    Treatment != "hose"), 
    aes(x = Date, y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white"
  ) +
  scale_fill_identity() +
  geom_ribbon(data = filter(heatmap_data_sorted, Treatment != "hose" &
                              Tree.Species == "Quercus robur"),
              aes(x = as.POSIXct(Date),
               ymin =(SWC_mean_treat - SWC_sd_treat) * vwc_scaling_factor + 1,
              ymax = (SWC_mean_treat + SWC_sd_treat) * vwc_scaling_factor + 1),
              alpha = 0.5) +
  geom_line(
    data = filter(heatmap_data_sorted, Treatment != "hose" &
                    Tree.Species == "Quercus robur"),
    aes(
      x = as.POSIXct(Date),
      y = SWC_mean_treat * vwc_scaling_factor + 1,  
      group = 1
    ),
    color = "black",
    linewidth = 1
  ) +
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2, Tree.Species == "Quercus robur" &
                    Treatment != "hose"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted, 
                            Tree.Species == "Quercus robur" & 
                              Treatment != "hose")$Tree.ID_num_cat),
    labels = 
      unique(filter(heatmap_data_sorted, 
                    Tree.Species == "Quercus robur" & 
                      Treatment != "hose")$Tree.ID_num_cat),
    sec.axis = sec_axis(
      ~ (. - 1) / vwc_scaling_factor,  
      name = "SWC [%]"
       )
  ) +
  geom_label(
    data = filter(heatmap_data_sorted, 
                  Tree.Species == "Quercus robur" & Treatment != "hose") %>%
      mutate(
        label_short = paste0(
          case_when(
            grepl("Buche", Tree.ID) ~ "B",
            grepl("Eiche", Tree.ID) ~ "E",
            TRUE ~ "X"  
          ),
          gsub("^D(\\d+)_.*", "\\1", Tree.ID),
          case_when(
            grepl("control", Tree.ID) ~ "U",
            grepl("drip", Tree.ID) ~ "I",
            TRUE ~ "X"
          )
        )
      ),
    aes(
      x = as.POSIXct("2023-05-25"),
      y = Tree.ID_num_cat,
      label = label_short
    ),
    color = "black",
    fill = "white",
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.2, "lines")
  ) +
  facet_grid(Treatment ~ Tree.Species, scales = "free", space = "free",
             switch = 'y',
             labeller = labeller(Tree.Species = as_labeller(c(
               "Fagus sylvatica" = "F. sylvatica",
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

heatmap_plot_E 

combined_plot_all <- heatmap_plot_B /  heatmap_plot_E / VPD_plot_1 / VPD_plot_2 +
  plot_layout(ncol = 2, nrow = 2, heights = c(1, 0.15), widths = c(1, 1))

combined_plot_all

write.csv(heatmap_data_sorted, 
    "./data/EB/results/Irrigation_thresholds/Irrigation_thres_EB_2023.csv")

write.csv(results, 
 "./data/EB/results/Irrigation_thresholds/Irrigation_results_EB_2023.csv")
