# =============================================================================
# For the growing season 2024 on EB:
# - Combine dendrometer, climate and swc data
# - Create overview plot of SWC from the rain-sheltered sensors (S5)
# - Calculate ITlow and ITup as well as the irrigation demand [%]
# - Create heatmaps for visualizing the irrigation demand in combination with 
#   environmental data (for fig. 5)
# - Explore irrigation effect on soil moisture, TWD and irrigation demand
#   (fig. 6, S6)
# =============================================================================

library(readr)
library(readxl)
library(purrr)
library(ggplot2)
library(patchwork)
library(hrbrthemes)
library(tidyr)
library(stringr)
library(forcats)
library(ggtext)
library(lubridate)
library(dplyr)

# Load functions
source("./code/Functions.R")

model_results <- read_csv(
  "./data/GH/results/log_model_results/Log_model_results.csv", 
  col_select = -1)

threshold_perc <- read_csv(
  "./data/GH/results/log_model_results/p12thres_percentages.csv", 
  col_select = -1)

# dendrometer data
rel_TWD_df <- read_csv("./data/EB/input_data/ZG_phase_rel_logTWD_2024.csv",
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

daily_df <- read_csv("./data/EB/input_data/daily_data_DF_2024.csv", 
                     col_select = -1)
daily_df$DATE <- as.Date(daily_df$DATE)
colnames(daily_df)[1] <- "Date"

daily_df  <- daily_df %>%
  filter(Tree.Species == "Fagus sylvatica" & Treatment != "hose" & 
           Treatment != "chip" |
           Tree.Species == "Quercus robur" & Treatment != "hose" & 
           Treatment != "chip")

# inforamtion about rainshelter installation
rs_meta <- read_csv("./data/EB/input_data/Overview_Rainshelters2024.csv")

rs_meta$Installation_Date <- as.POSIXct(rs_meta$Installation_Date,
                                         format = "%m/%d/%Y %H:%M",
                                         tz = "Europe/Berlin")
rs_meta$Deinstallation_Date <- as.POSIXct(rs_meta$Deinstallation_Date,
                                         format = "%m/%d/%Y %H:%M",
                                         tz = "Europe/Berlin")

# climate data
climate <- read_csv("./data/EB/input_data/Meteo_2024_daily.csv",
  col_select = -1)

climate <- filter(climate, 
                  Date > "2024-05-01" & Date < "2024-10-27") 

ggplot(climate)+
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/20), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, ymin = VPD_kPA_min, ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), color = "black", 
            linetype = "dashed") +   
  scale_y_continuous(
    name="Daily mean VPD [kPa]",  
    sec.axis = sec_axis(~.*20, name= "Precipitation [mm/day]") 
  )+
  labs(x = "Date") +
  xlim(as.Date("2024-05-01"), as.Date("2024-10-27"))+
  theme_bw()+ 
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
         "./graphics/EB/climate/2024/",
       filename = "Mean_VPD+Precipitation_2024.png",
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
             aes(xintercept = as.Date("2024-08-27")),
             color = "chartreuse4", linewidth = 1) +
  
  # only pedunculate oak
  geom_vline(data = climate,
             aes(xintercept = as.Date("2024-09-05")),
             color = "chartreuse4", linewidth = 1) +
  scale_y_continuous(
    name="Daily mean VPD [kPa]",  
    sec.axis = sec_axis(~.*20, name= "Precipitation [mm/day]") 
  ) +
  labs(x = "Date") +
  xlim(as.Date("2024-05-01"), as.Date("2024-10-27"))+
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
         "./graphics/EB/climate/2024/",
       filename = "Mean_VPD+Precipitation_2024+Irrig.dates.png",
       width = 2560, height = 1489, units = "px")


# vol. soil water content 
vwc_orig <- read_csv("./data/EB/input_data/vSWC_EB_2023-2025_clean.csv",
  col_select = -1)

# reformat date to CEST
vwc_orig$TIME <- sub(" UTC$", " CEST", vwc_orig$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
vwc_orig$TIME <- ifelse(nchar(vwc_orig$TIME) == 10,
                       paste(vwc_orig$TIME, "00:00:00"),
                       vwc_orig$TIME)
vwc_orig$TIME <- as.POSIXct(vwc_orig$TIME,
                           format="%Y-%m-%d %H:%M:%S",
                           tz = "Europe/Berlin")
vwc_orig <- filter(vwc_orig,
                  TIME >= as.POSIXct("2024-05-01 00:00:00",
                                              format="%Y-%m-%d %H:%M:%S",
                                              tz = "Europe/Berlin") &
                    TIME <= as.POSIXct("2024-10-01 01:00:00",
                                       format="%Y-%m-%d %H:%M:%S",
                                       tz = "Europe/Berlin"))

vwc_orig <- vwc_orig %>%
  select(TIME, contains("drip"), contains("control"))

vwc_orig <- vwc_orig %>%
  rename(
    D100_Buche_drip_ns_vSWC = D25_Buche_drip_ns_vSWC,
    D101_Buche_drip_ns_vSWC = D26_Buche_drip_ns_vSWC, # soil sensor not under rainshield, but tree
    D102_Buche_drip_ns_vSWC = D27_Buche_drip_ns_vSWC, # soil sensor not under rainshield, but tree
    D103_Buche_drip_ns_vSWC = D29_Buche_drip_ns_vSWC,
    D104_Buche_drip_ns_vSWC = D30_Buche_drip_ns_vSWC,
    D105_Buche_control_ns_vSWC = D32_Buche_control_ns_vSWC,
    D113_Eiche_drip_ns_vSWC = D4_Eiche_drip_ns_vSWC,
    D114_Eiche_drip_ns_vSWC = D6_Eiche_drip_ns_vSWC
  )

vwc_data <- vwc_orig %>%
  pivot_longer(cols = 2: ncol(vwc_orig), values_to = "SWC", 
               names_to = "Sensor.ID") %>%
  mutate(
    Treatment = case_when(
      str_detect(Sensor.ID, "drip") ~ "drip",
      str_detect(Sensor.ID, "control") ~ "control",
      TRUE ~ "unknown"
    ),
    Tree.Species = case_when(
      str_detect(Sensor.ID, "Buche") ~ "Buche",
      str_detect(Sensor.ID, "Eiche") ~ "Eiche",
      TRUE ~ "unknown"
    ),
    RainShelter = if_else(str_detect(Sensor.ID, "_rs_"), "rs", "ns"),
    Tree.ID = str_remove(Sensor.ID, "_vSWC")
  )

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
    RainShelter = if_else(str_detect(Tree.ID, "_rs"), "rs", "ns")
  )
colnames(vwc_daily)[1] <- "Date"


# calculate mean VWC per treatment and tree species from daily sensor means
vwc_daily_treat <- vwc_daily %>%
  # remove D101 and D102 since here, the swc sensor is not under the rainshelter
  filter(Tree.ID != "D101_Buche_drip_ns_vSWC" & 
           Tree.ID != "D102_Buche_drip_ns_vSWC" ) %>%
  group_by(Date, Treatment, Tree.Species, RainShelter) %>%
  mutate(SWC_mean_treat = mean(SWC_mean, na.rm = TRUE),
         SWC_min_treat = min(SWC_mean, na.rm = TRUE),
         SWC_max_treat = max(SWC_mean, na.rm = TRUE),
         SWC_sd_treat = sd(SWC_mean, na.rm = TRUE))

ggplot(vwc_daily_treat) +
  geom_ribbon(aes(x = Date, ymin = SWC_min_treat, ymax = SWC_max_treat)) +
  geom_line(aes(x = Date, y = SWC_mean_treat)) +
  facet_grid(Tree.Species ~ RainShelter ~ Treatment)

ggplot(vwc_daily_treat) +
  geom_ribbon(aes(x = Date, ymin = SWC_mean_treat - SWC_sd_treat, 
                  ymax = SWC_mean_treat + SWC_sd_treat), alpha = 0.4) +
  geom_line(aes(x = Date, y = SWC_mean_treat)) +
  ylab("Mean SWC [%]") +
  facet_grid(Tree.Species ~ RainShelter ~ Treatment)

ggsave(path =
         "./graphics/EB/SWC/",
       filename = "Mean_SWC_from_sensormeans_treat_2024.png",
       width = 2560, height = 1489, units = "px")

# create overview plot for Rainshield sensors
vwc_daily_rs <- vwc_daily %>%
  # remove D101 and D102 since here, the swc sensor is not under the rainshelter
  filter(RainShelter == "rs")

vwc_daily_rs$Treatment <- vwc_daily_rs$Treatment %>%
  fct_relevel("control", "drip") %>%  
  fct_recode(
    "Unirrigated" = "control",
    "Drip irrigated" = "drip"
  )

vwc_daily_rs$Tree.Species <- vwc_daily_rs$Tree.Species %>%
  fct_recode(
    "Fagus sylvatica" = "Buche",
    "Quercus robur" = "Eiche"
  )

# Create overview plot for rain sheltered sensors (S5)
ggplot(vwc_daily_rs %>%
         mutate(
           Tree.ID.short = paste0(
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
         )
) + 
  geom_vline(xintercept = as.Date("2024-08-01", tz = "Europe/Berlin"), 
             color = "chocolate4", size = 1) +
  geom_line(aes(x = Date, y = SWC_mean, color = Tree.ID.short)) +
  facet_grid(Tree.Species ~ Treatment) +
  xlab("2024") + ylab("Daily mean SWC [%]") +
  theme_bw()

ggsave(path =
         "./graphics/EB/SWC/",
       filename = "S5_SWC_from_sensormeans_rs_2024.png",
       width = 2560, height = 1489, units = "px")

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

# create new column with site
rel_TWD_df$Site <-  "Eilaberg"
climate_daily <- climate
climate_daily$Site <- "Eilaberg"

# combine dataframes
combined <- left_join(rel_TWD_df, climate_daily[, c(1, 2, 4:8)], 
                      by = c("Date", "Site"))
combined <- left_join(combined, vwc_daily_treat[, c(1:5, 9, 12)],
                      by = c("Date", "Tree.ID"))

# add min. rel TWD per day and Tree.ID 
combined  <- combined  %>%
 group_by(Tree.ID, Date) %>%
  mutate(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) %>%
  ungroup()

# calculate irrigation thresholds
# ITlow - calculate when min TWD > 0 
MinTWD_not0 <- combined %>%
  group_by(Date = Date, Tree.ID, Treatment, Tree.Species) %>%
  dplyr::summarise(Min_TWD = min(TWD, na.rm = TRUE)) %>%
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

FirstDayAfterLastGRO_df <- lastGRO[, c(6, 26, 26)]
colnames(FirstDayAfterLastGRO_df)[3] <- "Date"

# merge results in one dataframe
results <- reduce(list(combined, MinTWD_not0[, c(1, 2, 6)],
                       FirstDayAfterLastGRO_df,
                       p12_thres_crossed[, c(2, 3, 5)]),
                  left_join, by = c("Date", "Tree.ID")) 


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

# create new percentage column, where 0 = no IT crossed,
# 1-99 ITlow is crossed, 100 = ITup is reached and 101 = ITup is crossed
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

x_limits <- c(as.POSIXct("2024-05-01 00:00:00", tz = "Europe/Berlin"), 
              as.POSIXct("2024-09-30 23:50:00", tz = "Europe/Berlin"))

for (i in trees){
  Tree.ID_Nr <- i
  print(Tree.ID_Nr)
  filtered_data <- filter(results, Tree.ID == Tree.ID_Nr)
  selected_tree_species <- filtered_data$Tree.Species[[1]]
  rel_TWD_max_thres <- ifelse(selected_tree_species == 
                                "Quercus robur", thres_E,
                                            thres_B)
  plot_irrigation_warnings(filtered_data, 
   "./graphics/EB/Irrigation_Thresholds/Irrigation_Thresholds_Singletree/2024/")
}

# select daily irrigation warnings
warnings_daily <- results[, c(6:9, 13:32)] %>%
  group_by(Date, Tree.Species, Treatment)
warnings_daily_unique <- warnings_daily %>% distinct()

climate_daily_EB <- filter(climate_daily, 
                           Site == "Eilaberg" &
                           Date >= as.Date("2024-05-01", tz = "Europe/Berlin"))

climate_daily_EB$Date <- as.POSIXct(climate_daily_EB$Date, tz = "Europe/Berlin")

## Heatmaps
# create new column with information that no threshold was reached 
# (but we would have had data to differentiate NAs from thresholds not reached)
warnings_daily_unique <- warnings_daily_unique %>%
  mutate(Date_nothres = case_when(
    is.na(MinTWDnot0_Date) & 
      is.na(FirstDayAfterLastGRO) & 
      is.na(Date_rel_minTWD_thres_crossed) ~  
      as.Date(Date, tz = "Europe/Berlin"),
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

x_axis_limits <- c(as.POSIXct("2024-04-22", tz = "Europe/Berlin"),
                   as.POSIXct("2024-10-06", tz = "Europe/Berlin"))

x_axis_settings <- scale_x_datetime(
  date_labels = "%b",
  date_breaks = "1 month", 
  limits = x_axis_limits
)

combined$Date <- as.POSIXct(combined$Date, tz = "Europe/Berlin")
climate_daily$Date <- as.POSIXct(climate_daily$Date, tz = "Europe/Berlin")

# use data only until before 2024-10-01
climate_daily_s_EB <- climate_daily_EB %>%
  filter(Date < as.Date("2024-10-01", tz = "Europe/Berlin"))

VPD_plot_1 <- ggplot(filter(climate_daily_s_EB, 
                            Date >= as.POSIXct("2024-05-01", 
                                               tz = "Europe/Berlin"))) +
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/10), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, 
                  ymin = VPD_kPA_min,
                  ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  # irrigation dates
  geom_vline(data = climate_daily_s_EB,
             aes(xintercept = as.POSIXct("2024-08-27", 
                                         tz = "Europe/Berlin")),
             color = "chartreuse4", size = 1) +
  # rainshield installation dates
  geom_vline(xintercept = as.POSIXct("2024-08-01", tz = "Europe/Berlin"),
             color = "chocolate4", size = 1) +
  geom_line(aes(x = Date, y = VPD_kPA_mean), 
            color = "grey35", linewidth = 1, linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name="VPD [kPa]",
    limits = c(0, 4),
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


 VPD_plot_2 <- ggplot(filter(climate_daily_s_EB, 
                             Date >= as.POSIXct("2024-05-01", 
                                                tz = "Europe/Berlin"))) +
  geom_bar(aes(x = Date, y = Precip_mm_10min_sum/10), stat = "identity",   
           position = "stack", fill = "dodgerblue3") +
  geom_ribbon(aes(x = Date, 
                  ymin = VPD_kPA_min,
                  ymax = VPD_kPA_max), 
              fill ="black", alpha = 0.2) +
  # add irrigation dates
   geom_vline(data = climate_daily_s_EB,
              aes(xintercept = as.POSIXct("2024-08-27")),
              color = "chartreuse4", size = 1) +
   # only pedunculate oak
   geom_vline(data = climate_daily_s_EB,
              aes(xintercept = as.POSIXct("2024-09-05", tz = "Europe/Berlin")),
              color = "chartreuse4", size = 1) +
   # rainshield installation dates
   geom_vline(xintercept = as.POSIXct("2024-08-01", tz = "Europe/Berlin"), 
              color = "chocolate4", size = 1) +
   # rainshield deinstallation date
   geom_line(aes(x = Date, y = VPD_kPA_mean), 
             color = "grey35", linewidth = 1, linetype = "dashed") +
  x_axis_settings +
  scale_y_continuous(
    name="VPD [kPa]", 
    limits = c(0, 4),
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
                                      heatmap_data_sorted[, c(1, 22, 27)], 
                                      by = "Tree.ID")

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
    "Drip irrigated" = "drip"
  )

vwc_daily_treat_sel <- vwc_daily_treat %>%
  filter(RainShelter == "ns")

vwc_daily_treat_sel$Treatment <- vwc_daily_treat_sel$Treatment %>%
  fct_relevel("control", "drip") %>%  
  fct_recode(
    "Unirrigated" = "control",
    "Drip irrigated" = "drip"
  )

vwc_daily_treat_sel$Tree.Species <- vwc_daily_treat_sel$Tree.Species %>%
  fct_recode(
    "Fagus sylvatica" = "Buche",
    "Quercus robur" = "Eiche"
  )

# use data only until before 2024-10-01
heatmap_data_sorted_s <- heatmap_data_sorted %>%
  filter(Date < as.Date("2024-10-01", tz = "Europe/Berlin"))

FirstDayAfterLastGRO_df2_s <- FirstDayAfterLastGRO_df2 %>%
  filter(Date < as.Date("2024-10-01", tz = "Europe/Berlin"))

vwc_daily_treat_sel <- filter(vwc_daily_treat_sel,
                 Date >= as.Date("2024-05-01", tz = "Europe/Berlin") & 
                 Date <= as.Date("2024-09-30", tz = "Europe/Berlin"))

heatmap_plot_B <- ggplot() +
  geom_tile(
    data = filter(heatmap_data_sorted_s, Tree.Species == "Fagus sylvatica"), 
    aes(x = Date, y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white"
  ) +
  scale_fill_identity() +
  geom_ribbon(data = filter(vwc_daily_treat_sel,
                              Tree.Species == "Fagus sylvatica"),
              aes(x = as.POSIXct(Date),
                  ymin = (SWC_mean_treat - SWC_sd_treat) * vwc_scaling_factor + 1,
                  ymax = (SWC_mean_treat + SWC_sd_treat) * vwc_scaling_factor + 1),
              alpha = 0.5) +
  geom_line(
    data = filter(vwc_daily_treat_sel,
                    Tree.Species == "Fagus sylvatica"),
    aes(
      x = as.POSIXct(Date),
      y = SWC_mean_treat * vwc_scaling_factor + 1,  
      group = 1
    ),
    color = "black",
    linewidth = 1
  ) +
  # add FirstDayAfterLastGRO 
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2_s, 
                  Tree.Species == "Fagus sylvatica"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted_s, 
                            Tree.Species == "Fagus sylvatica" &
                              Treatment != "hose")$Tree.ID_num_cat),
    labels = 
      unique(filter(heatmap_data_sorted_s, 
                    Tree.Species == "Fagus sylvatica")$Tree.ID_num_cat),
  ) +
  geom_label(
    data = filter(heatmap_data_sorted_s, Tree.Species == "Fagus sylvatica" & 
                    Treatment != "hose"),
    aes(x = as.POSIXct("2024-04-23"), y = Tree.ID_num_cat, 
        label = gsub("^(D\\d+).*", "\\1", Tree.ID)),
    color = "black",
    fill = "white",  
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"), 
    label.r = unit(0.2, "lines")  
  ) +
  geom_label(
    data = filter(heatmap_data_sorted_s, Tree.Species == "Fagus sylvatica" &
                    Treatment != "hose" & grepl("_rs$", Tree.ID)),
    aes(x = as.POSIXct("2024-10-05"), 
        y = Tree.ID_num_cat, 
        label = gsub("^.*_(rs)$", "\\1", Tree.ID) 
    ),
    color = "black",
    fill = "white", 
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),  
    label.r = unit(0.2, "lines") 
  ) +
  facet_grid(Treatment ~ Tree.Species, scales = "free", space = "free",
             switch = 'y') +
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
    data = filter(heatmap_data_sorted_s, Tree.Species == "Quercus robur"), 
    aes(x = Date, y = Tree.ID_num_cat, fill = Percentage_color), 
    color = "white"
  ) +
  scale_fill_identity() +
  geom_ribbon(data = filter(vwc_daily_treat_sel,
                            Tree.Species == "Quercus robur"),
              aes(x = as.POSIXct(Date),
              ymin =(SWC_mean_treat - SWC_sd_treat) * vwc_scaling_factor + 1,
              ymax = (SWC_mean_treat + SWC_sd_treat) * vwc_scaling_factor + 1),
              alpha = 0.5) +
  geom_line(
    data = filter(vwc_daily_treat_sel, 
                  Tree.Species == "Quercus robur"),
    aes(
      x = as.POSIXct(Date),
      y = SWC_mean_treat * vwc_scaling_factor + 1,  
      group = 1
    ),
    color = "black",
    linewidth = 1
  ) +
  # FirstDayAfterLastGRO hinzufügen
  geom_point(
    data = filter(FirstDayAfterLastGRO_df2_s, Tree.Species == "Quercus robur"),
    aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
        y = Tree.ID_num_cat),
    shape = 19,
    color = "#00008b",
    size = 2
  ) +
  scale_y_continuous(
    breaks =  unique(filter(heatmap_data_sorted_s, 
                            Tree.Species == "Quercus robur")$Tree.ID_num_cat),
    labels = 
      unique(filter(heatmap_data_sorted_s, 
                    Tree.Species == "Quercus robur")$Tree.ID_num_cat),
    sec.axis = sec_axis(
      ~ (. - 1) / vwc_scaling_factor,  
      name = "SWC [%]"
     )
  ) +
  geom_label(
    data = filter(heatmap_data_sorted_s, Tree.Species == "Quercus robur"),
    aes(x = as.POSIXct("2024-04-23"), y = Tree.ID_num_cat, 
        label = gsub("^(D\\d+).*", "\\1", Tree.ID)),
    color = "black",
    fill = "white",  
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),  
    label.r = unit(0.2, "lines")  
  ) +
  geom_label(
    data = filter(heatmap_data_sorted_s, Tree.Species == "Quercus robur" &
                    Treatment != "hose" & grepl("_rs$", Tree.ID)),
    aes(x = as.POSIXct("2024-10-05"), 
        y = Tree.ID_num_cat, 
        label = gsub("^.*_(rs)$", "\\1", Tree.ID) 
        ),
    color = "black",
    fill = "white",  
    size = 4,
    fontface = "plain",
    label.padding = unit(0.2, "lines"),  
    label.r = unit(0.2, "lines")  
  ) +
  facet_grid(Treatment ~ Tree.Species, scales = "free", space = "free",
             switch = 'y') +
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

combined_plot_all <- heatmap_plot_B / heatmap_plot_E / VPD_plot_1 / VPD_plot_2 +
  plot_layout(ncol = 2, nrow = 2, heights = c(1, 0.15), widths = c(1, 1))

combined_plot_all

write.csv(heatmap_data_sorted_s, 
"./data/EB/results/Irrigation_thresholds/Irrigation_thres_EB_2024.csv")
write.csv(results, 
"./data/EB/results/Irrigation_thresholds/Irrigation_results_EB_2024.csv")

# Explore irrigation effect on soilmoisture, TWD and irrigation demand
climate_highres <- read_csv(
  "./data/EB/input_data/Meteo_2024.csv",
  col_select = -1)

# reformat date to CEST
climate_highres$TIME <- sub(" UTC$", " CEST", climate_highres$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
climate_highres$TIME <- ifelse(nchar(climate_highres$TIME) == 10,
                               paste(climate_highres$TIME, "00:00:00"),
                               climate_highres$TIME)
climate_highres$TIME <- as.POSIXct(climate_highres$TIME,
                                   format="%Y-%m-%d %H:%M:%S",
                                   tz = "Europe/Berlin")
climate_highres <- filter(climate_highres,
                          TIME >= as.POSIXct("2023-05-15 00:00:00",
                                             format="%Y-%m-%d %H:%M:%S",
                                             tz = "Europe/Berlin") &
                            TIME <= as.POSIXct("2023-10-01 01:00:00",
                                               format="%Y-%m-%d %H:%M:%S",
                                               tz = "Europe/Berlin"))

climate_highres_h <- climate_highres %>%
  mutate(hour = floor_date(TIME, unit = "hour")) %>%  
  group_by(hour) %>%
  summarise(VPD_kPA_mean = mean(VPD_kPA, na.rm = TRUE))

vwc_data$Tree.Species <- vwc_data$Tree.Species %>%
  fct_recode(
    "Fagus sylvatica" = "Buche",
    "Quercus robur" = "Eiche"
  )

vwc_data_means <- vwc_data %>%
  filter(RainShelter == "ns") %>%
  group_by(Treatment, Tree.Species,  hour = floor_date(TIME, unit = "hour")) %>%
  summarise(SWC_mean_treat = mean(SWC, na.rm = TRUE),
            SWC_median_treat = median(SWC, na.rm = TRUE),
            SWC_sd_treat = sd(SWC, na.rm = TRUE))

# check when SWC is back to pre irrigation level
vwc_data_means_10min <- vwc_data %>%
  filter(RainShelter == "ns") %>%
  group_by(Treatment, Tree.Species,  TIME) %>%
  summarise(SWC_mean_treat = mean(SWC, na.rm = TRUE),
            SWC_median_treat = median(SWC, na.rm = TRUE),
            SWC_sd_treat = sd(SWC, na.rm = TRUE))

vwc_data_means_check <- vwc_data_means_10min %>%
  filter(TIME >= as.POSIXct("2024-08-27 00:00:00", tz = "Europe/Berlin") & 
           TIME <= as.POSIXct("2024-08-27 23:50:00", tz = "Europe/Berlin") &
           Treatment == "drip") 

min_drip_2708_E <- vwc_data_means_10min %>%
  filter(Treatment == "drip" & Tree.Species == "Quercus robur" &
         TIME == as.POSIXct("2024-08-27 15:10:00", tz = "Europe/Berlin"))

min_drip_2708_B <- vwc_data_means_10min %>%
  filter(Treatment == "drip" & Tree.Species == "Fagus sylvatica" &
           TIME == as.POSIXct("2024-08-27 17:00:00", tz = "Europe/Berlin")) 

recovery_timeE <- vwc_data_means_10min %>%
  filter(Treatment == "drip", Tree.Species == "Quercus robur", 
         TIME > as.POSIXct("2024-08-27 15:10:00", tz = "Europe/Berlin")) %>%
  filter(SWC_mean_treat <= min_drip_2708_E$SWC_mean_treat) %>%
  slice_head(n = 1)

recovery_timeB <- vwc_data_means_10min %>%
  filter(Treatment == "drip", Tree.Species == "Fagus sylvatica", 
         TIME > as.POSIXct("2024-08-27 17:00:00", tz = "Europe/Berlin")) %>%
  filter(SWC_mean_treat <= min_drip_2708_B$SWC_mean_treat) %>%
  slice_head(n = 1)

recovery_timeE_duration <- min_drip_2708_E$TIME - recovery_timeE$TIME
recovery_timeB_duration <- as.numeric(min_drip_2708_B$TIME - recovery_timeB$TIME, 
                                      units = "hours")

max_swc_time <- vwc_data_means_check %>%
  mutate(date = as.Date(TIME, tz = "Europe/Berlin")) %>%  
  group_by(date, Tree.Species) %>%                
  slice_max(SWC_mean_treat, with_ties = FALSE) %>%  
  select(date, TIME, SWC_mean_treat)

mean_rel_TWD_log <- results %>%
  filter(RainShelter == "ns") %>%
  group_by(hour = floor_date(TIME, unit = "hour"), Treatment, Tree.Species) %>%
  summarise(rel_TWD_log_mean = mean(rel_TWD_log, na.rm = TRUE),
            rel_TWD_log_median = median(rel_TWD_log, na.rm = TRUE))

mean_percentage_w_cros <- results %>%
  filter(RainShelter == "ns") %>%
  group_by(hour = floor_date(TIME, unit = "hour"), Treatment, Tree.Species) %>%
  summarise(Percentage_w_cros_mean = mean(Percentage_w_cros, na.rm = TRUE),
            Percentage_w_cros_median = median(Percentage_w_cros, na.rm = TRUE))

mean_percentage_w_cros_2708 <- mean_percentage_w_cros %>%
  filter(hour>=as.POSIXct("2024-08-27 00:00:00", tz = "Europe/Berlin") & 
           hour<as.POSIXct("2024-08-28 00:00:00", tz = "Europe/Berlin")) %>%
  group_by(Treatment, Tree.Species)%>%
  summarise(daily_mean_perc_w_cros=mean(Percentage_w_cros_mean))

mean_percentage_w_cros_2608 <- mean_percentage_w_cros %>%
  filter(hour>=as.POSIXct("2024-08-26 00:00:00", tz = "Europe/Berlin") & 
           hour<as.POSIXct("2024-08-27 00:00:00", tz = "Europe/Berlin")) %>%
  group_by(Treatment, Tree.Species)%>%
  summarise(daily_mean_perc_w_cros=mean(Percentage_w_cros_mean))

# Irrigation 27.08. (Beech & Oak) & 05.09. (Oak)
Irrig_plot_20240827_relTWD <- ggplot() +
  geom_vline(data = filter(results, Date >= as.Date("2024-08-23") & 
                             Date <= as.Date("2024-09-03")),
             aes(xintercept = as.POSIXct(ifelse(
               Tree.Species == "Quercus robur", 
               "2024-08-27 14:15:00", "2024-08-27 16:10:00"))), 
             color = "chartreuse4", linewidth = 1) +
  
  geom_line(data = filter(mean_rel_TWD_log, 
                          hour >= as.POSIXct("2024-08-23 00:00:00") &
                            hour <= as.POSIXct("2024-09-03 00:00:00")),
            aes(x = as.POSIXct(hour), y = rel_TWD_log_mean,
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("2024") +
  ylab("Rel. TWD [µm]") +
  facet_grid(. ~ Tree.Species) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

Irrig_plot_20240827_Irrig_perc <- ggplot() +
  geom_vline(data = filter(results, Date >= as.Date("2024-08-23") & 
                             Date <= as.Date("2024-09-03")),
             aes(xintercept = as.POSIXct(ifelse(
               Tree.Species == "Quercus robur", 
               "2024-08-27 14:15:00", "2024-08-27 16:10:00"))), 
             color = "chartreuse4", linewidth = 1) +
    geom_line(data = filter(mean_percentage_w_cros, 
                          hour >= as.POSIXct("2024-08-23 00:00:00") &
                            hour <= as.POSIXct("2024-09-03 00:00:00")),
            aes(x = as.POSIXct(hour), y = Percentage_w_cros_mean,
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  ylab("Irrigation demand [%]") +
  xlab("2024") +
  facet_grid(. ~ Tree.Species) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")


Env_plot_20240827 <- ggplot() +
  geom_vline(data = filter(results, Date >= as.Date("2024-08-23") & 
                             Date <= as.Date("2024-09-03")),
             aes(xintercept = as.POSIXct(ifelse(
               Tree.Species == "Quercus robur", 
               "2024-08-27 14:15:00", "2024-08-27 16:10:00"))), 
             color = "chartreuse4", linewidth = 1) +
    geom_line(data = filter(vwc_data_means, 
                          hour >= as.POSIXct("2024-08-23 00:00:00") & 
                            hour <= as.POSIXct("2024-09-03 00:00:00")),
            aes(x = as.POSIXct(hour), 
                y = SWC_mean_treat , 
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  facet_grid(. ~ Tree.Species) +
  ylab("SWC [%]") + xlab("2024") +
  theme_bw() +
  theme(strip.text = element_blank())

Irrig_plot_20240827_combined <-  Irrig_plot_20240827_Irrig_perc /
  Irrig_plot_20240827_relTWD / 
  Env_plot_20240827  +
  plot_layout(ncol = 1, nrow = 3, heights = c(1, 1), widths = c(1, 1))

Irrig_plot_20240827_combined

ggsave(path = "./graphics/EB/Irrigation_Thresholds/Irrigation_Overview/",
       filename = "Irrig_plot_20240827_mean.png",
       width = 2560 / 300,     
       height = 1489 / 300,   
       units = "in",
       dpi = 300)

Irrig_plot_20240905_relTWD <- ggplot() +
  geom_vline(data = filter(results, Date >= as.Date("2024-09-03") & 
                             Date <= as.Date("2024-09-08") &
                             Tree.Species == "Quercus robur"),
             aes(xintercept = as.POSIXct("2024-09-05 13:20:00")), 
             color = "chartreuse4", linewidth = 1) +
  geom_line(data = filter(mean_rel_TWD_log, 
                          Tree.Species == "Quercus robur" &
                          hour >= as.POSIXct("2024-09-03 00:00:00") &
                            hour <= as.POSIXct("2024-09-08 00:00:00")),
            aes(x = as.POSIXct(hour), y = rel_TWD_log_mean,
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("2024") +
  ylab("Rel. TWD [µm]") +
  facet_grid(. ~ Tree.Species) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

Irrig_plot_20240905_Irrig_perc <- ggplot() +
  geom_vline(data = filter(results, Tree.Species == "Quercus robur" &
                           Date >= as.Date("2024-09-03") & 
                             Date <= as.Date("2024-09-08")),
             aes(xintercept = as.POSIXct("2024-09-05 13:20:00")), 
             color = "chartreuse4", linewidth = 1) +
  geom_line(data = filter(mean_percentage_w_cros,
                          Tree.Species == "Quercus robur" &
                          hour >= as.POSIXct("2024-09-03 00:00:00") &
                            hour <= as.POSIXct("2024-09-08 00:00:00")),
            aes(x = as.POSIXct(hour), y = Percentage_w_cros_mean,
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  ylab("Irrigation demand [%]") +
  xlab("2024") +
  facet_grid(. ~ Tree.Species) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")


Env_plot_20240905 <- ggplot() +
  geom_vline(data = filter(results, Tree.Species == "Quercus robur" &
                           Date >= as.Date("2024-09-03") & 
                             Date <= as.Date("2024-09-08")),
             aes(xintercept = as.POSIXct("2024-09-05 13:20:00")), 
             color = "chartreuse4", linewidth = 1) +
  geom_line(data = filter(vwc_data_means, 
                          Tree.Species == "Quercus robur" &
                          hour >= as.POSIXct("2024-09-03 00:00:00") & 
                            hour <= as.POSIXct("2024-09-08 00:00:00")),
            aes(x = as.POSIXct(hour), 
                y = SWC_mean_treat , 
                color = Treatment)) +
  scale_color_manual(values = c("red", "blue")) +
  facet_grid(. ~ Tree.Species) +
  ylab("SWC [%]") + xlab("2024") +
  theme_bw() +
  theme(strip.text = element_blank())

Irrig_plot_20240905_combined <- Irrig_plot_20240905_Irrig_perc /
  Irrig_plot_20240905_relTWD / Env_plot_20240905  +
  plot_layout(ncol = 1, nrow = 3, heights = c(1, 1), widths = c(1, 1))

Irrig_plot_20240905_combined

ggsave(path = "./graphics/EB/Irrigation_Thresholds/Irrigation_Overview/",
       filename = "Irrig_plot_20240905_mean.png",
       width = 2560 / 300,     
       height = 1489 / 300,   
       units = "in",
       dpi = 300)

