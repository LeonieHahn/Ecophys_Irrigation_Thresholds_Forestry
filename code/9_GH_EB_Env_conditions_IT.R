# =============================================================================
# - Compare environmental conditions when the irrigation thresholds 
#   where reached for GH, EB in both years
# =============================================================================
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(zoo)
library(readxl)

source("./code/Functions.R")

irrig_thres_GH <- read_csv(
  "./data/GH/results/Irrigation_thresholds/Irrigation_thres_GH_2023.csv",
  col_select = -1)

colnames(irrig_thres_GH)[8] <- "SWC"

# use data only before last date of irrigation 
irrig_thres_GH <- filter(irrig_thres_GH, Date < as.Date("2023-09-07")) 

stats_Watered_perc_GH <- irrig_thres_GH %>%
  filter(Treatment == "Watered") %>%
  group_by(Tree.Species) %>%
  summarise(min_Perc_w_cros = min(Percentage_w_cros, na.rm = TRUE),
            max_Perc_w_cros = max(Percentage_w_cros, na.rm = TRUE),
            mean_Perc_w_cros = mean(Percentage_w_cros, na.rm = TRUE))

irrig_thres_EB_2023 <- read_csv(
  "./data/EB/results/Irrigation_thresholds/Irrigation_thres_EB_2023.csv", 
  col_select = -1)

irrig_thres_EB_2024 <- read_csv(
  "./data/EB/results/Irrigation_thresholds/Irrigation_thres_EB_2024.csv", 
  col_select = -1)

irrig_thres_GH <- irrig_thres_GH %>%
  mutate(criterium = ifelse(Percentage_w_cros == 0, 
                            "no irrig. threshold reached",
                            ifelse(Percentage_w_cros == 101, "B crossed",
                                   ifelse(Percentage_w_cros == 100, 
                                          "B", "A"))),
         Location = "Greenhouse")%>%
  select(Location, Treatment, Tree.Species, Date, criterium, SWC,
         VPD_adj_dailyavg_kPa)

colnames(irrig_thres_GH)[6] <- "SWC_mean"
colnames(irrig_thres_GH)[7] <- "VPD_kPA_mean"

irrig_thres_EB_2023 <- irrig_thres_EB_2023 %>%
  mutate(criterium = ifelse(Percentage_w_cros == 0, 
                            "no irrig. threshold reached",
                            ifelse(Percentage_w_cros == 101, "B crossed",
                                   ifelse(Percentage_w_cros == 100, 
                                          "B", "A"))),
         Location = "EB2023")  %>%
  filter(Treatment == "Unirrigated") %>%
  select(Location, Treatment, Tree.Species, Date, criterium, SWC_mean,
         VPD_kPA_mean)

irrig_thres_EB_2024 <- irrig_thres_EB_2024 %>%
  mutate(criterium = ifelse(Percentage_w_cros == 0, 
                            "no irrig. threshold reached",
                            ifelse(Percentage_w_cros == 101, "B crossed",
                                   ifelse(Percentage_w_cros == 100, 
                                          "B", "A"))),
         Location = "EB2024") %>%
  filter(RainShelter == "ns") %>%
  select(Location, Treatment, Tree.Species, Date, criterium, SWC_mean,
         VPD_kPA_mean)

irrig <- rbind(irrig_thres_GH, irrig_thres_EB_2023,irrig_thres_EB_2024) %>%
  mutate(Treatment = case_when(Treatment == "Unirrigated" ~ "Droughted",
                               Treatment == "Drip irrigated" ~ "Watered",
                               TRUE ~ Treatment))

irrig_stats_overall <- irrig %>%
  group_by(criterium, Location, Tree.Species) %>%
  summarise(min_SWC_wI = round(min(SWC_mean, na.rm = TRUE), 2),
            mean_SWC_wI = round(mean(SWC_mean, na.rm = TRUE), 2),
            max_SWC_wI = round(max(SWC_mean, na.rm = TRUE), 2),
            sd_SWC_wI = round(sd(SWC_mean, na.rm = TRUE), 2),
            min_meanVPD_wI = round(min(VPD_kPA_mean, na.rm = TRUE), 2),
            max_meanVPD_wI = round(max(VPD_kPA_mean, na.rm = TRUE), 2),
            mean_meanVPD_wI = round(mean(VPD_kPA_mean, na.rm = TRUE), 2),
            sd_meanVPD_wI = round(sd(VPD_kPA_mean, na.rm = TRUE), 2))

irrig_stats_treatment <- irrig %>%
  group_by(criterium, Location, Tree.Species, Treatment) %>%
  summarise(min_SWC_wI = round(min(SWC_mean, na.rm = TRUE), 2),
            mean_SWC_wI = round(mean(SWC_mean, na.rm = TRUE), 2),
            max_SWC_wI = round(max(SWC_mean, na.rm = TRUE), 2),
            sd_SWC_wI = round(sd(SWC_mean, na.rm = TRUE), 2),
            min_meanVPD_wI = round(min(VPD_kPA_mean, na.rm = TRUE), 2),
            max_meanVPD_wI = round(max(VPD_kPA_mean, na.rm = TRUE), 2),
            mean_meanVPD_wI = round(mean(VPD_kPA_mean, na.rm = TRUE), 2),
            sd_meanVPD_wI = round(sd(VPD_kPA_mean, na.rm = TRUE), 2))

write.csv(irrig_stats_overall,
"./data/EB/results/Irrigation_thresholds/GH_FS_Irrig_Env_Comparison.csv")
write.csv(irrig_stats_treatment,
"./data/EB/results/Irrigation_thresholds/GH_FS_Irrig_Env_Comparison_treatment.csv")

# Overview boxplots
ggplot() +
  geom_boxplot(data = filter(irrig),
               aes(x = Location, y = VPD_kPA_mean, 
                   color = criterium))+
  scale_color_manual(values=c("no irrig. threshold reached" = "lightblue", 
                              "A" = "goldenrod2",
                              "B crossed" = "darkmagenta",
                              "B"="mediumorchid3"))+
  facet_grid(Tree.Species ~ Treatment ~ criterium) +
  theme_bw()+
  ggtitle("Mean VPD on days where irrigation criterium is fullfilled
          (dried trees only)")+
  ylab("Mean VPD [kPA]")+xlab("Irrigation criterium")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")

ggplot() +
  geom_boxplot(data = filter(irrig),
               aes(x = Location, y = SWC_mean, 
                   color = criterium))+
  scale_color_manual(values=c("no irrig. threshold reached" = "lightblue", 
                              "A" = "goldenrod2",
                              "B crossed" = "darkmagenta",
                              "B"="mediumorchid3"))+
  facet_grid(Tree.Species ~ Treatment ~ criterium) +
  theme_bw() +
  ylab("Mean SWC [%]")+xlab("")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="none")
