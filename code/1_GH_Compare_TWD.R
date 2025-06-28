# =============================================================================
# - Calculate mean, min. and max. per timestep, treatment and treatment of 
#   dendrometer values
# - Create overview plot (figure 3)
# =============================================================================

library(readr)
library(dplyr)
library(ggplot2)

# table with output from DendRoAnalyst (phase.zg) 
ZG_phase_DF <- read_csv("./data/GH/input_data/ZG_phase_DF.csv", 
                        col_select = -1)

# change tree species and treatment labels
ZG_phase_DF$Tree.Species <- factor(ZG_phase_DF$Tree.Species,
                                   levels = c("B", "E", "D"),
                                   labels = c("B" = "Fagus sylvatica",
                                              "E" = "Quercus robur",
                                              "D" = "Pseudotsuga menziesii"))

ZG_phase_DF$Treatment <- factor(ZG_phase_DF$Treatment, 
                                labels = c("D" = "Droughted", 
                                           "W" = "Watered"))

# calculate min, max and mean of dm, TWD and GRO per day, treatment and 
# tree  species
ZG_phase_DF_stats <- ZG_phase_DF %>%
  group_by(TIME, Treatment, Tree.Species) %>%
  summarise(Min_TSD = min(dm),
            Max_TSD = max(dm),
            Mean_TSD = mean(dm),
            Min_TWD = min(TWD),
            Max_TWD = max(TWD),
            Mean_TWD = mean(TWD))

# plot daily mean TWD and its ranges
ggplot(ZG_phase_DF_stats) +
  # add last day of irrigation
  geom_vline(aes(xintercept = as.POSIXct("2023-06-26", tz = "Europe/Berlin")),
             color = "chartreuse4", linewidth = 1) +
  geom_ribbon(aes(x = TIME, ymin = Min_TWD, ymax = Max_TWD, 
                  fill = Treatment), alpha = 0.3)+
  geom_line(aes(x = TIME,
                y = Mean_TWD, 
                color = Treatment)) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  facet_grid(Tree.Species ~ .,
             labeller = labeller(Tree.Species = as_labeller(c(
               "Fagus sylvatica" = "F. sylvatica",
               "Pseudotsuga menziesii" = "P. menziesii",
               "Quercus robur" = "Q. robur"
             )))) + 
  theme_bw() +
  xlab("2023") + ylab("TWD [µm]") +
  scale_x_datetime(date_labels = "%d %b") +
  theme_bw(base_size = 14) + 
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
ggsave(
  path = "./graphics/GH/Dendro_Metrics_DA/TWD/DA-TWD_Overview/",
  filename = "Fig.3_Overview_MeanMinMax_TWD.JPEG",
  width = 7.48,
  height = 5.5,
  dpi = 300,
  units = "in")

