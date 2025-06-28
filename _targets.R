library(targets)

# Load functions
tar_source("code/Functions.R")

# Define targets
list(
  tar_target(step1, {
    source("code/1_GH_Compare_TWD.R")
    TRUE
  }),
  tar_target(step2, {
    source("code/2_GH_Calculate_ITup.R")
    TRUE
  }),
  tar_target(step3, {
    source("code/3_GH_Apply_ITs.R")
    TRUE
  }),
  tar_target(step4, {
    source("code/4_GH_IT_Env_Statistic_Analysis.R")
    TRUE
  }),
  tar_target(step5, {
    source("code/5_EB_Apply_ITs_2023.R")
    TRUE
  }),
  tar_target(step6, {
    source("code/6_EB_Apply_ITs_2024.R")
    TRUE
  }),
  tar_target(step7, {
    source("code/7_EB_IT_Env_Statistic_Analysis_2023.R")
    TRUE
  }),
  tar_target(step8, {
    source("code/8_EB_IT_Env_Statistic_Analysis_2024.R")
    TRUE
  }),
  tar_target(step9, {
    source("code/9_GH_EB_Env_conditions_IT.R")
    TRUE
  })
)
