# Function for calculating confidence intervals with bootstrapping for logistic 
# regression model to predict TWD from WP
calc_CI_boot <- function(data, predictor, response, n_bootstrap = 1000, 
                         pred_range = seq(-9, 0, by = 0.01)) {
  if (!(predictor %in% names(data)) || !(response %in% names(data))) {
    stop("Predictor- or targetvariable colum don't exist in the data")
  }
  
  # Remove incomplete cases
  filtered_data <- data[complete.cases(data[[predictor]], data[[response]]), ]

  bootstrap_params <- data.frame(Asym = numeric(n_bootstrap),
                                 xmid = numeric(n_bootstrap),
                                 scal = numeric(n_bootstrap))
  pred_values <- data.frame(predictor = pred_range)
  colnames(pred_values) <- predictor
  bootstrap_preds <- matrix(NA, nrow = n_bootstrap, ncol = nrow(pred_values))
  
  # Apply bootstrapping 
  for (i in 1:n_bootstrap) {
    # take sample
    boot_data <- filtered_data[sample(1:nrow(filtered_data), replace = TRUE), ]
    
    model_boot <- tryCatch(
      nls(as.formula(paste(
        response, "~ SSlogis(", predictor, ", Asym, xmid, scal)")), 
        data = boot_data),
      error = function(e) NA
    )
    
    # save model, if it converges
    if (inherits(model_boot, "nls")) {
      bootstrap_params[i, ] <- coef(model_boot)
      bootstrap_preds[i, ] <- predict(model_boot, newdata = pred_values)
    }
  }
  
  # Remove missing values from the bootstrapping results
  bootstrap_params <- na.omit(bootstrap_params)
  bootstrap_preds <- bootstrap_preds[complete.cases(bootstrap_preds), ]
  
  #  Calculate confidence intervals for the parameters (90% CI)
  param_CI <- bootstrap_params %>%
    reframe(across(everything(), ~ quantile(.x, probs = c(0.05, 0.90)), 
                   .names = "CI_{.col}"))
  
  # Calculate confidence intervals for the predictions (90% CI)
  pred_values$lower_CI <- apply(bootstrap_preds, 2, quantile, probs = 0.05)
  pred_values$upper_CI <- apply(bootstrap_preds, 2, quantile, probs = 0.90)
  
  pred_values$mean_pred <- colMeans(bootstrap_preds)
  
  list(
    param_CI = param_CI,
    pred_values = pred_values
  )
}

# Function for creating a logistic model to predict rel. twd from wp,
# compare the values predicted with the model to the observed data,
# extract rel. twd values for p12 according to tree species and
# plot model predictions as well as derived rel. TWD values for p12 
log_pred_twd_wp <- function(data, tree_species){
  # Filter according to tree species
  species_data <- subset(data, Tree.Species == tree_species)
  
  # model for the prediction of TWD by WP
  model_MinTWD_PDWP <- tryCatch(
    nls(Min_rel_TWD ~ SSlogis(WP_predawn, Asym, xmid, scal), 
        data = species_data),
    error = function(e) NA
  )
  model_MaxTWD_PDWP <- tryCatch(
    nls(Max_rel_TWD ~ SSlogis(WP_predawn, Asym, xmid, scal), 
        data = species_data),
    error = function(e) NA
  )
  model_MinTWD_MDWP <- tryCatch(
    nls(Min_rel_TWD ~ SSlogis(WP_midday, Asym, xmid, scal), 
        data = species_data),
    error = function(e) NA
  )
  model_MaxTWD_MDWP <- tryCatch(
    nls(Max_rel_TWD ~ SSlogis(WP_midday, Asym, xmid, scal), 
        data = species_data),
    error = function(e) NA
  )
  
  model_summaries <- list(
    MinTWD_PDWP = if (inherits(model_MinTWD_PDWP, "nls")) 
      summary(model_MinTWD_PDWP) else NA,
    MaxTWD_PDWP = if (inherits(model_MaxTWD_PDWP, "nls")) 
      summary(model_MaxTWD_PDWP) else NA,
    MinTWD_MDWP = if (inherits(model_MinTWD_MDWP, "nls")) 
      summary(model_MinTWD_MDWP) else NA,
    MaxTWD_MDWP = if (inherits(model_MaxTWD_MDWP, "nls")) 
      summary(model_MaxTWD_MDWP) else NA
  )
  
  model_info <- list(model_info, model_summaries)
  
  # add predicted TWD values to subset of dataframe
  species_data[["logistic_pred_MinTWD_PDWP"]] <- 
    if (inherits(model_MinTWD_PDWP, "nls")) 
      predict(model_MinTWD_PDWP, newdata = species_data) else NA
  species_data[["logistic_pred_MaxTWD_PDWP"]] <- 
    if (inherits(model_MaxTWD_PDWP, "nls")) 
      predict(model_MaxTWD_PDWP, newdata = species_data) else NA
  species_data[["logistic_pred_MinTWD_MDWP"]] <- 
    if (inherits(model_MinTWD_MDWP, "nls")) 
      predict(model_MinTWD_MDWP, newdata = species_data) else NA
  species_data[["logistic_pred_MaxTWD_MDWP"]] <- 
    if (inherits(model_MaxTWD_MDWP, "nls")) 
      predict(model_MaxTWD_MDWP, newdata = species_data) else NA
  
  model_results <- rbind(model_results, species_data)
  
  # compare logistically predicted TWD with measured observations with a linear
  # regression (see Dietrich et al., 2018)
  model_log_pred_MinTWD_PDWP <- tryCatch(
    lm(logistic_pred_MinTWD_PDWP ~ Min_rel_TWD, data = species_data),
    error = function(e) NA
  )
  model_log_pred_MaxTWD_PDWP <- tryCatch(
    lm(logistic_pred_MaxTWD_PDWP ~ Max_rel_TWD, data = species_data),
    error = function(e) NA
  )
  model_log_pred_MinTWD_MDWP <- tryCatch(
    lm(logistic_pred_MinTWD_MDWP ~ Min_rel_TWD, data = species_data),
    error = function(e) NA
  )
  model_log_pred_MaxTWD_MDWP <- tryCatch(
    lm(logistic_pred_MaxTWD_MDWP ~ Max_rel_TWD, data = species_data),
    error = function(e) NA
  )
  
  model_comp_observ_summaries <- list(
    MinTWD_PDWP = if (inherits(model_log_pred_MinTWD_PDWP, "lm")) summary(
      model_log_pred_MinTWD_PDWP) else NA,
    MaxTWD_PDWP = if (inherits(model_log_pred_MaxTWD_PDWP, "lm")) summary(
      model_log_pred_MaxTWD_PDWP) else NA,
    MinTWD_MDWP = if (inherits(model_log_pred_MinTWD_MDWP, "lm")) summary(
      model_log_pred_MinTWD_MDWP) else NA,
    MaxTWD_MDWP = if (inherits(model_log_pred_MaxTWD_MDWP, "lm")) summary(
      model_log_pred_MaxTWD_MDWP) else NA
  )
  model_info_comp_observ <- c(model_info_comp_observ, 
                              model_comp_observ_summaries)
  
  MinTWD_PDWP_r2 = if (inherits(model_log_pred_MinTWD_PDWP, "lm")) summary(
    model_log_pred_MinTWD_PDWP)$r.squared else NA
  MaxTWD_PDWP_r2 = if (inherits(model_log_pred_MaxTWD_PDWP, "lm")) summary(
    model_log_pred_MaxTWD_PDWP)$r.squared else NA
  MinTWD_MDWP_r2 = if (inherits(model_log_pred_MinTWD_MDWP, "lm")) summary(
    model_log_pred_MinTWD_MDWP)$r.squared else NA
  MaxTWD_MDWP_r2 = if (inherits(model_log_pred_MaxTWD_MDWP, "lm")) summary(
    model_log_pred_MaxTWD_MDWP)$r.squared else NA
  tree_species
  
  R2_df <- rbind(R2_df, 
                 setNames(data.frame(MinTWD_PDWP_r2, MaxTWD_PDWP_r2,
                                     MinTWD_MDWP_r2, MaxTWD_MDWP_r2, 
                                     tree_species),
                          c("MinTWD_PDWP_r2", "MaxTWD_PDWP_r2",
                            "MinTWD_MDWP_r2", "MaxTWD_MDWP_r2",
                            "Tree.Species")))
  
  # Predict TWD for a set of created WP values with the created model 
  # (needed for plotting)
  wp_values <- setNames(data.frame(seq(0, -9, by = -0.01)), "WP_created")
  wp_values$logistic_pred_MinTWD_PDWP <- if (
    inherits(model_MinTWD_PDWP, "nls")) {
    predict(model_MinTWD_PDWP, 
            newdata = data.frame(WP_predawn = seq(0, -9, by = -0.01)))
  } else {
    rep(NA, length(seq(0, -9, by = -0.01)))  
  }
  wp_values$logistic_pred_MaxTWD_PDWP <- if (
    inherits(model_MaxTWD_PDWP, "nls")) {
    predict(model_MaxTWD_PDWP, 
            newdata = data.frame(WP_predawn = seq(0, -9, by = -0.01)))
  } else {
    rep(NA, length(seq(0, -9, by = -0.01)))  
  }
  wp_values$logistic_pred_MinTWD_MDWP <- if (
    inherits(model_MinTWD_MDWP, "nls")) {
    predict(model_MinTWD_MDWP, 
            newdata = data.frame(WP_midday = seq(0, -9, by = -0.01)))
  } else {
    rep(NA, length(seq(0, -9, by = -0.01)))  
  }
  wp_values$logistic_pred_MaxTWD_MDWP <- if (
    inherits(model_MaxTWD_MDWP, "nls")) {
    predict(model_MaxTWD_MDWP, 
            newdata = data.frame(WP_midday = seq(0, -9, by = -0.01)))
  } else {
    rep(NA, length(seq(0, -9, by = -0.01)))  
  }
  
  wp_created_df <- rbind(wp_created_df, wp_values)
  wp_created_df$Tree.Species <- tree_species
  
  # Logistically predict TWD corresponding to p12 (threshold_p12)
  p12_df <- setNames(data.frame(
    get(paste0("threshold_p12_", tree_species)),
    get(paste0("threshold_p12_", tree_species))), 
    c("WP_predawn", "WP_midday"))
  threshold_p12_MinTWD_PDWP <- if (inherits(
    model_MinTWD_PDWP, "nls")) predict(model_MinTWD_PDWP, 
                                       newdata = p12_df) else NA
  threshold_p12_MaxTWD_PDWP <- if (inherits(
    model_MaxTWD_PDWP, "nls")) predict(model_MaxTWD_PDWP, 
                                       newdata = p12_df) else NA
  threshold_p12_MinTWD_MDWP <- if (inherits(
    model_MinTWD_MDWP, "nls")) predict(model_MinTWD_MDWP, 
                                       newdata = p12_df) else NA
  threshold_p12_MaxTWD_MDWP <- if (inherits(
    model_MaxTWD_MDWP, "nls")) predict(model_MaxTWD_MDWP, 
                                       newdata = p12_df) else NA
  
  p12_thres <- rbind(p12_thres, 
                     setNames(data.frame(
                       threshold_p12_MinTWD_PDWP, 
                       threshold_p12_MaxTWD_PDWP,
                       threshold_p12_MinTWD_MDWP,
                       threshold_p12_MaxTWD_MDWP,
                       tree_species),
                       c("threshold_p12_MinTWD_PDWP",
                         "threshold_p12_MaxTWD_PDWP",
                         "threshold_p12_MinTWD_MDWP",
                         "threshold_p12_MaxTWD_MDWP",
                         "Tree.Species")))
  
  # calculate confidence intervals with bootstrapping
  result_CI_MinTWD_PDWP <- calc_CI_boot(
    data = species_data,
    predictor = "WP_predawn",
    response = "Min_rel_TWD"
  )
  
  result_CI_MaxTWD_PDWP <- calc_CI_boot(
    data = species_data,
    predictor = "WP_predawn",
    response = "Max_rel_TWD"
  )
  
  result_CI_MinTWD_MDWP <- calc_CI_boot(
    data = species_data,
    predictor = "WP_midday",
    response = "Min_rel_TWD"
  )
  
  result_CI_MaxTWD_MDWP <- calc_CI_boot(
    data = species_data,
    predictor = "WP_midday",
    response = "Max_rel_TWD"
  )
  
  CI_summary_MinTWD_PDWP <- data.frame(
    Model = "MinTWD_PDWP",
    WP = result_CI_MinTWD_PDWP[[2]]$WP_predawn,
    Lower_CI = result_CI_MinTWD_PDWP[[2]]$lower_CI,
    Upper_CI =  result_CI_MinTWD_PDWP[[2]]$upper_CI,
    Tree.Species = tree_species 
  )
  CI_summary_MaxTWD_PDWP <- data.frame(
    Model = "MaxTWD_PDWP",
    WP = result_CI_MaxTWD_PDWP[[2]]$WP_predawn,
    Lower_CI = result_CI_MaxTWD_PDWP[[2]]$lower_CI,
    Upper_CI =  result_CI_MaxTWD_PDWP[[2]]$upper_CI,
    Tree.Species = tree_species 
  )
  CI_summary_MinTWD_MDWP <- data.frame(
    Model = "MinTWD_MDWP",
    WP = result_CI_MinTWD_MDWP[[2]]$WP_midday,
    Lower_CI = result_CI_MinTWD_MDWP[[2]]$lower_CI,
    Upper_CI =  result_CI_MinTWD_MDWP[[2]]$upper_CI,
    Tree.Species = tree_species 
  )
  CI_summary_MaxTWD_MDWP <- data.frame(
    Model = "MaxTWD_MDWP",
    WP = result_CI_MaxTWD_MDWP[[2]]$WP_midday,
    Lower_CI = result_CI_MaxTWD_MDWP[[2]]$lower_CI,
    Upper_CI =  result_CI_MaxTWD_MDWP[[2]]$upper_CI,
    Tree.Species = tree_species 
  )
  
  CI_summary <- reduce(list(CI_summary_MinTWD_PDWP,
                            CI_summary_MaxTWD_PDWP,
                            CI_summary_MinTWD_MDWP,
                            CI_summary_MaxTWD_MDWP), rbind)
  
  CI_df <- rbind(CI_df, CI_summary)
  
  return(list(model_info = model_info, 
              model_results = model_results, 
              model_info_comp_observ = model_info_comp_observ,
              wp_created_df = wp_created_df, 
              p12_thres = p12_thres,
              R2_df = R2_df,
              CI_df = CI_df))
} 

# function for recalculating rel. TWD percentage values from 
# 0 (= no irrigation threshold reached), 1 (= ITlow reached) 
# to 100 (= ITup reached)
calc_thres_perc <- function(df_results, df_base_perc, 
                            treespecies_l, treespecies_s){
  df_results <- results %>%
    filter(Tree.Species == treespecies_l) %>%
    mutate(
      rounded_rel_TWD = sapply(min_rel_TWD, function(x) {
        thres_values <- df_base_perc %>%
          filter(Tree.Species == treespecies_s) %>%  
          pull(p12_thres_scaled)
        
        # Find closest value
        rounded_value <- thres_values[which.min(abs(thres_values - x))]
        round(rounded_value, 8) 
      })
    )
  
  # Join percentages to rounded values
  results_df <- left_join(
    df_results,
    df_base_perc %>%
      filter(Tree.Species == treespecies_s) %>%  
      select(Percentage, rounded_rel_TWD
      ), 
    by = c("rounded_rel_TWD" = "rounded_rel_TWD") 
  )
}

# Function to plot logistic regression and TWD values derived for p12
plot_log_models_p12_thres <- function(model_data_df, CI_data, tree_species,
                                      p12_data, w, i,
                                      logreg_path){
  
  # Filter according to tree species and model type
  model_data_df <- filter(model_data_df, Tree.Species == tree_species)
  
  wp_df_continous <- filter(w, Tree.Species == tree_species)
  p12_thres <- filter(p12_data, Tree.Species == tree_species)
  p12_twd_value <- p12_thres[[i, 3]]
  cont_log_pred_wp_name <- colnames(wp_df_continous)[i + 1]
  response_var <- model_combinations[[i]][[1]] # WP
  predictor_var <- model_combinations[[i]][[2]] # TWD
  p12 <-  get(paste0("threshold_p12_", tree_species))
  
  if(!is.na(p12_twd_value)){
    
    TWD_type <- ifelse(model_combinations[[i]][1] == "Min_rel_TWD", "MinTWD", 
                       ifelse(model_combinations[[i]][1] == "Max_rel_TWD", 
                              "MaxTWD", NA))
    
    WP_type <- ifelse(model_combinations[[i]][2] == "WP_predawn", "PDWP", 
                      ifelse(model_combinations[[i]][2] == "WP_midday", 
                             "MDWP", NA))
    model_type <- paste(TWD_type, WP_type, sep = "_")
    CI_df <- filter(CI_data, Tree.Species == tree_species & Model == model_type)
    
    wp_df_continous <- wp_df_continous %>%
      mutate(fill_gradient = ifelse(WP_created > p12 - 0.01,
                                    scales::rescale(WP_created,
                                                    from = c(0, p12)),
                                    NA))
    # Scatterplot with TWD and WP, logistic regression line and derived TWD 
    # values for p12 -> ITup
    ggplot() +
      # Add confidence intervals
      geom_ribbon(data = CI_df,
                  aes(x = WP, ymin = Lower_CI, ymax = Upper_CI),
                  fill = "grey60", alpha = 0.3) +
      # add logistic regression
      geom_line(data = wp_df_continous,
                aes(x = WP_created,
                    y = .data[[cont_log_pred_wp_name]]),
                linewidth = 1) +
      
      
      # add measurements
      geom_point(data = model_data_df,
                 aes(x = .data[[predictor_var]],
                     y = .data[[response_var]],
                     color = Treatment, shape = Treatment), size = 3) +
      scale_color_manual(values=c(D = "red", W = "blue")) +
      
      labs(y = paste(sub("^(.*?)_.*", "\\1", response_var), ". rel. TWD"),
           x =  paste0("Ψ", sub(".*_(.*)", "\\1", predictor_var), " [MPa]")) +
      ylim(0, 1) +
      xlim(-9, 0) +
      
      # add p12-thresholds
      annotate("segment",
               x = -Inf, xend = get(paste0("threshold_p12_", tree_species)),
               y = p12_twd_value,
               yend = p12_twd_value,
               linewidth = 1, linetype = "dashed", colour = "darkmagenta") +
      annotate("segment",
               x = get(paste0("threshold_p12_", tree_species)),
               xend = get(paste0("threshold_p12_", tree_species)), y = -Inf,
               yend = p12_twd_value,
               linewidth = 1, linetype = "dashed", colour = "darkmagenta") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size=18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12)) 
    
  } else {
    # Scatterplot with TWD and WP if logistic model could not be fitted
    ggplot(model_data_df, aes(x = .data[[predictor_var]], 
                              y = .data[[response_var]])) +
      geom_point(size = 3, aes(color = Treatment, shape = Treatment)) +
      scale_color_manual(values = c("red", "blue")) +
      labs(y = paste(sub("^(.*?)_.*", "\\1", response_var), ". rel. TWD"),
           x = paste0("Ψ", sub(".*_(.*)", "\\1", predictor_var), " [MPa]")) +
      ylim(0, 1) +
      xlim(-9, 0) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size=18),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12)) 
  }
  
  ggsave(path = logreg_path,
         filename = paste0(
           "Log_reg_", response_var, predictor_var, "_", tree_species, ".png"),
         width = 2560, height = 1489, units = "px")
}


# function for plotting the irrigation thresholds with environmental 
# measurements per tree, add also original TSD for the greenhouse trees
plot_irrigation_thresholds <- function(df, plotpath){
  
  # plot with VPD and SWC
  env_plot <- ggplot()+ 
    geom_vline(aes(xintercept = as.POSIXct("2023-06-26")),
               color = "chartreuse4", size = 1) +
    geom_ribbon(data = df, aes(x = Date, 
                               ymin = VPD_adj_dailymin_kPa,
                               ymax = VPD_adj_dailymax_kPa), 
                fill ="black", alpha = 0.2) +
    geom_line(data = df, aes(x = Date, 
                             y = VPD_adj_dailyavg_kPa), 
              linetype = "dashed") +
    geom_line(data = df, 
              aes(x = Date, 
                  y = VWC_mean / 10), color = "blue") + 
    scale_y_continuous(name = "VPD [kPa]", limits = c(0, 6), 
                       breaks = seq(0, 6, 2), 
                       sec.axis = sec_axis(~ . * 10, name = "SWC [%]", 
                                           breaks = c(0, 20, 40, 60))) + 
    xlab("2023") +  
    scale_x_datetime(limits = c(as.POSIXct("2023-06-21"), 
                                as.POSIXct("2023-09-06")),
                     date_labels = "%d %b") +
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = "black"),
      axis.text.y.left = element_text(color = "black"),
      axis.title.y.right = element_text(color = "blue"),
      axis.text.y.right = element_text(color = "blue"),
      axis.ticks.y.left = element_line(color = "black"),
      axis.ticks.y.right = element_line(color = "blue")
    )
  
  # Plot with jump-free and interpolated tree stem diameter 
  TSD_plot <- ggplot(df)+
    geom_vline(aes(xintercept = as.POSIXct("2023-06-26")),
               color = "chartreuse4", size = 1) +
    geom_line(aes(x = as.POSIXct(TIME), y = dm))+
    scale_x_datetime(limits = c(as.POSIXct("2023-06-21"), 
                                as.POSIXct("2023-09-06")), 
                     date_labels = "%d %b")+
    labs(title = paste0(Tree.ID_Nr, " (", df$Tree.Species[1], 
                        " - ", df$Treatment[1], ")"),
         x = "Date", y = "TSD [µm]") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Plot with ITlow and growth stop
  TWD_plot <- ggplot(df) +
    geom_vline(aes(xintercept = as.POSIXct("2023-06-26")),
               color = "chartreuse4", size = 1) +
    # Highlight MinTWDnot0_Date (=ITlow)
    geom_rect(data = filtered_data %>% filter(!is.na(MinTWDnot0_Date)),
              aes(xmin = as.POSIXct(TIME), xmax = as.POSIXct(TIME + 1440),
                  ymin = 0, ymax = Inf),
              fill = "darkgoldenrod2") +
    # Highlight first day after growth stopped
    geom_vline(data = filtered_data %>% 
                 filter(!is.na(FirstDayAfterLastGRO)) %>% slice(1),
               aes(xintercept = as.POSIXct(TIME)),
               color = "#4a2d9f", size = 2) +
    # Add TWD
    geom_line(aes(x = as.POSIXct(TIME), y = TWD)) +
    scale_x_datetime(limits = c(as.POSIXct("2023-06-21"), 
                                as.POSIXct("2023-09-06")), 
                     date_labels = "%d %b") +
    labs(x = "Date",
         y = "TWD [µm]") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # Plot with rel. TWD and percentage of reaching ITup,
  # Irrigation threshold ITlow and ITup are highlighted as well
  perc_plot <- ggplot(df) +
    
    # add rectangle with percentage of reaching ITup
    geom_rect(
      aes(xmin = Date,
          xmax = Date + lubridate::days(1),
          ymin = -Inf, ymax = Inf, fill = Percentage_color),
      alpha = 0.4
    ) +
    scale_fill_identity() +
    
    # # add 12% Xylem conductivity loss dates (ITup)
    geom_rect(data = df %>% filter(!is.na(Date_rel_minTWD_thres_crossed)),
              aes(xmin = as.POSIXct( TIME), xmax = as.POSIXct(TIME + 1440),
                  ymin = -Inf, ymax = Inf),
              fill = "darkmagenta") +
    
    geom_vline(aes(xintercept = as.POSIXct("2023-06-26")),
               color = "chartreuse4", size = 1) +
    
    # Add 12% Xylem conductivity loss threshold (ITup) line
    geom_hline(yintercept = rel_TWD_max_thres, linetype = "dotted",
               color = "black") +
    
    # Add rel. TWD
    geom_line(aes(x = as.POSIXct(TIME), y = rel_TWD_log)) +
    scale_x_datetime(limits = c(as.POSIXct("2023-06-21"),
                                as.POSIXct("2023-09-06")), 
                     date_labels = "%d %b") +
    labs(x = "Date",
         y = "Rel. TWD") + 
    theme_bw() +
    ylim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # combine plots
  combined_plot_area <-  TSD_plot / TWD_plot /  perc_plot / env_plot + 
    plot_layout(ncol = 1, nrow = 4, heights = c(1, 1))
  
  combined_plot_area
  ggsave(path = plotpath,
         filename = paste0(Tree.ID_Nr, "_Irrigation_Thresholds.png"),
         width = 2560, height = 1489, units = "px") 
}

# Function for finding max. percentage_w_cros where each tree did still recover 
# to an irrigation demand of 0
find_threshold_change <- function(df) {
  df <- df %>% arrange(Date) 
  
  # Index for when the irrigation demand returned to 0 for the last time
  zero_fall_indices <- which(df$Percentage_w_cros == 0 & 
                               lag(df$Percentage_w_cros, default = 1) > 0)
  
  if (length(zero_fall_indices) == 0) {
    return(data.frame(Tree.ID = df$Tree.ID[1], 
                      Max_Before_Zero = NA, 
                      Last_Zero_Value = NA, 
                      First_Increasing_Value = NA))
  }
  
  zero_fall_index <- max(zero_fall_indices, na.rm = TRUE)
  
  if (zero_fall_index <= 1) {
    return(data.frame(Tree.ID = df$Tree.ID[1], 
                      Max_Before_Zero = NA, 
                      Last_Zero_Value = 0, 
                      First_Increasing_Value = NA))
  }
  
  # Max. irrigation demand before last 0 value
  max_before_zero <- max(df$Percentage_w_cros[1:(zero_fall_index - 1)], 
                         na.rm = TRUE)
  
  # First increasing value after last 0-value
  change <- diff(df$Percentage_w_cros)
  increasing_indices <- which(change[zero_fall_index:length(change)] > 0)
  
  if (length(increasing_indices) == 0) {
    first_increasing_value <- NA
  } else {
    first_increasing_index <- zero_fall_index + min(increasing_indices, 
                                                    na.rm = TRUE)
    first_increasing_value <- df$Percentage_w_cros[first_increasing_index]
  }
  
  return(data.frame(Tree.ID = df$Tree.ID[1], 
                    Max_Before_Zero = max_before_zero, 
                    Last_Zero_Value = df$Percentage_w_cros[zero_fall_index], 
                    First_Increasing_Value = first_increasing_value))
}

# Function for calculating rolling statistics (mean, min, max) from 
# daily mean, min and max VPD over different timesteps (2-7 days)
calculate_rolling_stats <- function(data, columns, windows, funcs) {
  data %>%
    arrange(Tree.ID, Date) %>%
    group_by(Tree.ID) %>%
    mutate(across(
      all_of(columns),
      .fns = 
        lapply(windows, function(w) {
          lapply(funcs, function(f) {
            function(x) rollapply(x, width = w, FUN = f, align = "right", 
                                  fill = NA)
          })
        }) %>% 
        unlist(recursive = FALSE) %>%
        setNames(
          unlist(lapply(windows, function(w) {
            sapply(names(funcs), function(f) paste(w, f, sep = "_"))
          }))
        ),
      .names = "{.col}_{.fn}"
    )) %>%
    ungroup()
}

# Function for calculating correlation coefficients and p-values
calculate_correlations <- function(data, tree_species_name = "Overall", 
                                   treatment_name = "Overall") {
  cor_results <- rcorr(as.matrix(data), type = "spearman")
  
  # extract correlations and p-values
  cor_matrix <- as.data.frame(cor_results$r) %>%
    mutate(Tree.Species = tree_species_name, Treatment = treatment_name) %>%
    select(Tree.Species, Treatment, everything())
  
  p_matrix <- as.data.frame(cor_results$P) %>%
    mutate(Tree.Species = tree_species_name, Treatment = treatment_name) %>%
    select(Tree.Species, Treatment, everything())
  
  list(cor_matrix = cor_matrix, p_matrix = p_matrix)
}

# Function for filtering data for the different combinations of tree species and 
# treatments to calculate correlations for GH
calculate_all_correlations <- function(data, tree_species = NULL, 
                                       treatment = NULL) {
  if (!is.null(tree_species)) {
    data <- data %>% filter(Tree.Species == tree_species)
  }
  if (!is.null(treatment)) {
    data <- data %>% filter(Treatment == treatment)
  }
  data %>% select(c(5:8, 12, 19:ncol(data))) %>% na.omit() 
}

# Function for filtering data for the different combinations of tree species and 
# treatments to calculate correlations for EB 2023
calculate_all_correlations_EB_2023 <- function(data, tree_species = NULL, 
                                       treatment = NULL) {
  if (!is.null(tree_species)) {
    data <- data %>% filter(Tree.Species == tree_species)
  }
  if (!is.null(treatment)) {
    data <- data %>% filter(Treatment == treatment)
  }
  data %>% select(c(5:9, 19, 26:ncol(data))) %>% na.omit()
}

# Function for filtering data for the different combinations of tree species and 
# treatments to calculate correlations for EB 2024
calculate_all_correlations_EB_2024 <- function(data, tree_species = NULL, 
                                               treatment = NULL) {
  if (!is.null(tree_species)) {
    data <- data %>% filter(Tree.Species == tree_species)
  }
  if (!is.null(treatment)) {
    data <- data %>% filter(Treatment == treatment)
  }
  data %>% select(c(7:12, 21, 28:ncol(data))) %>% na.omit()
}

# Function for evaluating the models for explaining the irrigation demand
evaluate_irrigation_models_by_species_and_treatment <- function(data) {
  results <- data %>%
    group_by(Tree.Species, Treatment) %>%
    group_split() %>%
    map_dfr(~ {
      species <- unique(.x$Tree.Species)  # Extract species name
      treatment <- unique(.x$Treatment)  # Extract treatment name
      model_results <- evaluate_irrigation_models(.x) # apply function on subset
      
      # Save results as tibble
      tibble(
        Tree.Species = species,
        Treatment = treatment,
        AR2_Additive = model_results$R2_Additive,
        p_additive_swc = model_results$p_additive_swc, 
        p_additive_vpd = model_results$p_additive_vpd,
        AR2_Partial_VPD = model_results$R2_Partial_VPD,
        p_partial_VPD = model_results$p_partial_VPD,
        AR2_Partial_SWC = model_results$R2_Partial_SWC,
        p_partial_SWC = model_results$p_partial_SWC,
        VPD_Extra_Explanation = model_results$VPD_Extra_Explanation,
        SWC_Extra_Explanation = model_results$SWC_Extra_Explanation,
        p_value_model_additive = model_results$p_value_model_additive, 
        p_value_model_part_VPD = model_results$p_value_model_part_VPD,
        p_value_model_part_SWC = model_results$p_value_model_part_SWC,
        VIF_SWC = model_results$VIF_SWC,
        VIF_VPD = model_results$VIF_VPD
      )
    })
  return(results)
}

# Function for calculating multiple and partial linear regression between 
# irrigation demand and environmental drought stress variables
evaluate_irrigation_models <- function(data) {
  # multiple linear regression model
  model_perc <- lm(Percentage_w_cros ~ SWC_5_mean + VPD_adj_dailyavg_kPa_7_max, 
                   data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_swc <- summary(model_perc)$coefficients["SWC_5_mean", 4]  
  pvalue_additive_vpd <- 
    summary(model_perc)$coefficients["VPD_adj_dailyavg_kPa_7_max", 4]  
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive, 
                               df1 = summary(model_perc)$fstatistic[2], 
                               df2 = summary(model_perc)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # Partialmodel: Effect of VPD on irrig. demand without the influence of SWC
  model_perc_swc <- lm(Percentage_w_cros ~ SWC_5_mean, data = data)
  e.Perc_SWC <- residuals(model_perc_swc)
  
  model_VPD_SWC <- lm(VPD_adj_dailyavg_kPa_7_max ~ SWC_5_mean, data = data)
  e.VPD_SWC <- residuals(model_VPD_SWC)
  
  model_partial_VPD <- lm(e.Perc_SWC ~ e.VPD_SWC - 1)
  r2_partial_VPD <- summary(model_partial_VPD)$adj.r.squared
  pvalue_partial_VPD <- summary(model_partial_VPD)$coefficients[1, 4]  
  
  p_model_part_VPD <- summary(model_partial_VPD)$fstatistic[1]
  p_value_model_part_VPD <- pf(p_model_part_VPD, 
                               df1 = summary(model_partial_VPD)$fstatistic[2], 
                               df2 = summary(model_partial_VPD)$fstatistic[3], 
                               lower.tail = FALSE)
  
  
  # Partialmodel: Effect of SWC on irrig. demand without the influence of VPD
  model_perc_VPD <- lm(Percentage_w_cros ~ VPD_adj_dailyavg_kPa_7_max, 
                       data = data)
  e.Perc_VPD <- residuals(model_perc_VPD)
  
  model_SWC_VPD <- lm(SWC_5_mean ~ VPD_adj_dailyavg_kPa_7_max, data = data)
  e.SWC_VPD <- residuals(model_SWC_VPD)
  
  model_partial_SWC <- lm(e.Perc_VPD ~ e.SWC_VPD - 1)
  r2_partial_SWC <- summary(model_partial_SWC)$adj.r.squared
  pvalue_partial_SWC <- summary(model_partial_SWC)$coefficients[1, 4]  
  p_model_part_SWC <- summary(model_partial_SWC)$fstatistic[1]
  p_value_model_part_SWC <- pf(p_model_part_SWC, 
                               df1 = summary(model_partial_SWC)$fstatistic[2], 
                               df2 = summary(model_partial_SWC)$fstatistic[3], 
                               lower.tail = FALSE)
  
  
  # Difference between multiple linear regression  and partial regression 
  # for VPD
  VPD_explain_Perc <- r2_additive - r2_partial_SWC
  SWC_explain_Perc <- r2_additive - r2_partial_VPD
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_swc = pvalue_additive_swc, 
    p_additive_vpd = pvalue_additive_vpd,
    R2_Partial_VPD = r2_partial_VPD,
    p_partial_VPD = pvalue_partial_VPD,
    R2_Partial_SWC = r2_partial_SWC,
    p_partial_SWC = pvalue_partial_SWC,
    VPD_Extra_Explanation = VPD_explain_Perc,
    SWC_Extra_Explanation = SWC_explain_Perc,
    p_value_model_additive = p_value_model_additive, 
    p_value_model_part_VPD = p_value_model_part_VPD,
    p_value_model_part_SWC = p_value_model_part_SWC,
    VIF_VPD = vif_model_perc["VPD_adj_dailyavg_kPa_7_max"],
    VIF_SWC = vif_model_perc["SWC_5_mean"]
  ))
}

# Function for plotting overview of irrigation thresholds and environmental 
# conditions
plot_irrigation_warnings <- function(df, plotpath){
  # select if a rainshelter was present, if so, select rainshelter installation
  # and deinstallation date for plotting
  rs <- grepl("rs", df$Tree.ID[1])
  if(rs == TRUE){
    tree_id_short <- str_extract(df$Tree.ID[1], "D\\d+")
    filtered_meta <- rs_meta %>%
      filter(Dendro.ID == tree_id_short)
    rs_installation_date <- as.Date(filtered_meta$Installation_Date)
    rs_deinstallation_date <- as.Date(filtered_meta$Deinstallation_Date)
  }
  else{
    # invent a date for the installation of rs, 
    # that is not plotted to avoid errors, since no rs was installed
    rs_installation_date <- as.Date("2024-12-31")+1
    rs_deinstallation_date <- as.Date("2024-12-31")+1
  }
  
  # select irrigation dates according to treatment, location and year
  treat <- as.character(df$Treatment[1])
  location <- as.character(df$Site[1])
  year <- format(df$Date[1], "%Y") 
  # only Quercus robur hose was irrigated on 2023-07-12
  if (treat == "hose" & 
      location == "Eilaberg" & year == "2023" &
      filtered_data$Tree.Species[[1]] == "Quercus robur"){
    irrigation_date1 <- as.Date("2023-07-11")+1
    irrigation_date2 <- as.Date("2023-07-12")+1
    irrigation_date3 <- as.Date("2023-09-08")+1
    irrigation_date4 <- as.Date("2023-09-28")+1
  }
  else if (treat %in% c("drip", "hose") & 
           location == "Eilaberg" & year == "2023"){
    irrigation_date1 <- as.Date("2023-07-11")+1
    irrigation_date2 <- as.Date("2023-12-31")+1
    irrigation_date3 <- as.Date("2023-09-08")+1
    irrigation_date4 <- as.Date("2023-09-28")+1
  }
  
  else if (treat %in% c("drip", "hose") & 
           location == "Eilaberg" & year == "2024"){
    irrigation_date1 <- as.Date("2024-08-27")+1
    # only one irrigation date, so invent a date, that is not 
    # plotted to avoid errors
    irrigation_date2 <- as.Date("2024-12-31")+1
    irrigation_date3 <- as.Date("2024-12-31")+1
    irrigation_date4 <- as.Date("2024-12-31")+1
  }
  else if (treat == "drip" & 
           location == "Heiligenholz" & year == "2024"){
    irrigation_date1 <- as.Date("2024-08-12")+1
    irrigation_date2 <- as.Date("2024-08-29")+1
    # only two irrigation dates, so invent a date, that is not 
    # plotted to avoid errors
    irrigation_date3 <- as.Date("2024-12-31")+1
    irrigation_date4 <- as.Date("2024-12-31")+1
  }
  else if (treat == "hose" & 
           location == "Heiligenholz" & year == "2024"){
    irrigation_date1 <- as.Date("2024-08-29")+1
    # only one irrigation date, so invent dates, that is not 
    # plotted to avoid errors
    irrigation_date2 <- as.Date("2024-12-31")+1
    irrigation_date3 <- as.Date("2024-12-31")+1
    irrigation_date4 <- as.Date("2024-12-31")+1    
  }
  
  else {
    # no irrigation, so invent dates, that is not 
    # plotted to avoid errors
    irrigation_date1 <- as.Date("2024-12-31")+1
    irrigation_date2 <- as.Date("2024-12-31")+1
    irrigation_date3 <- as.Date("2024-12-31")+1
    irrigation_date4 <- as.Date("2024-12-31")+1
  }
  
  if (treat %in% c("drip", "hose") & 
      filtered_data$Tree.Species[[1]] %in% 
      c("Quercus robur", "Pseudotsuga menziesii")) {
    irrigation_date5 <- as.Date("2024-09-05")+1
  }
  
  else {
    irrigation_date5 <- as.Date("2024-12-31")+1
  }
  
  # plot with VPD and vol. SWC
  env_plot <- ggplot(df) +
    # add irrigation dates
    geom_vline(xintercept = as.POSIXct(irrigation_date1),
               color = "chartreuse4", size = 1) +
    geom_vline(xintercept = as.POSIXct(irrigation_date2),
               color = "chartreuse4", size = 1) +
    geom_vline(xintercept = as.POSIXct(irrigation_date3),
               color = "chartreuse4", size = 1) +
    geom_vline(xintercept = as.POSIXct(irrigation_date4),
               color = "chartreuse4", size = 1) +
    geom_vline(xintercept = as.POSIXct(irrigation_date5),
               color = "chartreuse4", size = 1) +
    
    # add rainshelter installation and deinstallation dates
    geom_vline(xintercept = as.POSIXct(rs_installation_date),
               color = "chocolate4", size = 1) +
    geom_vline(xintercept = as.POSIXct(rs_deinstallation_date),
               color = "chocolate4", size = 1) +
    
    geom_col(data = filter(climate_daily, Site == location), 
             aes(x = as.POSIXct(Date), y = Precip_mm_10min_sum/10),
             position = "stack", fill = "dodgerblue3") +
    geom_ribbon(data = filter(climate_daily, Site == location), 
                aes(x = as.POSIXct(Date),
                    ymin = VPD_kPA_min,
                    ymax = VPD_kPA_max),
                fill ="black", alpha = 0.2) +
    geom_line(data = filter(climate_daily, Site == location), 
              aes(x = as.POSIXct(Date), y = VPD_kPA_mean), 
              linetype = "dashed") +
    
    scale_y_continuous(
      name="VPD [kPA]",  
      limits = c(0, 6),
      sec.axis = sec_axis(~.*10, 
      name = "<span style='color:dodgerblue3;'>P [mm/day]</span> / <span style='color:blue;'>SWC [%]</span>")
    ) +
    geom_line(data = df,
              aes(x = as.POSIXct(Date), y = SWC_mean/10), color = "blue") +
    xlab(year) +
    scale_x_datetime(limits = x_limits, 
                     date_labels = "%d %b") +
    theme_bw() +
    theme(
      axis.title.y.right = element_markdown(size = 12),
      axis.text.y.right = element_text(color = "dodgerblue3")) 
  
  # Plot with jump-free and interpolated tree stem diameter 
  TSD_plot <- ggplot(df)+
    geom_line(aes(x = as.POSIXct(TIME), y = dm))+
    scale_x_datetime(limits = x_limits, 
                     date_labels = "%d %b") +
    labs(title = paste0(Tree.ID_Nr, " (", df$Tree.Species[1], 
                        " - ", df$Treatment[1], ")"),
         x = year, y = "TSD [µm]") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Plot with irrigation threshold ITlow and Growth stop
  TWD_plot <- ggplot(df) +
    # Highlight MinTWDnot0_Date (ITlow)
    geom_rect(data = filtered_data %>% filter(!is.na(MinTWDnot0_Date)),
              aes(xmin = as.POSIXct(TIME), xmax = as.POSIXct(TIME + 1440),
                  ymin = 0, ymax = Inf),
              fill = "darkgoldenrod2") +
    # Highlight first day after growth stopped
    geom_vline(data = filtered_data %>% 
                 filter(!is.na(FirstDayAfterLastGRO)) %>% slice(1),
               aes(xintercept = as.POSIXct(TIME)),
               color = "#4a2d9f", size = 2) +
    # Add TWD
    geom_line(aes(x = as.POSIXct(TIME), y = TWD)) +
    scale_x_datetime(limits = x_limits, 
                     date_labels = "%d %b") +
    labs(
      x = year,
      y = "TWD [µm]") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Plot with rel. TWD and percentage of reaching irrigation ITlow,
  # ITlow and ITup are highlighted as well
  perc_plot <- ggplot(df) +
    
    # add rectangle with percentage of reaching threshold B
    geom_rect(
      aes(xmin = as.POSIXct(Date),
          xmax = as.POSIXct(Date) + lubridate::days(1),
          ymin = -Inf, ymax = Inf, fill = Percentage_color),
      alpha = 0.4
    ) +
    scale_fill_identity() + 
    
    # Add 12% Xylem conductivity loss (ITup) threshold line
    geom_hline(yintercept = rel_TWD_max_thres, linetype = "dotted",
               color = "black") +
    
    # Add rel. TWD
    geom_line(aes(x = as.POSIXct(TIME), y = rel_TWD_log)) +
    scale_x_datetime(limits = x_limits, 
                     date_labels = "%d %b") +
    labs(x = "Date",
         y = "Rel. TWD") + 
    theme_bw() +
    ylim(0, 1) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  # combine plots
  combined_plot_area <-  TSD_plot / TWD_plot / perc_plot / env_plot + 
    plot_layout(ncol = 1, nrow = 4, heights = c(1, 1))
  
  combined_plot_area
  ggsave(path = plotpath,
         filename = paste0(Tree.ID_Nr, "_Irrigation_Thresholds.png"),
         width = 2560, height = 1489, units = "px") 
}

# function for evaluating the effect of Precip_mm_10min_sum_2_mean and
# VPD_kPA_mean on the irrigation demand of F. sylvatica (Part 1)
evaluate_irrigation_models_part_1_2023 <- function(data) {
  # Additive model
  model_perc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean + 
                     VPD_kPA_mean, data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_precip <- summary(model_perc)$coefficients["Precip_mm_10min_sum_2_mean", 4]  # p-value for precip
  pvalue_additive_vpd <- summary(model_perc)$coefficients["VPD_kPA_mean", 4]  # p-value for VPD
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive, 
                               df1 = summary(model_perc)$fstatistic[2], 
                               df2 = summary(model_perc)$fstatistic[3], 
                               lower.tail = FALSE)
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # partial model: effect of precip on Percentage_w_cros without effect of VPD 
  model_perc_from_vpd <- lm(Percentage_w_cros ~ VPD_kPA_mean , data = data)
  e.perc_from_vpd <- residuals(model_perc_from_vpd)
  
  model_precip_from_VPD<- lm(Precip_mm_10min_sum_2_mean ~ VPD_kPA_mean, data = data)
  e.precip_from_VPD <- residuals(model_precip_from_VPD)
  
  model_partial_precip <- lm(e.perc_from_vpd ~ e.precip_from_VPD- 1)
  r2_partial_precip <- summary(model_partial_precip)$adj.r.squared
  p_value_partial_precip <- summary(model_partial_precip)$coefficients[1, 4]  # p-value of VPD in partial model
  
  p_model_part_precip <- summary(model_partial_precip)$fstatistic[1]
  p_value_model_part_precip <- pf(p_model_part_precip, 
                                  df1 = summary(model_partial_precip)$fstatistic[2], 
                                  df2 = summary(model_partial_precip)$fstatistic[3], 
                                  lower.tail = FALSE)
  
  # partial model: effect of VPD on Percentage_w_cros without effect of precip 
  # and SWC
  model_perc_from_precip <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean, data = data)
  e.perc_from_precip <- residuals(model_perc_from_precip)
  
  model_VPD_from_precip <- lm(VPD_kPA_mean ~ Precip_mm_10min_sum_2_mean, data = data)
  e.VPD_from_precip <- residuals(model_VPD_from_precip)
  
  model_partial_vpd <- lm(e.perc_from_precip ~ e.VPD_from_precip - 1)
  r2_partial_vpd <- summary(model_partial_vpd)$adj.r.squared
  p_value_partial_vpd <- summary(model_partial_vpd)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_vpd <- summary(model_partial_vpd)$fstatistic[1]
  p_value_model_part_vpd <- pf(p_model_part_vpd, 
                               df1 = summary(model_partial_vpd)$fstatistic[2], 
                               df2 = summary(model_partial_vpd)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # Calculate difference between additive and partial model for VPD
  VPD_explain_Perc <- r2_additive - r2_partial_precip
  precip_explain_Perc <- r2_additive - r2_partial_vpd
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_precip = pvalue_additive_precip, 
    p_additive_vpd = pvalue_additive_vpd,
    R2_Partial_VPD = r2_partial_vpd,
    p_partial_VPD = p_value_partial_vpd,
    R2_Partial_precip = r2_partial_precip,
    p_partial_precip = p_value_partial_precip,
    VPD_Extra_Explanation = VPD_explain_Perc,
    Precip_Extra_Explanation = precip_explain_Perc,
    p_value_model_additive = p_value_model_additive, 
    p_value_model_part_VPD = p_value_model_part_vpd,
    p_value_model_part_precip = p_value_model_part_precip,
    VIF_Precip = vif_model_perc["Precip_mm_10min_sum_2_mean"],
    VIF_VPD = vif_model_perc["VPD_kPA_mean"]
  ))
}

# function for applying evaluate_irrigation_models_part_1_2023 on a subset of
# data 
evaluate_irrigation_models_by_species_part_1_2023 <- function(data) {
  results <- data %>%
    group_by(Tree.Species) %>%
    group_split() %>%
    map_dfr(~ {
      species <- unique(.x$Tree.Species)  
      # apply function on subset
      model_results <- evaluate_irrigation_models_part_1_2023(.x)  
      tibble(
        Tree.Species = species,
        Treatment = "Overall",
        AR2_Additive = model_results$R2_Additive,
        p_additive_precip = model_results$p_additive_precip, 
        p_additive_vpd = model_results$p_additive_vpd,
        AR2_Partial_VPD = model_results$R2_Partial_VPD,
        p_partial_VPD = model_results$p_partial_VPD,
        AR2_Partial_precip = model_results$R2_Partial_precip,
        p_partial_precip = model_results$p_partial_precip,
        VPD_Extra_Explanation = model_results$VPD_Extra_Explanation,
        Precip_Extra_Explanation = model_results$Precip_Extra_Explanation,
        p_value_model_additive = model_results$p_value_model_additive, 
        p_value_model_part_VPD = model_results$p_value_model_part_VPD,
        p_value_model_part_precip = model_results$p_value_model_part_precip,
        VIF_Precip = model_results$VIF_Precip,
        VIF_VPD = model_results$VIF_VPD
      )
    })
  return(results)
}

# function for evaluating the effect of Precip_mm_10min_sum_2_mean and
# VPD_kPA_mean_2_min on the irrigation demand of Q. robur (Part 2)
evaluate_irrigation_models_part_2_2023 <- function(data) {
  # additives model
  model_perc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean + VPD_kPA_mean_2_min, data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_precip <- summary(model_perc)$coefficients["Precip_mm_10min_sum_2_mean", 4]  # p-value for precip
  pvalue_additive_vpd <- summary(model_perc)$coefficients["VPD_kPA_mean_2_min", 4]  # p-value for VPD
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive,
                               df1 = summary(model_perc)$fstatistic[2],
                               df2 = summary(model_perc)$fstatistic[3],
                               lower.tail = FALSE)
  
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # partial model: effect of precip Percentage_w_cros without the effect of VPD
  model_perc_precip <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean, data = data)
  e.Perc_precip <- residuals(model_perc_precip)
  
  model_VPD_precip <- lm(VPD_kPA_mean_2_min ~ Precip_mm_10min_sum_2_mean, data = data)
  e.VPD_precip <- residuals(model_VPD_precip)
  
  model_partial_VPD <- lm(e.Perc_precip ~ e.VPD_precip - 1)
  r2_partial_VPD <- summary(model_partial_VPD)$adj.r.squared
  pvalue_partial_VPD <- summary(model_partial_VPD)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_VPD <- summary(model_partial_VPD)$fstatistic[1]
  p_value_model_part_VPD <- pf(p_model_part_VPD,
                               df1 = summary(model_partial_VPD)$fstatistic[2],
                               df2 = summary(model_partial_VPD)$fstatistic[3],
                               lower.tail = FALSE)
  
  # Partial model: effect of VPD on Percentage_w_cros without precip
  model_perc_VPD <- lm(Percentage_w_cros ~ VPD_kPA_mean_2_min, data = data)
  e.Perc_VPD <- residuals(model_perc_VPD)
  
  model_precip_VPD <- lm(Precip_mm_10min_sum_2_mean ~ VPD_kPA_mean_2_min, data = data)
  e.precip_VPD <- residuals(model_precip_VPD)
  
  model_partial_precip <- lm(e.Perc_VPD ~ e.precip_VPD - 1)
  r2_partial_precip <- summary(model_partial_precip)$adj.r.squared
  pvalue_partial_precip <- summary(model_partial_precip)$coefficients[1, 4]  # p-value for precip in partial model
  p_model_part_precip <- summary(model_partial_precip)$fstatistic[1]
  p_value_model_part_precip <- pf(p_model_part_precip,
                                  df1 = summary(model_partial_precip)$fstatistic[2],
                                  df2 = summary(model_partial_precip)$fstatistic[3],
                                  lower.tail = FALSE)
  
  # Difference between additive and partial model for VPD
  VPD_explain_Perc <- r2_additive - r2_partial_precip
  precip_explain_Perc <- r2_additive - r2_partial_VPD
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_precip = pvalue_additive_precip,
    p_additive_vpd = pvalue_additive_vpd,
    R2_Partial_VPD = r2_partial_VPD,
    p_partial_VPD = pvalue_partial_VPD,
    R2_Partial_precip = r2_partial_precip,
    p_partial_precip = pvalue_partial_precip,
    VPD_Extra_Explanation = VPD_explain_Perc,
    Precip_Extra_Explanation = precip_explain_Perc,
    p_value_model_additive = p_value_model_additive,
    p_value_model_part_VPD = p_value_model_part_VPD,
    p_value_model_part_precip = p_value_model_part_precip,
    VIF_Precip = vif_model_perc["Precip_mm_10min_sum_2_mean"],
    VIF_VPD = vif_model_perc["VPD_kPA_mean_2_min"]
  ))
}

# function for applying evaluate_irrigation_models_part_2_2023 on a subset of
# data
evaluate_irrigation_models_by_species_and_treatment_part_2_2023 <- function(data) {
  results <- data %>%
    group_by(Tree.Species) %>%
    group_split() %>%
    map_dfr(~ {
      species <- unique(.x$Tree.Species)  
      # Apply function on subset
      model_results <- evaluate_irrigation_models_part_2_2023(.x)  
      
      tibble(
        Tree.Species = species,
        Treatment = "Overall",
        AR2_Additive = model_results$R2_Additive,
        p_additive_precip = model_results$p_additive_precip,
        p_additive_vpd = model_results$p_additive_vpd,
        AR2_Partial_VPD = model_results$R2_Partial_VPD,
        p_partial_VPD = model_results$p_partial_VPD,
        AR2_Partial_precip = model_results$R2_Partial_precip,
        p_partial_precip = model_results$p_partial_precip,
        VPD_Extra_Explanation = model_results$VPD_Extra_Explanation,
        Precip_Extra_Explanation = model_results$Precip_Extra_Explanation,
        p_value_model_additive = model_results$p_value_model_additive,
        p_value_model_part_VPD = model_results$p_value_model_part_VPD,
        p_value_model_part_precip = model_results$p_value_model_part_precip,
        VIF_Precip = model_results$VIF_Precip,
        VIF_VPD = model_results$VIF_VPD
      )
    })
  return(results)
}


# function for evaluating the effect of Precip_mm_10min_sum_2_mean and
# VPD_kPA_mean_2_min on the irrigation demand of both tree species (Part 3)
evaluate_irrigation_models_part_3_2023 <- function(data) {
  # Additive model
  model_perc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean + VPD_kPA_mean_2_min, data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_precip <- summary(model_perc)$coefficients["Precip_mm_10min_sum_2_mean", 4]  # p-value for precip
  pvalue_additive_vpd <- summary(model_perc)$coefficients["VPD_kPA_mean_2_min", 4]  # p-value for VPD
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive,
                               df1 = summary(model_perc)$fstatistic[2],
                               df2 = summary(model_perc)$fstatistic[3],
                               lower.tail = FALSE)
  
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # partial model: effect of precip on Percentage_w_cros without influence of VPD
  model_perc_precip <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean, data = data)
  e.Perc_precip <- residuals(model_perc_precip)
  
  model_VPD_precip <- lm(VPD_kPA_mean_2_min ~ Precip_mm_10min_sum_2_mean, data = data)
  e.VPD_precip <- residuals(model_VPD_precip)
  
  model_partial_VPD <- lm(e.Perc_precip ~ e.VPD_precip - 1)
  r2_partial_VPD <- summary(model_partial_VPD)$adj.r.squared
  pvalue_partial_VPD <- summary(model_partial_VPD)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_VPD <- summary(model_partial_VPD)$fstatistic[1]
  p_value_model_part_VPD <- pf(p_model_part_VPD,
                               df1 = summary(model_partial_VPD)$fstatistic[2],
                               df2 = summary(model_partial_VPD)$fstatistic[3],
                               lower.tail = FALSE)
  
  
  # partial model: effect of VPD on Percentage_w_cros without the effect of 
  # precip
  model_perc_VPD <- lm(Percentage_w_cros ~ VPD_kPA_mean_2_min, data = data)
  e.Perc_VPD <- residuals(model_perc_VPD)
  
  model_precip_VPD <- lm(Precip_mm_10min_sum_2_mean ~ VPD_kPA_mean_2_min, data = data)
  e.precip_VPD <- residuals(model_precip_VPD)
  
  model_partial_precip <- lm(e.Perc_VPD ~ e.precip_VPD - 1)
  r2_partial_precip <- summary(model_partial_precip)$adj.r.squared
  pvalue_partial_precip <- summary(model_partial_precip)$coefficients[1, 4]  # p-value for precip in partial model
  p_model_part_precip <- summary(model_partial_precip)$fstatistic[1]
  p_value_model_part_precip <- pf(p_model_part_precip,
                                  df1 = summary(model_partial_precip)$fstatistic[2],
                                  df2 = summary(model_partial_precip)$fstatistic[3],
                                  lower.tail = FALSE)
  
  
  # difference between additive and partial model for VPD
  VPD_explain_Perc <- r2_additive - r2_partial_precip
  precip_explain_Perc <- r2_additive - r2_partial_VPD
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_precip = pvalue_additive_precip,
    p_additive_vpd = pvalue_additive_vpd,
    R2_Partial_VPD = r2_partial_VPD,
    p_partial_VPD = pvalue_partial_VPD,
    R2_Partial_precip = r2_partial_precip,
    p_partial_precip = pvalue_partial_precip,
    VPD_Extra_Explanation = VPD_explain_Perc,
    Precip_Extra_Explanation = precip_explain_Perc,
    p_value_model_additive = p_value_model_additive,
    p_value_model_part_VPD = p_value_model_part_VPD,
    p_value_model_part_precip = p_value_model_part_precip,
    VIF_Precip = vif_model_perc["Precip_mm_10min_sum_2_mean"],
    VIF_VPD = vif_model_perc["VPD_kPA_mean_2_min"]
  ))
}

# function for applying evaluate_irrigation_models_part_3_2023 on a subset of
# data 
evaluate_irrigation_models_overall_part_3_2023 <- function(data) {
  model_results <- evaluate_irrigation_models_part_3_2023(data)    
  tibble(
    Tree.Species = species,
    Treatment = "Overall",
    AR2_Additive = model_results$R2_Additive,
    p_additive_precip = model_results$p_additive_precip,
    p_additive_vpd = model_results$p_additive_vpd,
    AR2_Partial_VPD = model_results$R2_Partial_VPD,
    p_partial_VPD = model_results$p_partial_VPD,
    AR2_Partial_precip = model_results$R2_Partial_precip,
    p_partial_precip = model_results$p_partial_precip,
    VPD_Extra_Explanation = model_results$VPD_Extra_Explanation,
    Precip_Extra_Explanation = model_results$Precip_Extra_Explanation,
    p_value_model_additive = model_results$p_value_model_additive,
    p_value_model_part_VPD = model_results$p_value_model_part_VPD,
    p_value_model_part_precip = model_results$p_value_model_part_precip,
    VIF_Precip = model_results$VIF_Precip,
    VIF_VPD = model_results$VIF_VPD
  )
}

# function for evaluating the effect of Precip_mm_10min_sum_3_mean, SWC_mean
# and VPD_kPA_min_2_min on the irrigation demand of F. sylvatica (Part 1) 
# for 2024
evaluate_irrigation_models_part_1_2024 <- function(data) {
  # Additive model
  model_perc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_3_mean + 
                     VPD_kPA_min_2_min + SWC_mean, data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_precip <- summary(model_perc)$coefficients["Precip_mm_10min_sum_3_mean", 4]  # p-value for precip
  pvalue_additive_vpd <- summary(model_perc)$coefficients["VPD_kPA_min_2_min", 4]  # p-value for VPD
  pvalue_additive_swc <- summary(model_perc)$coefficients["SWC_mean", 4]
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive, 
                               df1 = summary(model_perc)$fstatistic[2], 
                               df2 = summary(model_perc)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # partial model: effect of precip on Percentage_w_cros without effect of VPD and SWC
  model_perc_from_vpd_swc <- lm(Percentage_w_cros ~ 
                                  VPD_kPA_min_2_min + SWC_mean, data = data)
  e.perc_from_vpd_swc <- residuals(model_perc_from_vpd_swc)
  
  model_precip_from_VPD_swc <- lm(Precip_mm_10min_sum_3_mean ~ 
                                    VPD_kPA_min_2_min + SWC_mean, data = data)
  e.precip_from_VPD_swc <- residuals(model_precip_from_VPD_swc)
  
  model_partial_precip <- lm(e.perc_from_vpd_swc ~ e.precip_from_VPD_swc - 1)
  r2_partial_precip <- summary(model_partial_precip)$adj.r.squared
  p_value_partial_precip <- summary(model_partial_precip)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_precip <- summary(model_partial_precip)$fstatistic[1]
  p_value_model_part_precip <- pf(p_model_part_precip, 
                                  df1 = summary(model_partial_precip)$fstatistic[2], 
                                  df2 = summary(model_partial_precip)$fstatistic[3], 
                                  lower.tail = FALSE)
  
  # partial model: effect of VPD on Percentage_w_cros without the effect of
  # precip and SWC
  model_perc_from_precip_swc <- lm(Percentage_w_cros ~ 
                                     Precip_mm_10min_sum_3_mean + SWC_mean, 
                                   data = data)
  e.perc_from_precip_swc <- residuals(model_perc_from_precip_swc)
  
  model_VPD_from_precip_swc <- lm(VPD_kPA_min_2_min ~ 
                                    Precip_mm_10min_sum_3_mean + SWC_mean, 
                                  data = data)
  e.VPD_from_precip_swc <- residuals(model_VPD_from_precip_swc)
  
  model_partial_vpd <- lm(e.perc_from_precip_swc ~ e.VPD_from_precip_swc - 1)
  r2_partial_vpd <- summary(model_partial_vpd)$adj.r.squared
  p_value_partial_vpd <- summary(model_partial_vpd)$coefficients[1, 4]  # p-value for VPD in partialmodel
  
  p_model_part_vpd <- summary(model_partial_vpd)$fstatistic[1]
  p_value_model_part_vpd <- pf(p_model_part_vpd, 
                               df1 = summary(model_partial_vpd)$fstatistic[2], 
                               df2 = summary(model_partial_vpd)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # Partial model: effect of SWC on Percentage_w_cros without the effects of
  # precip and VPD
  model_perc_from_vpd_precip <- lm(Percentage_w_cros ~ VPD_kPA_min_2_min + 
                                     Precip_mm_10min_sum_3_mean, data = data)
  e.perc_from_vpd_precip <- residuals(model_perc_from_vpd_precip)
  
  model_swc_from_VPD_precip <- lm(SWC_mean ~ Precip_mm_10min_sum_3_mean + 
                                    VPD_kPA_min_2_min , data = data)
  e.swc_from_VPD_precip <- residuals(model_swc_from_VPD_precip)
  
  model_partial_swc <- lm(e.perc_from_vpd_precip ~ e.swc_from_VPD_precip - 1)
  
  r2_partial_swc <- summary(model_partial_swc)$adj.r.squared
  p_value_partial_swc <- summary(model_partial_swc)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_swc <- summary(model_partial_swc)$fstatistic[1]
  p_value_model_part_swc <- pf(p_model_part_swc, 
                               df1 = summary(model_partial_swc)$fstatistic[2], 
                               df2 = summary(model_partial_swc)$fstatistic[3], 
                               lower.tail = FALSE)
  
  
  # Calculate difference between additive and partial model for VPD
  VPD_SWC_explain_Perc <- r2_additive - r2_partial_precip
  precip_SWC_explain_Perc <- r2_additive - r2_partial_vpd
  precip_VPD_explain_Perc <- r2_additive - r2_partial_swc
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_precip = pvalue_additive_precip, 
    p_additive_vpd = pvalue_additive_vpd,
    p_additive_swc = pvalue_additive_swc,
    R2_Partial_VPD = r2_partial_vpd,
    p_partial_VPD = p_value_partial_vpd,
    R2_Partial_precip = r2_partial_precip,
    p_partial_precip = p_value_partial_precip,
    R2_Partial_swc = r2_partial_swc,
    p_partial_swc = p_value_partial_swc,
    VPD_SWC_Extra_Explanation = VPD_SWC_explain_Perc,
    Precip_SWC_Extra_Explanation = precip_SWC_explain_Perc,
    Precip_VPD_Extra_Explanation = precip_VPD_explain_Perc,
    p_value_model_additive = p_value_model_additive, 
    p_value_model_part_VPD = p_value_model_part_vpd,
    p_value_model_part_precip = p_value_model_part_precip,
    p_value_model_part_swc = p_value_model_part_swc,
    VIF_Precip = vif_model_perc["Precip_mm_10min_sum_3_mean"],
    VIF_VPD = vif_model_perc["VPD_kPA_min_2_min"],
    VIF_SWC = vif_model_perc["SWC_mean"]
  ))
}

# function for applying evaluate_irrigation_models_part_1_2024 on a subset of
# data 
evaluate_irrigation_models_by_species_part_1_2024 <- function(data) {
  results <- data %>%
    group_by(Tree.Species) %>%
    group_split() %>%
    map_dfr(~ {
      species <- unique(.x$Tree.Species)  
      # apply function to subset
      model_results <- evaluate_irrigation_models_part_1_2024(.x)  
      
      tibble(
        Tree.Species = species,
        Treatment = treatment,
        AR2_Additive = model_results$R2_Additive,
        p_additive_precip = model_results$p_additive_precip, 
        p_additive_vpd = model_results$p_additive_vpd,
        p_additive_swc = model_results$p_additive_swc,
        AR2_Partial_VPD = model_results$R2_Partial_VPD,
        p_partial_VPD = model_results$p_partial_VPD,
        AR2_Partial_precip = model_results$R2_Partial_precip,
        p_partial_precip = model_results$p_partial_precip,
        AR2_Partial_swc = model_results$R2_Partial_swc,
        p_partial_swc = model_results$p_partial_swc,
        VPD_SWC_Extra_Explanation = model_results$VPD_SWC_Extra_Explanation,
        Precip_SWC_Extra_Explanation = model_results$Precip_SWC_Extra_Explanation,
        Precip_VPD_Extra_Explanation = model_results$Precip_VPD_Extra_Explanation ,
        p_value_model_additive = model_results$p_value_model_additive, 
        p_value_model_part_VPD = model_results$p_value_model_part_VPD,
        p_value_model_part_precip = model_results$p_value_model_part_precip,
        p_value_model_part_swc = model_results$p_value_model_part_swc,
        VIF_Precip = model_results$VIF_Precip,
        VIF_VPD = model_results$VIF_VPD,
        VIF_SWC = model_results$VIF_SWC
      )
    })
  return(results)
}

# function for evaluating the effect of Precip_mm_10min_sum_2_mean, SWC_mean
# and VPD_kPA_min_2_min on the irrigation demand of Q. robur (Part 1) 
# for 2024 
evaluate_irrigation_models_part_2_2024 <- function(data) {
  # Additive model
  model_perc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean + 
                     VPD_kPA_min_2_min + SWC_mean, data = data)
  r2_additive <- summary(model_perc)$adj.r.squared
  pvalue_additive_precip <- summary(model_perc)$coefficients["Precip_mm_10min_sum_2_mean", 4]  # p-value for precip
  pvalue_additive_vpd <- summary(model_perc)$coefficients["VPD_kPA_min_2_min", 4]  # p-value for VPD
  pvalue_additive_swc <- summary(model_perc)$coefficients["SWC_mean", 4]
  p_model_additive <- summary(model_perc)$fstatistic[1]
  p_value_model_additive <- pf(p_model_additive, 
                               df1 = summary(model_perc)$fstatistic[2], 
                               df2 = summary(model_perc)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # check multicolinearity -> Variance Inflation Factor (VIF)
  vif_model_perc <- vif(model_perc) 
  print(vif_model_perc)
  
  # Partial model: effect of precip on Percentage_w_cros without the effect of
  # VPD and SWC
  model_perc_from_vpd_swc <- lm(Percentage_w_cros ~ VPD_kPA_min_2_min + 
                                  SWC_mean, data = data)
  e.perc_from_vpd_swc <- residuals(model_perc_from_vpd_swc)
  
  model_precip_from_VPD_swc <- lm(Precip_mm_10min_sum_2_mean ~ 
                                    VPD_kPA_min_2_min + SWC_mean, data = data)
  e.precip_from_VPD_swc <- residuals(model_precip_from_VPD_swc)
  
  model_partial_precip <- lm(e.perc_from_vpd_swc ~ e.precip_from_VPD_swc - 1)
  r2_partial_precip <- summary(model_partial_precip)$adj.r.squared
  p_value_partial_precip <- summary(model_partial_precip)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_precip <- summary(model_partial_precip)$fstatistic[1]
  p_value_model_part_precip <- pf(p_model_part_precip, 
                                  df1 = summary(model_partial_precip)$fstatistic[2], 
                                  df2 = summary(model_partial_precip)$fstatistic[3], 
                                  lower.tail = FALSE)
  
  # Partial model: effect of VPD on Percentage_w_cros without the effect of
  # precip and SWC
  model_perc_from_precip_swc <- lm(Percentage_w_cros ~ Precip_mm_10min_sum_2_mean + SWC_mean, data = data)
  e.perc_from_precip_swc <- residuals(model_perc_from_precip_swc)
  
  model_VPD_from_precip_swc <- lm(VPD_kPA_min_2_min ~ Precip_mm_10min_sum_2_mean + SWC_mean, data = data)
  e.VPD_from_precip_swc <- residuals(model_VPD_from_precip_swc)
  
  model_partial_vpd <- lm(e.perc_from_precip_swc ~ e.VPD_from_precip_swc - 1)
  r2_partial_vpd <- summary(model_partial_vpd)$adj.r.squared
  p_value_partial_vpd <- summary(model_partial_vpd)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_vpd <- summary(model_partial_vpd)$fstatistic[1]
  p_value_model_part_vpd <- pf(p_model_part_vpd, 
                               df1 = summary(model_partial_vpd)$fstatistic[2], 
                               df2 = summary(model_partial_vpd)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # Partial model: effect of SWC on Percentage_w_cros without the effect of 
  # precip and VPD
  model_perc_from_vpd_precip <- lm(Percentage_w_cros ~ VPD_kPA_min_2_min + 
                                     Precip_mm_10min_sum_2_mean, data = data)
  e.perc_from_vpd_precip <- residuals(model_perc_from_vpd_precip)
  
  model_swc_from_VPD_precip <- lm(SWC_mean ~ Precip_mm_10min_sum_2_mean + 
                                    VPD_kPA_min_2_min , data = data)
  e.swc_from_VPD_precip <- residuals(model_swc_from_VPD_precip)
  
  model_partial_swc <- lm(e.perc_from_vpd_precip ~ e.swc_from_VPD_precip - 1)
  
  r2_partial_swc <- summary(model_partial_swc)$adj.r.squared
  p_value_partial_swc <- summary(model_partial_swc)$coefficients[1, 4]  # p-value for VPD in partial model
  
  p_model_part_swc <- summary(model_partial_swc)$fstatistic[1]
  p_value_model_part_swc <- pf(p_model_part_swc, 
                               df1 = summary(model_partial_swc)$fstatistic[2], 
                               df2 = summary(model_partial_swc)$fstatistic[3], 
                               lower.tail = FALSE)
  
  # Difference between additive model and partial model for VPD
  VPD_SWC_explain_Perc <- r2_additive - r2_partial_precip
  precip_SWC_explain_Perc <- r2_additive - r2_partial_vpd
  precip_VPD_explain_Perc <- r2_additive - r2_partial_swc
  
  return(list(
    R2_Additive = r2_additive,
    p_additive_precip = pvalue_additive_precip, 
    p_additive_vpd = pvalue_additive_vpd,
    p_additive_swc = pvalue_additive_swc,
    R2_Partial_VPD = r2_partial_vpd,
    p_partial_VPD = p_value_partial_vpd,
    R2_Partial_precip = r2_partial_precip,
    p_partial_precip = p_value_partial_precip,
    R2_Partial_swc = r2_partial_swc,
    p_partial_swc = p_value_partial_swc,
    VPD_SWC_Extra_Explanation = VPD_SWC_explain_Perc,
    Precip_SWC_Extra_Explanation = precip_SWC_explain_Perc,
    Precip_VPD_Extra_Explanation = precip_VPD_explain_Perc,
    p_value_model_additive = p_value_model_additive, 
    p_value_model_part_VPD = p_value_model_part_vpd,
    p_value_model_part_precip = p_value_model_part_precip,
    p_value_model_part_swc = p_value_model_part_swc,
    VIF_Precip = vif_model_perc["Precip_mm_10min_sum_2_mean"],
    VIF_VPD = vif_model_perc["VPD_kPA_min_2_min"],
    VIF_SWC = vif_model_perc["SWC_mean"]
  ))
}

# function for applying evaluate_irrigation_models_part_2_2024 on a subset of
# data
evaluate_irrigation_models_by_species_part_2_2024 <- function(data) {
  results <- data %>%
    group_by(Tree.Species) %>%
    group_split() %>%
    map_dfr(~ {
      species <- unique(.x$Tree.Species)
      # apply function on subset
      model_results <- evaluate_irrigation_models_part_2_2024(.x)  
      
      tibble(
        Tree.Species = species,
        Treatment = treatment,
        AR2_Additive = model_results$R2_Additive,
        p_additive_precip = model_results$p_additive_precip, 
        p_additive_vpd = model_results$p_additive_vpd,
        p_additive_swc = model_results$p_additive_swc,
        AR2_Partial_VPD = model_results$R2_Partial_VPD,
        p_partial_VPD = model_results$p_partial_VPD,
        AR2_Partial_precip = model_results$R2_Partial_precip,
        p_partial_precip = model_results$p_partial_precip,
        AR2_Partial_swc = model_results$R2_Partial_swc,
        p_partial_swc = model_results$p_partial_swc,
        VPD_SWC_Extra_Explanation = model_results$VPD_SWC_Extra_Explanation,
        Precip_SWC_Extra_Explanation = model_results$Precip_SWC_Extra_Explanation,
        Precip_VPD_Extra_Explanation = model_results$Precip_VPD_Extra_Explanation ,
        p_value_model_additive = model_results$p_value_model_additive, 
        p_value_model_part_VPD = model_results$p_value_model_part_VPD,
        p_value_model_part_precip = model_results$p_value_model_part_precip,
        p_value_model_part_swc = model_results$p_value_model_part_swc,
        VIF_Precip = model_results$VIF_Precip,
        VIF_VPD = model_results$VIF_VPD,
        VIF_SWC = model_results$VIF_SWC
      )
    })
  return(results)
}

# function for applying evaluate_irrigation_models_part_2_2024 on a subset of
# data but both species
evaluate_irrigation_models_overall_part_2_2024 <- function(data) {
  model_results <- evaluate_irrigation_models_part_2_2024(data)
  tibble(
    Tree.Species = species,
    Treatment = treatment,
    AR2_Additive = model_results$R2_Additive,
    p_additive_precip = model_results$p_additive_precip, 
    p_additive_vpd = model_results$p_additive_vpd,
    p_additive_swc = model_results$p_additive_swc,
    AR2_Partial_VPD = model_results$R2_Partial_VPD,
    p_partial_VPD = model_results$p_partial_VPD,
    AR2_Partial_precip = model_results$R2_Partial_precip,
    p_partial_precip = model_results$p_partial_precip,
    AR2_Partial_swc = model_results$R2_Partial_swc,
    p_partial_swc = model_results$p_partial_swc,
    VPD_SWC_Extra_Explanation = model_results$VPD_SWC_Extra_Explanation,
    Precip_SWC_Extra_Explanation = model_results$Precip_SWC_Extra_Explanation,
    Precip_VPD_Extra_Explanation = model_results$Precip_VPD_Extra_Explanation ,
    p_value_model_additive = model_results$p_value_model_additive, 
    p_value_model_part_VPD = model_results$p_value_model_part_VPD,
    p_value_model_part_precip = model_results$p_value_model_part_precip,
    p_value_model_part_swc = model_results$p_value_model_part_swc,
    VIF_Precip = model_results$VIF_Precip,
    VIF_VPD = model_results$VIF_VPD,
    VIF_SWC = model_results$VIF_SWC
  )
}
