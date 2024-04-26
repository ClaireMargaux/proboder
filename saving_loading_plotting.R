#' Load Data
#'
#' Load simulated or real data based on the specified type, region, and frequency.
#'
#' @param type Character, type of data to load ('simulated_LSODA', 'simulated_HETTMO' or 'real').
#' @param region Character, region for which to load the data ('BE' or 'GE').
#' @param daily_or_weekly Character, frequency of the data ('daily' or 'weekly').
#' @param directory Character, directory path where the data is stored.
#'
#' @return A list containing the loaded data including observations, population, and real beta (if available).
#'
#' @examples
#' load_data("simulated_LSODA", "BE", "daily", "/path/to/data/directory/LSODA")
#' load_data("real", "GE", "weekly", "/path/to/data/directory/real")
#'
#' @export
load_data <- function(type, region, daily_or_weekly, directory) {
  # Check if type is valid
  if (!(type %in% c('simulated_LSODA', 'simulated_HETTMO', 'real'))) {
    stop("Error: Invalid type. Type must be 'simulated' or 'real'.")
  }
  
  # Check if region is valid
  if (!(region %in% c('BE', 'GE'))) {
    stop("Error: Invalid region. Region must be 'BE' or 'GE'.")
  }
  
  # Check if daily_or_weekly is valid
  if (!(daily_or_weekly %in% c('daily', 'weekly'))) {
    stop("Error: Invalid value for daily_or_weekly. Must be 'daily' or 'weekly'.")
  }
  
  # Load data based on type
  if (type == 'simulated_LSODA') {
    load(file.path(directory, "simulated_data_LSODA.Rdata"))
    obs <- df
    load(file.path(directory, "simulated_noisy_data_LSODA.Rdata"))
    obs_with_noise <- df_with_noise
    load(file.path(directory, "simulated_params_LSODA.Rdata"))
    params <- df_params
    real_beta <- readRDS(file.path(directory, "simulated_beta_LSODA.Rds"))
    return <- list(obs = obs, obs_with_noise = obs_with_noise, params = params, real_beta = real_beta)
  }
  else if (type == 'simulated_HETTMO') {
    load(file.path(directory, "simulated_data_HETTMO.Rdata"))
    obs <- observations_simulate
    load(file.path(directory, "simulated_params_HETTMO.Rdata"))
    params <- df_params
    real_beta <- readRDS(file.path(directory, "simulated_beta_HETTMO.Rds"))
    return <- list(obs = obs, params = params, real_beta = real_beta)
  }
  else if (type == 'real') {
    region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory, region_filename))
    obs <- observations
    population <- readRDS(file.path(directory, population_filename))
    real_beta <- NULL 
    return <- list(obs = obs, params = params, real_beta = real_beta)
  }

  # Return loaded data
  return(return)
}

#' Save Results as .RData Files
#'
#' This function saves X_values, U_values, P_X_values, and P_U_values as .Rdata files.
#'
#' @param X_values Vector to be saved as an .RData file.
#' @param U_values Vector to be saved as an .RData file.
#' @param P_X_values Vector to be saved as an .RData file.
#' @param P_U_values Vector to be saved as an .RData file.
#' @param directory Character string indicating the directory path where the data files will be saved.
#' @export
save_results_as_Rdata <- function(X_values, U_values, P_X_values, P_U_values, directory) {
  # Save X_values, U_values, P_X_values, and P_U_values as .Rdata files
  save(X_values, file = paste0(directory, "/X_values.Rdata"))
  save(U_values, file = paste0(directory, "/U_values.Rdata"))
  save(P_X_values, file = paste0(directory, "/P_X_values.Rdata"))
  save(P_U_values, file = paste0(directory, "/P_U_values.Rdata"))
}

#' Load and Process Data
#'
#' This function loads matrices from saved .Rdata files and processes them for visualization.
#'
#' @param directory_res Character string indicating the directory path where the .Rdata files are located.
#' @param time_grid Numeric vector representing the time grid.
#' @return A list containing processed data for visualization.
#' @export
load_and_process_data <- function(directory_res, time_grid) {
  # Load matrices from saved .Rdata files
  load(file.path(directory_res, "X_values.Rdata"))
  load(file.path(directory_res, "U_values.Rdata"))
  load(file.path(directory_res, "P_X_values.Rdata"))
  load(file.path(directory_res, "P_U_values.Rdata"))
  
  # Process U_values for visualization
  U_value <- U_values[1, ]
  U_scaled <- sigmoid(U_value)
  U_plot <- data.frame(time = time_grid, U_scaled = U_scaled, row.names = NULL)
  
  # Process P_U_values for visualization
  P_U <- P_U_values[1, 1, ]
  
  # Process X_values for visualization
  X <- t(X_values[1:5, ])
  Xval <- data.frame(time = time_grid, X=X, row.names = NULL)
  colnames(Xval) <- c("time_grid","S","E","I","R","D")
  
  # Generate values for error area
  ymin <- mapply(function(mu, sigma) sigmoid(qnorm(0.025, mean = mu, sd = sqrt(sigma))), U_value, P_U)
  ymax <- mapply(function(mu, sigma) sigmoid(qnorm(0.975, mean = mu, sd = sqrt(sigma))), U_value, P_U)
  
  # Create data frame for P_plot
  P_plot <- data.frame(time = time_grid, ymin = ymin, ymax = ymax, row.names = NULL)
  
  # Return processed data
  return(list(U_plot = U_plot, P_plot = P_plot, U_scaled = U_scaled, Xval = Xval))
}

#' Save Processed Data
#'
#' This function saves processed data as .Rdata and .csv files.
#'
#' @param U_plot Data frame containing the processed data for plotting the contact rate.
#' @param P_plot Data frame containing the processed data for plotting the confidence interval of the contact rate.
#' @param ymin Numeric vector containing the minimum values for the confidence interval.
#' @param ymax Numeric vector containing the maximum values for the confidence interval.
#' @param U_scaled Numeric vector containing the scaled contact rate data.
#' @param Xval Data frame containing the states data to be saved as an .RData file.
#' @param directory Character string indicating the directory path where the data files will be saved.
#' @return NULL
#' @export
save_processed_data <- function(U_plot, P_plot, ymin, ymax, U_scaled, Xval, directory){
  save(U_plot, file = paste0(directory, "/U_plot.Rdata"))
  save(P_plot, file = paste0(directory, "/P_plot.Rdata"))
  save(Xval, file = paste0(directory, "/Xval.Rdata"))
  beta_with_CI <- data.frame(U_scaled, ymin, ymax)
  save(beta_with_CI, file = paste0(directory, "/yminmax.Rdata"))
  write.csv(beta_with_CI, file = paste0(directory, "/beta_with_CI.csv"), row.names = FALSE)
}

################################################################################

#' Plot Simulated Contact Rate
#'
#' This function generates a plot of the simulated contact rate
#'
#' @param real_beta_df Data frame containing the real contact rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_sim_contact_rate <- function(real_beta_df = NULL, latency_rate, recovery_rate, fatality_rate) {
  ggplot() +
    geom_line(data = real_beta_df, aes(x = time, y = real_beta, color = "Simulated contact rate"), linetype = "dashed") +
    labs(x = "Time", y = "Simulated contact rate", title = "Simulated contact rate",
         color = "Legend") +  
    scale_color_manual(values = c("Simulated contact rate" = "#E69F00"),
                       labels = c("Simulated contact rate"), name = "Lines") + 
    theme_minimal() +
    theme(legend.position = "top") +
    guides(color = guide_legend(order = 1)) +
    annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Latency Rate:", latency_rate), hjust = 0, vjust = -3, size = 3) +
    annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Recovery Rate:", recovery_rate), hjust = 0, vjust = -1, size = 3) +
    annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Fatality Rate:", fatality_rate), hjust = 0, vjust = 1, size = 3)
}

#' Plot Contact Rate: simulated data
#'
#' This function generates a plot of the estimated contact rate along with its 95%-confidence interval.
#'
#' @param U_plot Data frame containing the estimated contact rate values.
#' @param ymin Numeric vector containing the lower bounds of the confidence interval.
#' @param ymax Numeric vector containing the upper bounds of the confidence interval.
#' @param U_scaled Numeric vector containing the scaled estimated contact rate values.
#' @param real_beta_df Data frame containing the real contact rate values (optional).
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate.
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_contact_rate_sim <- function(U_plot, ymin, ymax, U_scaled, real_beta_df = NULL, latency_rate, recovery_rate, fatality_rate, lengthscale) {
    ggplot() +
      geom_ribbon(data = data.frame(time = U_plot$time, ymin = ymin, ymax = ymax), aes(x = time, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
      geom_line(data = real_beta_df, aes(x = time, y = real_beta, color = "Simulated contact rate"), linetype = "dashed") +
      geom_line(data = U_plot, aes(x = time, y = U_scaled, color = "Inferred contact rate"), linewidth = 1) +
      labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
           color = "Legend") +  
      scale_color_manual(values = c("Inferred contact rate" = "#009E73", "Simulated contact rate" = "#E69F00"),
                         labels = c("Inferred contact rate", "Simulated contact rate"), name = "Lines") + 
      scale_fill_manual(values = c("Error Area" = "#66FFCC"),
                        labels = "95%-confidence interval", name = "Ribbon") +
      theme_minimal() +
      theme(legend.position = "top") +
      guides(color = guide_legend(order = 1),
             fill = guide_legend(order = 2)) +
      annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Latency Rate:", latency_rate), hjust = 0, vjust = -3, size = 3) +
      annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Recovery Rate:", recovery_rate), hjust = 0, vjust = -1, size = 3) +
      annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Fatality Rate:", fatality_rate), hjust = 0, vjust = 1, size = 3) +
      annotate("text", x = min(U_plot$time), y = min(U_scaled), label = paste("Length Scale:", lengthscale), hjust = 0, vjust = 3, size = 3)
}

#' Plot Data: simulated data
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param obs Data frame containing the simulated compartment counts.
#' @param obs_with_noise Data frame containing the noisy compartment counts.
#' @param Xval Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the simulated, noisy and inferred compartment counts.
#' @export
plot_data_sim <- function(obs, obs_with_noise, Xval) {
  plot_title <- "Compartments counts"
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  ggplot() +
    geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 0.3) +
    geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 0.3) +
    geom_point(data = obs_with_noise, aes(x = t, y = S, color = "S", shape = "Noisy data")) +
    geom_point(data = obs_with_noise, aes(x = t, y = E, color = "E", shape = "Noisy data")) +
    geom_point(data = obs_with_noise, aes(x = t, y = I, color = "I", shape = "Noisy data")) +
    geom_point(data = obs_with_noise, aes(x = t, y = R, color = "R", shape = "Noisy data")) +
    geom_point(data = obs_with_noise, aes(x = t, y = D, color = "D", shape = "Noisy data")) +
    geom_line(data = Xval, aes(x = time_grid, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
    labs(x = "Time", y = "Count (log-scale)", title = plot_title) +
    theme_minimal() +
    scale_y_continuous(trans = 'log10') + 
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "solid", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Simulated data", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey"))) +
    scale_shape_manual(values = c("Noisy data" = 1),
                       name = "Points",
                       labels = c("Noisy data"),
                       guide = guide_legend(override.aes = list(color = "grey")))
}


################################################################################

#' Plot Contact Rate: real data
#'
#' This function generates a plot of the estimated contact rate along with its 95%-confidence interval.
#'
#' @param U_plot Data frame containing the estimated contact rate values.
#' @param ymin Numeric vector containing the lower bounds of the confidence interval.
#' @param ymax Numeric vector containing the upper bounds of the confidence interval.
#' @param U_scaled Numeric vector containing the scaled estimated contact rate values.
#' @param latency_rate Numeric value indicating the latency rate.
#' @param recovery_rate Numeric value indicating the recovery rate.
#' @param fatality_rate Numeric value indicating the fatality rate.
#' @param lengthscale Numeric value indicating the length scale.
#' @return A ggplot object displaying the contact rate plot.
#' @export
plot_contact_rate_real <- function(U_plot, ymin, ymax, U_scaled, latency_rate, recovery_rate, fatality_rate, lengthscale) {
  ggplot() +
    geom_line(data = U_plot, aes(x = time, y = U_scaled, color = "Inferred contact rate"), linewidth = 1) +
    geom_ribbon(data = data.frame(time = U_plot$time, ymin = ymin, ymax = ymax), aes(x = time, ymin = ymin, ymax = ymax, fill = "Error Area"), alpha = 0.5) +
    labs(x = "Time", y = "Contact rate", title = "Contact rate with 95%-confidence interval",
         color = "Legend") +  
    scale_color_manual(values = c("Inferred contact rate" = "darkgreen"),
                       labels = "Inferred contact rate", name = "Line") +  
    scale_fill_manual(values = c("Error Area" = "lightgreen"),
                      labels = "95%-confidence interval", name = "Ribbon") +
    theme_minimal() +
    theme(legend.position = "top")+
    guides(color = guide_legend(order = 1),
           fill = guide_legend(order = 2)) +
    annotate("text", x = max(U_plot$time), y = max(U_scaled), label = paste("Latency Rate:", latency_rate), hjust = 1, vjust = -3, size = 3) +
    annotate("text", x = max(U_plot$time), y = max(U_scaled), label = paste("Recovery Rate:", recovery_rate), hjust = 1, vjust = -1, size = 3) +
    annotate("text", x = max(U_plot$time), y = max(U_scaled), label = paste("Fatality Rate:", fatality_rate), hjust = 1, vjust = 1, size = 3) +
    annotate("text", x = max(U_plot$time), y = max(U_scaled), label = paste("Length Scale:", lengthscale), hjust = 1, vjust = 3, size = 3)
}

#' Plot Data: real data
#'
#' This function generates a plot of the observed and inferred compartment counts.
#'
#' @param obs Data frame containing the observed compartment counts.
#' @param Xval Data frame containing the inferred compartment counts.
#' @return A ggplot object displaying the observed and inferred compartment counts.
#' @export
plot_data_real <- function(obs, Xval) {
  plot_title <- "Compartments counts"
  
  cols <- c("S" = "#E69F00", 
            "E" = "#56B4E9",
            "I" = "#009E73", 
            "R" = "#F0E442",
            "D" = "#0072B2")
  
  ggplot() +
    geom_line(data = obs, aes(x = t, y = S, color = "S", linetype = "Data"), linewidth = 1) +
    geom_line(data = obs, aes(x = t, y = E, color = "E", linetype = "Data"), linewidth = 1) +
    geom_line(data = obs, aes(x = t, y = I, color = "I", linetype = "Data"), linewidth = 1) +
    geom_line(data = obs, aes(x = t, y = R, color = "R", linetype = "Data"), linewidth = 1) +
    geom_line(data = obs, aes(x = t, y = D, color = "D", linetype = "Data"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = S, color = "S", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = E, color = "E", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = I, color = "I", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = R, color = "R", linetype = "Inferred"), linewidth = 1) +
    geom_line(data = Xval, aes(x = time_grid, y = D, color = "D", linetype = "Inferred"), linewidth = 1) +
    labs(x = "Time", y = "Count (log-scale)", title = plot_title) +
    theme_minimal() +
    scale_y_continuous(trans = 'log10') + 
    scale_color_manual(values = cols,
                       name = "Compartments") +
    scale_linetype_manual(values = c("Data" = "dotted", "Inferred" = "solid"),
                          name = "Lines",
                          labels = c("Real data (simulated)", "Inferred data"),
                          guide = guide_legend(override.aes = list(color = "grey")))
}