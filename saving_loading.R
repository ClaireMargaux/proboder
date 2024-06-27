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
  if (!(type %in% c('simulated_LSODA_sin', 'simulated_LSODA_log', 'simulated_HETTMO', 'real'))) {
    stop("Error: Invalid type.")
  }
  
  # Check if region is valid
  if (!(region %in% c('BE', 'GE', ''))) {
    stop("Error: Invalid region. Region must be 'BE' or 'GE'.")
  }
  
  # Check if daily_or_weekly is valid
  if (!(daily_or_weekly %in% c('daily', 'weekly', ''))) {
    stop("Error: Invalid value for daily_or_weekly. Must be 'daily' or 'weekly'.")
  }
  
  # Load data based on type
  if (type == 'simulated_LSODA_sin') {
    load(file.path(directory, "simulated_data_LSODA.Rdata"))
    obs <- df
    load(file.path(directory, "simulated_noisy_data_LSODA.Rdata"))
    obs_with_noise <- df_with_noise
    load(file.path(directory, "simulated_params_LSODA.Rdata"))
    params <- df_params
    real_beta <- readRDS(file.path(directory, "simulated_beta_LSODA.Rds"))
    return <- list(obs = obs, obs_with_noise = obs_with_noise, params = params, real_beta = real_beta)
  }
  if (type == 'simulated_LSODA_log') {
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
    return <- list(obs = obs, population = population)
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

  # Process P_U_values for visualization
  P_U <- P_U_values[1, 1, ]
  P_U_scaled <- sigmoid(P_U)
  
  # Process X_values for visualization
  X <- t(X_values[1:5, ])
  X_val <- data.frame(time = time_grid, X=X, row.names = NULL)
  colnames(X_val) <- c("time_grid","S","E","I","R","D")
  
  # Process P_X_values for visualization
  P_X <- P_X_values[1:5,1:5,]
  
  # Generate values for error area
  
  calculate_y_bounds <- function(X, P_X, sigma) {
    # Replace negative values in P_X with 0
    P_X[P_X < 0] <- 0
    
    # Calculate ymin and ymax
    ymin <- mapply(function(mu, sigma) (qnorm(0.025, mean = mu, sd = sqrt(sigma))), X, P_X)
    ymax <- mapply(function(mu, sigma) (qnorm(0.975, mean = mu, sd = sqrt(sigma))), X, P_X)
    
    return(list(ymin = ymin, ymax = ymax))
  }
  
  ymin_S <- calculate_y_bounds(X_val$S,P_X[1,1,])[[1]]
  ymax_S <- calculate_y_bounds(X_val$S,P_X[1,1,])[[2]]
  
  ymin_E <- calculate_y_bounds(X_val$E,P_X[2,2,])[[1]]
  ymax_E <- calculate_y_bounds(X_val$E,P_X[2,2,])[[2]]
  
  ymin_I <- calculate_y_bounds(X_val$I,P_X[3,3,])[[1]]
  ymax_I <- calculate_y_bounds(X_val$I,P_X[3,3,])[[2]]
  
  ymin_R <- calculate_y_bounds(X_val$R,P_X[4,4,])[[1]]
  ymax_R <- calculate_y_bounds(X_val$R,P_X[4,4,])[[2]]
  
  ymin_D <- calculate_y_bounds(X_val$D,P_X[5,5,])[[1]]
  ymax_D <- calculate_y_bounds(X_val$D,P_X[5,5,])[[2]]
  
  minmax <- data.frame(
    ymin_S = ymin_S,
    ymax_S = ymax_S, 
    ymin_E = ymin_E,
    ymax_E = ymax_E,
    ymin_I = ymin_I,
    ymax_I = ymax_I,
    ymin_R = ymin_R,
    ymax_R = ymax_R,
    ymin_D = ymin_D,
    ymax_D = ymax_D
  )
  
  ymin_U <- mapply(function(mu, sigma) sigmoid(qnorm(0.025, mean = mu, sd = sqrt(sigma))), U_value, P_U)
  ymax_U <- mapply(function(mu, sigma) sigmoid(qnorm(0.975, mean = mu, sd = sqrt(sigma))), U_value, P_U)
  
  # Create a data frame for the plot of the contact rate
  U_plot <- data.frame(t = time_grid, U_scaled = U_scaled, ymin = ymin_U, ymax = ymax_U, P_U_scaled = P_U_scaled, row.names = NULL)
  
  # Create data frame for the plot of the compartments
  X_plot <- data.frame(X_val = X_val, minmax = minmax, row.names = NULL)
  colnames(X_plot) <- c('t','S','E','I','R','D',
                        'minS','maxS',
                        'minE','maxE',
                        'minI','maxI',
                        'minR','maxR',
                        'minD','maxD')
  
  # Return processed data
  return(list(U_plot = U_plot, X_plot = X_plot))
}

#' Save Processed Data
#'
#' This function saves processed data as .Rdata and .csv files.
#'
#' @param U_plot Data frame containing the processed data for plotting the contact rate.
#' @param X_plot Data frame containing the processed data for plotting the compartment counts.
#' @param directory Character string indicating the directory path where the data files will be saved.
#' @return NULL
#' @export
save_processed_data <- function(U_plot, X_plot, directory){
  save(U_plot, file = paste0(directory, "/U_plot.Rdata"))
  save(X_plot, file = paste0(directory, "/X_plot.Rdata"))
 }