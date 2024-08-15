#' Process Data
#'
#' This function processes the inferred results for visualization.
#'
#' @param inference_results List of the inferred results.
#' @param grids List of time grids used for inference.
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param pop Integer representing the total population.
#' @param lambda Numeric value representing the latency rate.
#' @param gamma Numeric value representing the recovery rate.
#' @param eta Numeric value representing the fatality rate.
#' @return A list containing processed data for visualization.
#' @export
process_data <- function(inference_results, grids,
                         model, pop, lambda, gamma, eta) {
  
  # Get results and time grid
  X_values <- inference_results$X_values
  U_values <- inference_results$U_values
  P_X_values <- inference_results$P_X_values
  P_U_values <- inference_results$P_U_values
  cond_numbers <- inference_results$cond_numbers
  time_grid <- grids$time_grid
  
  # Process U_values for visualization
  U_value <- U_values[1, ]
  U_scaled <- sigmoid(U_value)

  # Process P_U_values for visualization
  P_U <- P_U_values[1, 1, ]
  P_U <- pmax(P_U, 0) # Replace negative values with 0
  P_U_scaled <- sigmoid(P_U)
  
  # Process X_values for visualization
  X <- t(X_values[1:5, ])
  X_val <- data.frame(time = time_grid, X=X, row.names = NULL)
  colnames(X_val) <- c("time_grid","S","E","I","R","D")
  
  # Process P_X_values for visualization
  P_X <- P_X_values[1:5,1:5,]
  P_X <- pmax(P_X, 0) # Replace negative values with 0
  
  # Generate values for error area
  
  calculate_y_bounds <- function(X, P_X, sigma) {
    
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
  
  # Create a data frame for the plot of the transmission rate
  U_plot <- data.frame(t = time_grid, U_scaled = U_scaled, ymin = ymin_U, ymax = ymax_U, P_U_scaled = P_U_scaled, U_not_scaled = U_value, P_U_not_scaled = P_U, row.names = NULL)
  
  # Create data frame for the plot of the compartments
  X_plot <- data.frame(X_val = X_val, minmax = minmax, row.names = NULL)
  colnames(X_plot) <- c('t','S','E','I','R','D',
                        'minS','maxS',
                        'minE','maxE',
                        'minI','maxI',
                        'minR','maxR',
                        'minD','maxD')
  
  # Create a data frame for the plot of the reproduction number
  R_plot <- data.frame(t = time_grid, R = NA, ymin = NA, ymax = NA)
  for (i in 1:length(time_grid)) {
    t <- time_grid[i]
    S_t <- X_plot$S[i]
    beta_t <- U_plot$U_scaled[i]
    var_beta_t <- 
      if (model == 'SEIR') {
        R_plot$R[i] <- (beta_t / (gamma + lambda)) * (S_t / pop)
        R_plot$ymin[i] <- (U_plot$ymin[i] / (gamma + lambda)) * (X_plot$minS[i] / pop)
        R_plot$ymax[i] <- (U_plot$ymax[i] / (gamma + lambda)) * (X_plot$maxS[i] / pop)
      } else if (model == 'SEIRD') {
        R_plot$R[i] <- (beta_t / (gamma + lambda + eta)) * (S_t / pop)
        R_plot$ymin[i] <- (U_plot$ymin[i] / (gamma + lambda + eta)) * (X_plot$minS[i] / pop)
        R_plot$ymax[i] <- (U_plot$ymax[i] / (gamma + lambda + eta)) * (X_plot$maxS[i] / pop)
      }
  }
  
  # Create a data frame for the plot of the condition numbers
  cond_plot <- data.frame(time = time_grid[-1], cond_numbers = t(cond_numbers), row.names = NULL)
  colnames(cond_plot) <- c("time_grid","cond_S_obs","cond_S_ode","cond_P_X","cond_P_U")
  
  # Return processed data
  return(list(U_plot = U_plot, X_plot = X_plot, R_plot = R_plot, cond_plot = cond_plot))
}

#' Create a new directory.
#'
#' Creates a new directory with 
#' - the specified folder name, if provided
#' - the current date and time as its name if no folder name is provided
#' inside a specified base directory.
#'
#' @param base_directory The base directory where the new directory should be created.
#' @param folder_name Name for the new folder (optional).
#' @return The path of the newly created directory.
#' @export
create_directory <- function(base_directory, folder_name = NULL) {
  # Ensure the base directory ends with a slash
  if (substr(base_directory, nchar(base_directory), nchar(base_directory)) != "/") {
    base_directory <- paste0(base_directory, "/")
  }
  
  # Determine the new directory name
  if (is.null(folder_name)) {
    # Get the current date and time
    current_time <- Sys.time()
    
    # Format the current time to a string suitable for folder names
    formatted_time <- format(current_time, "%Y-%m-%d_%H-%M-%S")
    
    # Use formatted time as the folder name
    new_folder_name <- formatted_time
  } else {
    # Use the provided folder name
    new_folder_name <- folder_name
  }
  
  # Create the new directory path
  new_directory <- file.path(base_directory, new_folder_name)
  
  # Create the new directory
  dir.create(new_directory, recursive = TRUE)
  
  # Return the new directory path
  return(new_directory)
}

#' Save Processed Data and Plots
#'
#' This function saves processed data and ggplot objects as .Rdata and ggplots as .png files.
#' 
#' @param directory Character string indicating the directory path where the data files will be saved.
#' @param folder_name Character string indicating the name of the folder (optional).
#' @param inference_results List of data frames with the inferred values of X, U, P_X and P_U.
#' @param processed_data List of data frames containing the processed data for plotting the contact rate.
#' @param simulated_compartments ggplot object, plot of the simulated compartment counts.
#' @param simulated_beta ggplot object, plot of the simulated contact rate.
#' @param simulated_R ggplot object, plot of the simulated reproduction number.
#' @param real_compartments ggplot object, plot of the observed compartment counts.
#' @param compartments ggplot object, plot of the inferred compartment counts.
#' @param transmission_rate_with_CI ggplot object, plot of the inferred contact rate with 95% confidence interval.
#' @param reproduction_number_with_CI ggplot object, plot of the inferred reproduction number with approx. 95% confidence interval.
#' @param grid_plots_sep Grid of ggplot objects arranged using grid.arrange, plots of the inferred compartment counts separately.
#' @param scores_table knitr_kable of the mean scores.
#' @param plot_width Width of the saved plots.
#' @param plot_height Height of the saved plots.
#' @return NULL
#' @export
save_processed_data <- function(directory,
                                folder_name = NULL,
                                # Data
                                inference_results = NULL, 
                                processed_data = NULL, 
                                # Plots
                                simulated_compartments = NULL,
                                simulated_compartments_except_S = NULL,
                                real_compartments = NULL,
                                simulated_beta = NULL,
                                simulated_R = NULL,
                                compartments = NULL,
                                compartments_except_S = NULL,
                                transmission_rate_with_CI = NULL,
                                reproduction_number_with_CI = NULL,
                                grid_plots_sep = NULL,
                                condition_numbers = NULL,
                                scores_table = NULL,
                                # Random search
                                grid_of_mean_plots = NULL,
                                grid_of_median_plots = NULL,
                                best_params_table = NULL,
                                # Size
                                plot_width = 8,
                                plot_height = 6) {
  
  # Create directory
  if (is.null(folder_name)) {
    new_directory <- create_directory(base_directory = directory)
  } else {
    new_directory <- create_directory(base_directory = directory, 
                                      folder_name = folder_name)
  }
  
  # Get data frames to be saved, save data frames as .Rdata
  if (!is.null(inference_results)) {
    X_values <- inference_results$X_values
    U_values <- inference_results$U_values
    P_X_values <- inference_results$P_X_values
    P_U_values <- inference_results$P_U_values
    save(X_values, file = file.path(new_directory, "/X_values.Rdata"))
    save(U_values, file = file.path(new_directory, "/U_values.Rdata"))
    save(P_X_values, file = file.path(new_directory, "/P_X_values.Rdata"))
    save(P_U_values, file = file.path(new_directory, "/P_U_values.Rdata"))
  }
  if (!is.null(processed_data)) {
    U_plot <- processed_data$U_plot
    X_plot <- processed_data$X_plot
    R_plot <- processed_data$R_plot
    cond_plot <- processed_data$cond_plot
    save(U_plot, file = file.path(new_directory, "/U_plot.Rdata"))
    save(X_plot, file = file.path(new_directory, "/X_plot.Rdata"))
    save(R_plot, file = file.path(new_directory, "/R_plot.Rdata"))
    save(cond_plot, file = file.path(new_directory, "/cond_plot.Rdata"))
  }

  # Save ggplot images
  if (!is.null(simulated_compartments)) {
    ggsave(filename = file.path(new_directory, "simulated_compartments.png"), 
           plot = simulated_compartments, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(simulated_compartments_except_S)) {
    ggsave(filename = file.path(new_directory, "simulated_compartments_except_S.png"), 
           plot = simulated_compartments_except_S, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(real_compartments)) {
    ggsave(filename = file.path(new_directory, "real_compartments.png"), 
           plot = real_compartments, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(simulated_beta)) {
    ggsave(filename = file.path(new_directory, "simulated_beta.png"), 
           plot = simulated_beta, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(simulated_R)) {
    ggsave(filename = file.path(new_directory, "simulated_R.png"), 
           plot = simulated_R, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(real_compartments)) {
    ggsave(filename = file.path(new_directory, "real_compartments.png"), 
           plot = real_compartments, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(compartments)) {
    ggsave(filename = file.path(new_directory, "compartments.png"), 
           plot = compartments, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(compartments_except_S)) {
    ggsave(filename = file.path(new_directory, "compartments_except_S.png"), 
           plot = compartments_except_S, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(transmission_rate_with_CI)) {
    ggsave(filename = file.path(new_directory, "transmission_rate_with_CI.png"), 
           plot = transmission_rate_with_CI, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(reproduction_number_with_CI)) {
    ggsave(filename = file.path(new_directory, "reproduction_number_with_CI.png"), 
           plot = reproduction_number_with_CI, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(grid_plots_sep)) {
    ggsave(filename = file.path(new_directory, "grid_plots_sep.png"),
           plot = grid_plots_sep,
           bg = "white",
           width = plot_width,
           height = plot_height)
  }
  if (!is.null(condition_numbers)) {
    ggsave(filename = file.path(new_directory, "condition_numbers.png"),
           plot = condition_numbers,
           bg = "white",
           width = plot_width,
           height = plot_height)
  }
  if (!is.null(scores_table)) {
    file_path_html <- file.path(new_directory, "scoring_results.html")
    save_kable(scores_table, file = file_path_html, type = "html")
  }
  if (!is.null(grid_of_mean_plots)) {
    ggsave(filename = file.path(new_directory, "grid_of_mean_plots.png"),
           plot = grid_of_mean_plots,
           bg = "white",
           width = plot_width,
           height = plot_height)
  }
  if (!is.null(grid_of_median_plots)) {
    ggsave(filename = file.path(new_directory, "grid_of_median_plots.png"),
           plot = grid_of_median_plots,
           bg = "white",
           width = plot_width,
           height = plot_height)
  }
  if (!is.null(best_params_table)) {
    file_path_html <- file.path(new_directory, "best_params_table.html")
    save_kable(best_params_table, file = file_path_html, type = "html")
  }
}

#' Load real data
#'
#' Load real data based on the specified region and frequency.
#'
#' @param region Character, region for which to load the data ('BE' or 'GE').
#' @param daily_or_weekly Character, frequency of the data ('daily' or 'weekly').
#' @param augmented Boolean, TRUE if augmented data are wished.
#' @param directory_data Character, directory path where the data is stored.
#'
#' @return A list containing the loaded data including observations, population, and real beta (if available).
#'
#' @examples
#' load_data("GE", "weekly", "/path/to/data/directory/real")
#'
#' @export
load_data <- function(region, daily_or_weekly, augmented = FALSE, directory_data) {
  # Check if region is valid
  if (!(region %in% c('BE', 'GE', ''))) {
    stop("Error: Invalid region. Region must be 'BE' or 'GE'.")
  }
  
  # Check if daily_or_weekly is valid
  if (!(daily_or_weekly %in% c('daily', 'weekly', ''))) {
    stop("Error: Invalid value for daily_or_weekly. Must be 'daily' or 'weekly'.")
  }
  
  if (augmented == FALSE) {
    region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory_data, region_filename))
    obs <- observations
    population <- readRDS(file.path(directory_data, population_filename))
    res <- list(obs = obs, population = population)
  }
  else if (augmented == TRUE) {
    region_filename <- paste0("real_data_augmented_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_augmented_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory_data, region_filename))
    obs <- compartments
    population <- readRDS(file.path(directory_data, population_filename))
    res <- list(obs = obs, population = population)
  }
  
  # Return loaded data
  return(res)
}
