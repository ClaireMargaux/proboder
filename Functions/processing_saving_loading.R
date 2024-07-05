#' Process Data
#'
#' This function processes the inferred results for visualization.
#'
#' @param inference_results List of the inferred results.
#' @param grids List of time grids used for inference.
#' @return A list containing processed data for visualization.
#' @export
process_data <- function(inference_results, grids) {
  # Get results and time grid
  X_values <- inference_results$X_values
  U_values <- inference_results$U_values
  P_X_values <- inference_results$P_X_values
  P_U_values <- inference_results$P_U_values
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

#' Create a timestamped directory
#'
#' Creates a new directory with the current date and time as its name inside a specified base directory.
#'
#' @param base_directory The base directory where the new directory should be created.
#' @return The path of the newly created directory.
#' @export
#' @examples
#' \dontrun{
#' base_directory <- '~/Documents/GitHub/proboder/Results'
#' new_directory <- create_timestamped_directory(base_directory)
#' print(new_directory)
#' }
create_timestamped_directory <- function(base_directory) {
  # Ensure the base directory ends with a slash
  if (substr(base_directory, nchar(base_directory), nchar(base_directory)) != "/") {
    base_directory <- paste0(base_directory, "/")
  }
  
  # Get the current date and time
  current_time <- now()
  
  # Format the current time to a string suitable for folder names
  formatted_time <- format(current_time, "%Y-%m-%d_%H-%M-%S")
  
  # Create the new directory path
  new_directory <- file.path(base_directory, formatted_time)
  
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
#' @param inference_results List of data frames with the inferred values of X, U, P_X and P_U.
#' @param processed_data List of data frames containing the processed data for plotting the contact rate.
#' @param simulated_compartments ggplot object, plot of the simulated compartment counts.
#' @param simulated_beta ggplot object, plot of the simulated contact rate.
#' @param compartments ggplot object, plot of the inferred compartment counts.
#' @param contact_rate_with_CI ggplot object, plot of the inferred contact rate with 95% confidence interval.
#' @param grid_plots_sep Grid of ggplot objects arranged using grid.arrange, plots of the inferred compartment counts separately.
#' @param styled_table knitr_kable of the mean scores.
#' @param plot_width Width of the saved plots.
#' @param plot_height Height of the saved plots.
#' @return NULL
#' @export
save_processed_data <- function(directory,
                                # Data
                                inference_results, 
                                processed_data, 
                                # Plots
                                simulated_compartments = NULL,
                                simulated_beta = NULL,
                                compartments = NULL,
                                contact_rate_with_CI = NULL,
                                grid_plots_sep = NULL,
                                styled_table = NULL,
                                plot_width = 8,
                                plot_height = 6) {
  
  # Create directory
  new_directory <- create_timestamped_directory(directory)
  
  # Get data frames to be saved
  X_values <- inference_results$X_values
  U_values <- inference_results$U_values
  P_X_values <- inference_results$P_X_values
  P_U_values <- inference_results$P_U_values
  U_plot <- processed_data$U_plot
  X_plot <- processed_data$X_plot
  
  # Save data frames as .Rdata
  save(U_plot, file = file.path(new_directory, "/U_plot.Rdata"))
  save(X_plot, file = file.path(new_directory, "/X_plot.Rdata"))
  save(X_values, file = file.path(new_directory, "/X_values.Rdata"))
  save(U_values, file = file.path(new_directory, "/U_values.Rdata"))
  save(P_X_values, file = file.path(new_directory, "/P_X_values.Rdata"))
  save(P_U_values, file = file.path(new_directory, "/P_U_values.Rdata"))

  # Save ggplot images
  if (!is.null(simulated_compartments)) {
    ggsave(filename = file.path(new_directory, "simulated_compartments.png"), 
           plot = simulated_compartments, 
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
  if (!is.null(compartments)) {
    ggsave(filename = file.path(new_directory, "compartments.png"), 
           plot = compartments, 
           bg = "white",
           width = plot_width, 
           height = plot_height)
  }
  if (!is.null(contact_rate_with_CI)) {
    ggsave(filename = file.path(new_directory, "contact_rate_with_CI.png"), 
           plot = contact_rate_with_CI, 
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
  if (!is.null(styled_table)) {
    file_path_html <- file.path(new_directory, "scoring_results.html")
    save_kable(styled_table, file = file_path_html, type = "html")
  }
}

#' Load real data
#'
#' Load real data based on the specified region and frequency.
#'
#' @param region Character, region for which to load the data ('BE' or 'GE').
#' @param daily_or_weekly Character, frequency of the data ('daily' or 'weekly').
#' @param augmented Boolean, TRUE if augmented data are wished.
#' @param directory Character, directory path where the data is stored.
#'
#' @return A list containing the loaded data including observations, population, and real beta (if available).
#'
#' @examples
#' load_data("GE", "weekly", "/path/to/data/directory/real")
#'
#' @export
load_data <- function(region, daily_or_weekly, augmented = FALSE, directory) {
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
    load(file.path(directory, region_filename))
    obs <- observations
    population <- readRDS(file.path(directory, population_filename))
    return <- list(obs = obs, population = population)
  }
  else if (augmented == TRUE) {
    region_filename <- paste0("real_data_augmented_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_augmented_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory, region_filename))
    obs <- compartments
    population <- readRDS(file.path(directory, population_filename))
    return <- list(obs = obs, population = population)
  }
  
  # Return loaded data
  return(return)
}