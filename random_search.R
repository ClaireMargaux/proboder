##########################################
############## Random search #############
##########################################

# Import functions
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/inference.R')
source('~/Documents/GitHub/proboder/saving_loading.R')
source('~/Documents/GitHub/proboder/scoring.R')
source('~/Documents/GitHub/proboder/plotting.R')

# Necessary packages
library(Matrix) # for sparseMatrix() and expm()
library(numDeriv) # for jacobian()
library(matrixcalc) # for svd.inverse()
library(knitr) # for nice tables
library(kableExtra) # for storing nice tables
library(ggplot2) # for ggplot()
library(gridExtra) # for multiple plots
library(tictoc) # for benchmarking
library(dplyr) # for treating data
library(progress) # to track progress of grid search

# ----
# Data
# ----

# Choose data to be imported
type <- 'simulated_LSODA_sin' # set 'simulated_LSODA' for simulated data using LSODA, and 'simulated_HETTMO' for simulated data using HETTMO
region <- ''
daily_or_weekly <- ''

if(type == 'simulated_LSODA_sin'){
  directory_data <- "~/Documents/GitHub/proboder/Data/LSODA/sin" # directory of data
}else if(type == 'simulated_LSODA_log'){
  directory_data <- "~/Documents/GitHub/proboder/Data/LSODA/log" # directory of data
}else if(type == 'simulated_HETTMO'){
  directory_data <- "~/Documents/GitHub/proboder/Data/HETTMO" # directory of data
}

# Import data
data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- as.data.frame(data$obs)
obs_with_noise <- data$obs_with_noise # only with LSODA
params <- data$params
real_beta <- data$real_beta

# Sanity check.
head(obs_with_noise)

# Create data frame for real beta values
real_beta_df <- data.frame(time = obs$t, real_beta = real_beta)
colnames(real_beta_df) <- c('t','beta')

#------------
# Grid search
# -----------

param_grid <- expand.grid(
  #l = seq(1, 10, by = 0.05),
  l = 8.65,
  noise_wiener_X = seq(20, 1000, by = 25),
  #noise_wiener_X = 345,
  #noise_wiener_U = seq(0.005, 1, by = 0.005),
  noise_wiener_U = 0.015,
  #beta0prime = c(0, 0.1, 0.3, 0.5, 1)
  beta0prime = 0.1
)

set.seed(123) # for reproducibility

# Random search grid
num_param_sets <- 40
sampled_grid <- param_grid %>% sample_n(num_param_sets)

num_params <- dim(param_grid)[2]
num_predictions <- length(obs_with_noise$t)

# Matrix for results
results_grid_search <- matrix(
 nrow = num_param_sets,
 ncol = (num_params + 3*num_predictions)
)

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = nrow(sampled_grid), clear = FALSE, width = 60
)

# Start counting computation time
tic("Duration of the whole workflow")

for (i in 1:num_param_sets) {
  param_i <- sampled_grid[i, ]
  
  # Initialize with current set of parameters
  initial_params <- initialization(
    obs_with_noise, beta0 = real_beta[1], beta0prime = param_i$beta0prime,
    lambda = params$lambda, gamma = params$gamma, eta = params$eta,
    l = param_i$l, scale = 1, noise_obs = params$obs_noise,
    noise_X = sqrt(params$obs_noise), noise_U = 0.01,
    noise_wiener_X = param_i$noise_wiener_X, noise_wiener_U = param_i$noise_wiener_U,
    pop = params$pop, num_points_between = 0,
    num_initial_values = 1
  )
  
  # ---------
  # Inference
  # ---------
  
  # Data grid
  data_grid <- obs[,'t']
  
  # ODE grid
  ode_grid <- data_grid # no more points than observations
  
  # Adding more points than observations
  for (j in 1:(length(data_grid) - 1)) {
    # Generate equidistant points between the current and next data point
    equidistant_points <- seq(data_grid[j], data_grid[j + 1], length.out = initial_params$num_points_between + 2)[-c(1, initial_params$num_points_between + 2)]
    # Append the equidistant points to the ODE grid
    ode_grid <- c(ode_grid, equidistant_points)
    ode_grid <- sort(ode_grid)
  }
  
  # Overall time grid
  time_grid <- sort(unique(c(data_grid, ode_grid)))
  
  # Time steps
  steps <- ode_grid[2]-ode_grid[1]
  
  # Run inference
  inference_results <- inference(time_grid, data_grid, ode_grid, steps, obs, initial_params)
  
  X_values <- inference_results$X_values
  U_values <- inference_results$U_values
  P_X_values <- inference_results$P_X_values
  P_U_values <- inference_results$P_U_values
  
  # ------------
  # Save results
  # ------------
  
  # Specify directory for results
  directory_res <- "~/Documents/GitHub/proboder/Results/temp" # directory of data
  
  # Save results to the specified directory
  save_results_as_Rdata(X_values, U_values, P_X_values, P_U_values, directory_res)
  
  # ---------------------
  # Extract relevant data
  # ---------------------
  
  # Load and process data from the specified directory
  processed_data <- load_and_process_data(directory_res,time_grid)
  U_plot <- processed_data$U_plot
  X_plot <- processed_data$X_plot
  
  # Save processed data to the specified directory
  save_processed_data(U_plot, X_plot, directory_res)
  
  # --------
  # Scoring
  # --------
  
  # Compute the different scores
  score_SPE <- squared_prediction_error(U_plot,real_beta_df)
  score_NLPD <-negative_log_predictive_density(U_plot,real_beta_df)
  score_CRPS <- continuous_ranked_probability_score(U_plot,real_beta_df)
  
  # Store results (during grid search)
  results_grid_search[i, ] <- c(param_i$l, param_i$noise_wiener_X, param_i$noise_wiener_U, param_i$beta0prime, score_SPE, score_NLPD, score_CRPS)
  
  rm(score_SPE,score_NLPD,score_CRPS)
  
  # Update progress bar
  pb$tick()
  
} # End of grid search

# Total computation time
toc()

# ---------------
# Extract results
# ---------------

# Extract the parameter values
params <- results_grid_search[, 1:num_params]

# Extract the score vectors
SPE_scores <- results_grid_search[, (num_params + 1):(num_params + num_predictions)]
NLPD_scores <- results_grid_search[, (num_params + num_predictions + 1):(num_params + 2 * num_predictions)]
CRPS_scores <- results_grid_search[, (num_params + 2 * num_predictions + 1):(num_params + 3 * num_predictions)]

# Calculate the mean and median for each score
mean_SPE <- rowMeans(SPE_scores)
median_SPE <- apply(SPE_scores, 1, median)

mean_NLPD <- rowMeans(NLPD_scores)
median_NLPD <- apply(NLPD_scores, 1, median)

mean_CRPS <- rowMeans(CRPS_scores)
median_CRPS <- apply(CRPS_scores, 1, median)

# Combine into a new data frame
results_summary <- data.frame(
  l = params[, 1],
  noise_wiener_X = params[, 2],
  noise_wiener_U = params[, 3],
  beta0prime = params[, 4],
  mean_SPE = mean_SPE,
  median_SPE = median_SPE,
  mean_NLPD = mean_NLPD,
  median_NLPD = median_NLPD,
  mean_CRPS = mean_CRPS,
  median_CRPS = median_CRPS
)

# View the results summary
print(results_summary)

# --------------------
# Sensitivity analysis
# --------------------

# Function to plot mean scores for a given parameter
plot_mean_scores <- function(df, parameter) {
  df %>%
    ggplot(aes_string(x = parameter)) +
    geom_point(aes(y = mean_SPE, color = "SPE")) +
    geom_point(aes(y = mean_NLPD, color = "NLPD")) +
    geom_point(aes(y = mean_CRPS, color = "CRPS")) +
    labs(title = paste("Mean scores for parameter", parameter), x = parameter, y = "Mean score") +
    scale_color_manual(values = c("SPE" = "#E69F00", "NLPD" = "#56B4E9", "CRPS" = "#009E73")) +
    theme_minimal()
}

plot_median_scores <- function(df, parameter) {
  df %>%
    ggplot(aes_string(x = parameter)) +
    geom_point(aes(y = median_SPE, color = "SPE")) +
    geom_point(aes(y = median_NLPD, color = "NLPD")) +
    geom_point(aes(y = median_CRPS, color = "CRPS")) +
    labs(title = paste("Median scores for parameter", parameter), x = parameter, y = "Median score") +
    scale_color_manual(values = c("SPE" = "#E69F00", "NLPD" = "#56B4E9", "CRPS" = "#009E73")) +
    theme_minimal()
}

# Identify parameters (columns) to plot
parameters_to_plot <- c('l','noise_wiener_X','noise_wiener_U','beta0prime')

# Create a list to store the plots
plot_list_mean <- list()
plot_list_median <- list()

# Generate and store plots for each parameter
for (param in parameters_to_plot) {
  plot_list_mean[[param]] <- plot_mean_scores(results_summary, param)
  plot_list_median[[param]] <- plot_median_scores(results_summary, param)
}

# Print all plots
for (plot in plot_list_mean) {
  print(plot)
}
for (plot in plot_list_median) {
  print(plot)
}

# --------------------------------
# Find best set of hyperparameters
# --------------------------------

# Find the best set of hyperparameters for each score
best_spe <- results_summary[which.min(results_summary$mean_SPE), ]
best_nlpd <- results_summary[which.min(results_summary$mean_NLPD), ]
best_crps <- results_summary[which.min(results_summary$mean_CRPS), ]

# Check if there is a set minimizing all three scores
best_all <- results_summary %>%
  filter(mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS))

if (nrow(best_all) > 0) {
  best_params <- best_all[1, ]
} else {
  # Find sets minimizing at least two of the scores
  best_two <- results_summary %>%
    filter((mean_SPE == min(mean_SPE) & mean_NLPD == min(mean_NLPD)) |
             (mean_SPE == min(mean_SPE) & mean_CRPS == min(mean_CRPS)) |
             (mean_NLPD == min(mean_NLPD) & mean_CRPS == min(mean_CRPS)))
  
  if (nrow(best_two) > 0) {
    best_params <- best_two[1, ]
  } else {
    # If no sets minimize at least two scores, prioritize minimal CRPS
    best_params <- best_crps
  }
}

print(best_params)