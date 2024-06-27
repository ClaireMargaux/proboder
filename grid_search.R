########################################
############## Grid search #############
########################################

# Import functions
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/functions_for_inference.R')
source('~/Documents/GitHub/proboder/saving_loading.R')
source('~/Documents/GitHub/proboder/scoring.R')
source('~/Documents/GitHub/proboder/plotting_sim.R')

# Necessary packages
library(Matrix) # for sparseMatrix() and expm()
library(numDeriv) # for jacobian()
library(matrixcalc) # for svd.inverse()
library(knitr) # for nice tables
library(kableExtra) # for storing nice tables
library(ggplot2) # for ggplot()
library(gridExtra) # for multiple plots
library(tictoc) # for benchmarking
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
head(obs)

# Create data frame for real beta values
real_beta_df <- data.frame(time = obs$t, real_beta = real_beta)
colnames(real_beta_df) <- c('t','beta')

#------------
# Grid search
# -----------

param_grid <- expand.grid(
  l = seq(0.5, 15, by = 0.5),
  noise_wiener_X = c(1, 10, 100),
  noise_wiener_U = c(0.001, 0.01, 0.1),
  beta0prime = c(-1, -0.1, 0, 0.1, 1),
  num_points_between = c(0,1)
)

set.seed(123) # for reproducibility
sampled_grid <- param_grid %>% sample_n(100)

results_grid_search <- data.frame()

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = nrow(sampled_grid), clear = FALSE, width = 60
)

# Start counting computation time
tic("Duration of the whole workflow")

for (i in 1:nrow(sampled_grid)) {
  param_i <- sampled_grid[i, ]
  
  # Initialize with current set of parameters
  initial_params <- initialization(
    obs_with_noise, beta0 = real_beta[1], beta0prime = param_i$beta0prime,
    lambda = params$lambda, gamma = params$gamma, eta = params$eta,
    l = param_i$l, scale = 1, noise_obs = params$obs_noise,
    noise_X = sqrt(params$obs_noise), noise_U = 0.01,
    noise_wiener_X = param_i$noise_wiener_X, noise_wiener_U = param_i$noise_wiener_U,
    #noise_wiener_X = (params$obs_noise), noise_wiener_U = 0.01,
    pop = params$pop, num_points_between = param_i$num_points_between
  )
  
  # ---------
  # Inference
  # ---------
  
  # Data grid
  data_grid <- obs[,'t']
  
  # ODE grid
  ode_grid <- data_grid # no more points than observations
  
  # Adding more points than observations
  for (i in 1:(length(data_grid) - 1)) {
    # Generate equidistant points between the current and next data point
    equidistant_points <- seq(data_grid[i], data_grid[i + 1], length.out = initial_params$num_points_between + 2)[-c(1, initial_params$num_points_between + 2)]
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
  SPE <- squared_prediction_error(U_plot,real_beta_df)
  NLPD <-negative_log_predictive_density(U_plot,real_beta_df)
  CRPS <- continuous_ranked_probability_score(U_plot,real_beta_df)
  
  # Store results (during grid search)
  results_grid_search <- rbind(results_grid_search, cbind(param_i, SPE, NLPD, CRPS))
  
  # Update progress bar
  pb$tick()
  
} # End of grid search

# Total computation time
toc()

# --------------------------------
# Find best set of hyperparameters
# --------------------------------

# Find the best set of hyperparameters for each score
best_spe <- results_grid_search[which.min(results_grid_search$SPE), ]
best_nlpd <- results_grid_search[which.min(results_grid_search$NLPD), ]
best_crps <- results_grid_search[which.min(results_grid_search$CRPS), ]

# Check if there is a set minimizing all three scores
best_all <- results_grid_search %>%
  filter(SPE == min(SPE) & NLPD == min(NLPD) & CRPS == min(CRPS))

if (nrow(best_all) > 0) {
  best_params <- best_all[1, ]
} else {
  # Find sets minimizing at least two of the scores
  best_two <- results_grid_search %>%
    filter((SPE == min(SPE) & NLPD == min(NLPD)) |
             (SPE == min(SPE) & CRPS == min(CRPS)) |
             (NLPD == min(NLPD) & CRPS == min(CRPS)))
  
  if (nrow(best_two) > 0) {
    best_params <- best_two[1, ]
  } else {
    # If no sets minimize at least two scores, prioritize minimal CRPS
    best_params <- best_crps
  }
}

print(best_params)