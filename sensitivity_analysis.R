#################################################
############## Sensitivity analysis #############
#################################################

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
library(sensitivity) # for the sensitivity analysis

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

# Data grid
data_grid <- obs[,'t']

# ODE grid
ode_grid <- data_grid # no more points than observations

# Overall time grid
time_grid <- sort(unique(c(data_grid, ode_grid)))

# Time steps
steps <- ode_grid[2]-ode_grid[1]

param_grid <- expand.grid(
  l = seq(1, 10, by = 0.05),
  noise_wiener_X = seq(20, 1000, by = 25),
  noise_wiener_U = seq(0.005, 1, by = 0.005),
  beta0prime = seq(0, 1, by = 0.2)
)

set.seed(123) # for reproducibility
num_param_sets <- 1000
sampled_grid <- param_grid %>% sample_n(num_param_sets)

X <- as.matrix(sampled_grid)
rownames(X) <- NULL
X1 <- data.frame(X[1:(num_param_sets/2), ])
X2 <- data.frame(X[(num_param_sets/2+1):num_param_sets, ])

sobolfun <- function(X) {
  # x is a matrix where each row is a set of parameters
  
  CRPS <- numeric()
  
  for(i in 1:nrow(X)) {
    
  x <- X[i,]
    
  l <- x[,1]
  noise_wiener_X <- x[,2]
  noise_wiener_U <- x[,3]
  beta0prime <- x[,4]
  
  # Initialize with current set of parameters
  initial_params <- initialization(
    obs_with_noise, beta0 = real_beta[1], beta0prime = beta0prime,
    lambda = params$lambda, gamma = params$gamma, eta = params$eta,
    l = l, scale = 1, noise_obs = params$obs_noise,
    noise_X = sqrt(params$obs_noise), noise_U = 0.01,
    noise_wiener_X = noise_wiener_X, noise_wiener_U = noise_wiener_U,
    pop = params$pop, num_points_between = 0, 
    num_initial_values = 1
  )
  
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
  
  # Compute the mean CRPS
  score_CRPS <- mean(continuous_ranked_probability_score(U_plot,real_beta_df))
  
  CRPS <- rbind(CRPS,score_CRPS)
  
  }
  
  rownames(CRPS) <- NULL
  
  return(CRPS)
}

# Start counting computation time
tic("Duration of the whole workflow")

# Sobol design
sobol_design <- sobol2007(model = NULL, X1 = X1, X2 = X2, nboot = 100)
y <- sobolfun(sobol_design$X) 
sobol_design <- tell(x=sobol_design, y)

# Total computation time
toc()

print(sobol_design)
ggplot(sobol_design)
