#####################################
############# WORKFLOW ##############
#####################################

# Import functions
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/functions_for_inference.R')
source('~/Documents/GitHub/proboder/saving_loading.R')
source('~/Documents/GitHub/proboder/scoring.R')
source('~/Documents/GitHub/proboder/plotting.R')

# Necessary packages
library(Matrix) # for sparseMatrix() and expm()
library(numDeriv) # for jacobian()
library(matrixcalc) # for svd.inverse()
library(ggplot2) # for ggplot()
library(gridExtra) # for multiple plots

# Choose data to be imported
type <- 'simulated_LSODA' # set 'real' for real data, 'simulated_LSODA' for simulated data using LSODA, and 'simulated_HETTMO' for simulated data using HETTMO
region <- 'BE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly' (if 'real' data selected)
if(type == 'simulated_LSODA'){
  directory_data <- "~/Documents/GitHub/proboder/Data/LSODA" # directory of data
}else if(type == 'simulated_HETTMO'){
  directory_data <- "~/Documents/GitHub/proboder/Data/HETTMO" # directory of data
}else if(type == 'real'){
  directory_data <- "~/Documents/GitHub/proboder/Data/real" # directory of data
}

# Import data
data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- data$obs
obs_with_noise <- data$obs_with_noise
params <- data$params
if(type != 'real'){
  real_beta <- data$real_beta
}

# Sanity check.
head(obs)

#####################################
########## INITIALIZATION ###########
#####################################

# If using data simulated with LSODA:
# lambda = 0.6, gamma = 0.4, eta = 0.2, noise_obs = 0.1
# beta0 = 0.5, beta0prime = 0.3
# pop = 1000

initial_params <- 
  initialization(obs_with_noise, beta0 = 0.1, beta0prime = 0.3, 
                 lambda = 0.6, gamma = 0.4, eta = 0.2, 
                 l = 10, scale = 1, noise_obs = 1,
                 noise_X = 0.001, noise_U = 1,
                 noise_wiener_X = 1000, noise_wiener_U = 0.01,
                 pop = 1000)

# If using data simulated with HETTMO:
# lambda = 0.3703704, gamma = 0.3703704, eta = 0
# pop = 1e+05

# initial_params <- 
#   initialization(obs, beta0 = 0.9721224, beta0prime = -1.5, 
#                  lambda = 0.3703704, gamma = 0.3703704, eta = 0, 
#                  l = 1.2, scale = 3, noise_obs = 100,
#                  noise_wiener_X = 1e+03, noise_wiener_U = 0.1,
#                  pop = 1e+05)

#####################################
############# ALGORITHM #############
#####################################

# Data grid
data_grid <- obs[,'t']

# ODE grid
ode_grid <- data_grid # no more points than observations

# Adding more points than observations
#num_points_between <- 2
#for (i in 1:(length(data_grid) - 1)) {
  # Generate equidistant points between the current and next data point
  #equidistant_points <- seq(data_grid[i], data_grid[i + 1], length.out = num_points_between + 2)[-c(1, num_points_between + 2)]
  # Append the equidistant points to the ODE grid
  #ode_grid <- c(ode_grid, equidistant_points)
#}

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
directory_res = "~/Documents/GitHub/proboder/Results"
# Save results to the specified directory
save_results_as_Rdata(X_values, U_values, P_X_values, P_U_values, directory_res)

###############################################
########### SCORING & VISUALIZATION ###########
###############################################

# ---------------------
# Extract relevant data
# ---------------------

# Load and process data from the specified directory
processed_data <- load_and_process_data(directory_res,time_grid)
U_plot <- processed_data$U_plot
X_plot <- processed_data$X_plot

# Save processed data to the specified directory
save_processed_data(U_plot, X_plot, directory_res)

# Create data frame for real beta values (if available)
if(type!='real'){
  real_beta_df <- data.frame(time = data_grid, real_beta = real_beta)
  colnames(real_beta_df) <- c('t','beta')
}

# --------
# Scoring
# --------

squared_prediction_error(U_plot,real_beta_df)
negative_log_predictive_density(U_plot,real_beta_df)
continuous_ranked_probability_score(U_plot,real_beta_df)

# --------
# Plotting
# --------

setwd(directory_res)

# Plot compartment counts
pdf("SEIRD-counts.pdf", width = 8, height = 6)
plot_data_sim(obs,obs_with_noise,X_plot)
dev.off()
plot_data_sim(obs,obs_with_noise,X_plot)

# Plot compartment counts separately
plots <- plot_compartment(obs,obs_with_noise,X_plot)
for (i in 1:5) {
  pdf(paste0("SEIRD-counts-sep-with-CI-", i, ".pdf"), width = 8, height = 6)
  plot <- plots[[i]]
  print(plot)
  dev.off()
}
for (i in 1:5) {
  plot <- plots[[i]]
  print(plot)
}

# Get some fixed values to plot together with contact rate
lambda <- initial_params$lambda
gamma <- initial_params$gamma
eta <- initial_params$eta
l <- initial_params$l

# Plot simulated contact rate
pdf("sim-contact-rate.pdf", width = 8, height = 6)
plot_sim_contact_rate(real_beta_df, lambda, gamma, eta)
dev.off()
plot_sim_contact_rate(real_beta_df, lambda, gamma, eta)

# Plot inferred contact rate
pdf("inf-contact-rate-with-CI.pdf", width = 8, height = 6)
plot_contact_rate_with_CI(U_plot, real_beta_df, lambda, gamma, eta, l)
dev.off()
plot_contact_rate_with_CI(U_plot, real_beta_df, lambda, gamma, eta, l)