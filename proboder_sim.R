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
library(knitr) # for nice tables
library(kableExtra) # for storing nice tables
library(ggplot2) # for ggplot()
library(gridExtra) # for multiple plots
library(tictoc) # for benchmarking
library(progress) # to track progress of grid search
library(dplyr) # for treating data

# Start counting computation time
tic("Duration of the whole workflow")

#################################
############# DATA ##############
#################################

# Choose data to be imported
type <- 'simulated_LSODA_log' # set 'simulated_LSODA' for simulated data using LSODA, and 'simulated_HETTMO' for simulated data using HETTMO
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

#####################################
########## INITIALIZATION ###########
#####################################

if(exists("best_params")){
  
  initial_params <-
    initialization(obs_with_noise, beta0 = real_beta[1], beta0prime = best_params$beta0prime,
                   lambda = params$lambda, gamma = params$gamma, eta = params$eta,
                   l = best_params$l, scale = 1, noise_obs = params$obs_noise,
                   noise_X = sqrt(params$obs_noise), noise_U = 0.01,
                   noise_wiener_X = best_params$noise_wiener_X, noise_wiener_U = best_params$noise_wiener_U,
                   pop = params$pop,
                   num_points_between = 0)
  
}else if(type == 'simulated_LSODA_sin'){
  
  initial_params <-
    initialization(obs_with_noise, beta0 = real_beta[1], beta0prime = 0.3,
                   lambda = params$lambda, gamma = params$gamma, eta = params$eta,
                   l = 9.7, scale = 1, noise_obs = params$obs_noise,
                   noise_X = sqrt(params$obs_noise), noise_U = 0.01,
                   noise_wiener_X = 50, noise_wiener_U = 0.01,
                   pop = params$pop,
                   num_points_between = 0)
  
  # best_params
  #   l noise_wiener_X noise_wiener_U beta0prime         SPE     NLPD      CRPS
  # 9.7             50           0.01        0.3 0.004045279 1.382473 0.3722231
  
}else if(type == 'simulated_LSODA_log'){
  
  initial_params <-
    initialization(obs_with_noise, beta0 = real_beta[1], beta0prime = 0.5,
                   lambda = params$lambda, gamma = params$gamma, eta = params$eta,
                   l = 9.55, scale = 1, noise_obs = params$obs_noise,
                   noise_X = sqrt(params$obs_noise), noise_U = 0.01,
                   noise_wiener_X = 25, noise_wiener_U = 0.01,
                   pop = params$pop, 
                   num_points_between = 0)
  
  # best_params
  #    l noise_wiener_X noise_wiener_U beta0prime         SPE     NLPD      CRPS
  # 9.55             25           0.01        0.5 0.004370289 1.381363 0.3718669
  
}

# If using data simulated with HETTMO:
# initial_params <-
#   initialization(obs, beta0 = 0.9721224, beta0prime = -1.5,
#                  lambda = 0.3703704, gamma = 0.3703704, eta = 0,
#                  l = 1.2, scale = 3, noise_obs = 100,
#                  noise_wiener_X = 1e+03, noise_wiener_U = 0.1,
#                  pop = 1e+05)

#####################################
############# INFERENCE #############
#####################################

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
if(type == 'simulated_LSODA_sin'){
  directory_res <- "~/Documents/GitHub/proboder/Results/sin" # directory of data
}else if(type == 'simulated_LSODA_log'){
  directory_res <- "~/Documents/GitHub/proboder/Results/log" # directory of data
}else if(type == 'simulated_HETTMO'){
  directory_res <- "~/Documents/GitHub/proboder/Results/hettmo" # directory of data
}

# Save results to the specified directory
save_results_as_Rdata(X_values, U_values, P_X_values, P_U_values, directory_res)

###############################
########### SCORING ###########
###############################

# ---------------------
# Extract relevant data
# ---------------------

# Load and process data from the specified directory
processed_data <- load_and_process_data(directory_res,time_grid)
U_plot <- processed_data$U_plot
X_plot <- processed_data$X_plot

# Save processed data to the specified directory
save_processed_data(U_plot, X_plot, directory_res)

# Create data frame for real beta values
real_beta_df <- data.frame(time = data_grid, real_beta = real_beta)
colnames(real_beta_df) <- c('t','beta')

# --------
# Scoring
# --------

# Compute the different scores
SPE <- squared_prediction_error(U_plot,real_beta_df)
NLPD <-negative_log_predictive_density(U_plot,real_beta_df)
CRPS <- continuous_ranked_probability_score(U_plot,real_beta_df)

# Create a data frame with the results
results <- data.frame(
  Method = c("Squared Prediction Error", "Negative Log Predictive Density", "Continuous Ranked Probability Score"),
  Value = c(SPE, NLPD, CRPS)
)

# # Create nice table using knitr and kableExtra
table <- kable(results, align = "c", caption = "Scoring Methods Results")
(styled_table <- kableExtra::kable_styling(table, bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE))

# ------
# Saving
# ------

# Save the table as a HTML file
file_path_html <- file.path(directory_res, "scoring_results.html")
save_kable(styled_table, file = file_path_html, type = "html")

# To store the table as a png file, use the manual export option from the viewer.
# Save in size 400x200.

#####################################
########### VISUALIZATION ###########
#####################################

# --------
# Plotting
# --------

# Plot simulated compartment counts and noisy simulated observations.
file_path <- file.path(directory_res, "sim-counts.pdf")
pdf(file_path, width = 8, height = 6)
plot_sim(obs, obs_with_noise)
dev.off()
plot_sim(obs, obs_with_noise)

# Plot compartment counts inferred from simulated data
file_path <- file.path(directory_res, "SEIRD-counts.pdf")
pdf(file_path, width = 10, height = 6)
plot_data_sim(obs,obs_with_noise,X_plot)
dev.off()
plot_data_sim(obs,obs_with_noise,X_plot)

# Plot compartment counts separately
plots <- plot_compartment(obs,obs_with_noise,X_plot)
for (i in 1:5) {
  file_path <- file.path(directory_res, paste0("SEIRD-counts-sep-with-CI-", i, ".pdf"))
  pdf(file_path, width = 8, height = 6)
  plot <- plots[[i]]
  print(plot)
  dev.off()
}
for (i in 1:5) {
  plot <- plots[[i]]
  print(plot)
}

# Get some fixed values to plot together with contact rate
lambda <- round(initial_params$lambda, 4)
gamma <- round(initial_params$gamma, 4)
eta <- initial_params$eta
l <- initial_params$l

# Plot simulated contact rate
file_path <- file.path(directory_res, "sim-contact-rate.pdf")
pdf(file_path, width = 8, height = 6)
plot_sim_contact_rate(real_beta_df, lambda, gamma, eta)
dev.off()
plot_sim_contact_rate(real_beta_df, lambda, gamma, eta)

# Plot contact rate inferred from simulated data
file_path <- file.path(directory_res, "inf-contact-rate-with-CI.pdf")
pdf(file_path, width = 10, height = 6)
plot_contact_rate_with_CI(U_plot, real_beta_df, lambda, gamma, eta, l)
dev.off()
plot_contact_rate_with_CI(U_plot, real_beta_df, lambda, gamma, eta, l)

# ------------
# Benchmarking
# ------------

# Total computation time
toc()