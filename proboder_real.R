#####################################
############# WORKFLOW ##############
#####################################

# Import functions
source('~/Documents/GitHub/proboder/initialization.R')
source('~/Documents/GitHub/proboder/inference.R')
source('~/Documents/GitHub/proboder/saving_loading.R')
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

# Start counting computation time
tic("Duration of the whole workflow")

#################################
############# DATA ##############
#################################

# Choose data to be imported
type <- 'real' # set 'real' for real data, 'simulated_LSODA' for simulated data using LSODA, and 'simulated_HETTMO' for simulated data using HETTMO
region <- 'GE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'daily' # choose either 'daily' or 'weekly' (if 'real' data selected)

directory_data <- "~/Documents/GitHub/proboder/Data/real" # directory of data

# Import data
data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- as.data.frame(data$obs)
pop <- data$population

# Sanity check.
head(obs)

# Take only 30 first days.
obs <- obs[1:30,]

#####################################
########## INITIALIZATION ###########
#####################################

# If using daily data:
initial_params <-
  initialization(obs, beta0 = 0.8, beta0prime = 0.3,
                 lambda = 1/2.6, gamma = 1/2.6, eta = 0.024/15,
                 l = 7.7, scale = 1, noise_obs = 10,
                 noise_X = sqrt(10), noise_U = 0.1,
                 noise_wiener_X = 95, noise_wiener_U = 0.015,
                 pop = pop,
                 num_points_between = 0)

#####################################
############# INFERENCE #############
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
steps <- 1

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

##################################
########### PROCESSING ###########
##################################

# ---------------------
# Extract relevant data
# ---------------------

# Load and process data from the specified directory
processed_data <- load_and_process_data(directory_res,time_grid)
U_plot <- processed_data$U_plot
X_plot <- processed_data$X_plot

# Save processed data to the specified directory
save_processed_data(U_plot, X_plot, directory_res)

#####################################
########### VISUALIZATION ###########
#####################################

# --------
# Plotting
# --------

# Plot compartment counts inferred from real data
plot_data_real(obs, X_plot)

# Plot compartment counts separately
plots <- plot_compartment_real(obs,X_plot)
for (i in 1:5) {
  plot <- plots[[i]]
  print(plot)
}

# Get some fixed values to plot together with contact rate
lambda <- initial_params$lambda
gamma <- initial_params$gamma
eta <- initial_params$eta
l <- initial_params$l

# Plot contact rate inferred from real data
plot_contact_rate_with_CI_real(U_plot, lambda, gamma, eta, l) 

# ------------
# Benchmarking
# ------------

# Total computation time
toc()