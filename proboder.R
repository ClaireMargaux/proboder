#################################### PROBODER ##################################
################################ Claire Descombes ##############################
################################################################################

# Import functions
source('~/Documents/GitHub/proboder/functions_for_inference.R')
source('~/Documents/GitHub/proboder/saving_loading_plotting.R')

# Necessary packages
library(Matrix)
library(numDeriv)
library(matrixcalc)
library(greybox)
library(ggplot2)

# Choose data to be imported (in case 'real': date-S-I-D, in case 'simulated': date-S-I-R)
directory_data <- "~/Documents/GitHub/proboder/Data" # directory of data
type <- 'simulated' # set 'real' for real data, 'simulated' for simulated data
region <- 'BE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly' (if 'real' data selected)

# Import data
data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- data$observations
pop <- data$population
if(type == 'simulated'){
  real_beta <- data$real_beta
}

# Sanity check.
head(obs)
summary(obs)

#####################################
########## INITIALIZATION ###########
#####################################

# Initialize latent parameter (contact rate) and its first derivative
U <- as.vector(c(
  logit(0.99), # Initial value for contact rate
  0) # Initial values for 1st derivative
)

# Initialize solution of SIRD-ODE and its two first derivatives
X <- as.vector(c(
  c(obs[1,2], rep(0,3), # Initial values for S, I, and R or D
  rep(0.1,4), # Initial values for 1st derivatives
  rep(0,4))) # Initial values for 2nd derivatives
)

# Fixed parameters
gamma <- 0.4    # Recovery rate
eta <- 0.002    # Fatality rate
l <- 2          # Length scale

# Drift matrices
F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
F_X <- as.matrix(sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12)))

# Dispersion matrices
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- as.matrix(sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4)))

# Observation matrix (for observation of S,I, and R or D)
if(type == 'real'){
  H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,14)))
}else if(type == 'simulated'){
  H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,3), x = 1, dims = c(3,14)))
}else{
  print('Wrong type!')
}

# Observation noise
R <- cov(obs[,-1])

# Noise of priors
P_X <- matrix(100, nrow = 12, ncol = 12)
P_U <- matrix(100, nrow = 2, ncol = 2)

# Noise of Wiener process
noise_wiener_U <- diag(10, nrow = ncol(L_U), ncol = ncol(L_U))
noise_wiener_X <- diag(10, nrow = ncol(L_X), ncol = ncol(L_X))

#####################################
############# ALGORITHM #############
#####################################

# Data grid
data_grid <- obs[,'date']

# ODE grid
ode_grid <- data_grid # more points could be added

# Overall time grid
time_grid <- sort(unique(c(data_grid, ode_grid)))

# Run inference.
inference_results <- inference(time_grid, obs,
                               X, U, P_X, P_U, 
                               F_X, F_U, L_X, L_U, 
                               noise_wiener_X, noise_wiener_U,
                               H, pop, gamma, eta)

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

#####################################
########### VISUALIZATION ###########
#####################################

# ---------------------
# Extract relevant data
# ---------------------

# Load and process data from the specified directory
processed_data <- load_and_process_data(directory_res,time_grid)
U_plot <- processed_data$U_plot
P_plot <- processed_data$P_plot
ymin <- P_plot$ymin
ymax <- P_plot$ymax
U_scaled <- processed_data$U_scaled

# Save process data to the specified directory
save_processed_data(U_plot, P_plot, ymin, ymax, U_scaled, directory_res)

# Create data frame for real beta values (if available)
if(type=='simulated'){
  real_beta_df <- data.frame(time = time_grid, real_beta = real_beta)
}

# --------
# Plotting
# --------

# Plot data
plot_data(obs,type)

# Plot contact rate
plot_contact_rate(type, U_plot, ymin, ymax, U_scaled, real_beta_df, gamma, eta, l)
