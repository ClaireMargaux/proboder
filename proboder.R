#################################### PROBODER ##################################
################################ Claire Descombes ##############################
################################################################################

# Import functions.
source('~/Documents/GitHub/proboder/functions_for_inference.R')
source('~/Documents/GitHub/proboder/saving_loading_plotting.R')

# Necessary packages.
library(Matrix)
library(numDeriv)
library(matrixcalc)

# Import data (in any case: date-S-I-D data).
directory_data <- "~/Documents/GitHub/proboder/Data" # directory of data
type <- 'simulated' # set 'real' for real data, 'simulated' for simulated data
region <- 'BE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly' (if 'real' data selected)

data <- load_data(type,region,daily_or_weekly,directory_data)
obs <- data$observations
pop <- data$population
real_beta <- data$real_beta

# Sanity check.
head(obs)
summary(obs)

#####################################
########## INITIALIZATION ###########
#####################################

# Initialize solution of SIRD-ODE and its two first derivatives
X <- as.vector(c(data= c(obs[1,2],rep(0,11))))
# Initialize latent parameter (contact rate) and its first derivative
U <- as.vector(c(0,0))

# Fixed parameters
gamma <- 0.06  # Recovery rate
eta <- 0.002   # Fatality rate
l <- 14        # Length scale

# Drift matrices
F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
F_X <- as.matrix(sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12)))

# Dispersion matrices
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- as.matrix(sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4)))

# Observation matrix (for observation of S,I and D)
H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,14)))

# Observation noise
R <- matrix(0.001, nrow = 3, ncol = 3)

# Noise of priors
P_X <- matrix(0.001, nrow = 12, ncol = 12)
P_U <- matrix(0.001, nrow = 2, ncol = 2)

#####################################
############# ALGORITHM #############
#####################################

# Function parameters
# @param recovery_rate
# @param fatality_rate
# @param length_scale
# @return contact_rate

# Data grid
data_grid <- obs[,'date']

# ODE grid
ode_grid <- data_grid # more points could be added

# Overall time grid
time_grid <- sort(unique(c(data_grid, ode_grid)))

# Arrays to store values
X_values <- matrix(data = NA, nrow = 12, ncol = length(time_grid))
U_values <- matrix(data = NA, nrow = 2, ncol = length(time_grid))
P_X_values <- array(data = NA, dim = c(12, 12, length(time_grid)))
P_U_values <- array(data = NA, dim = c(2, 2, length(time_grid)))

# Run inference.
for (loc in time_grid){
  X_values[,loc] <- as.vector(X)
  U_values[,loc] <- as.vector(U)
  P_X_values[,,loc] <- as.matrix(P_X)
  P_U_values[,,loc] <- as.matrix(P_U)
  
  # Prediction step.
  U <- as.vector(prediction_U(m_U=U,P_U=P_U,F_U=F_U,L_U=L_U)[[1]])
  P_U <- as.matrix(prediction_U(m_U=U,P_U=P_U,F_U=F_U,L_U=L_U)[[2]])
  X <- as.vector(prediction_X(m_X=X,P_X=P_X,F_X=F_X,L_X=L_X)[[1]])
  P_X <- as.matrix(prediction_X(m_X=X,P_X=P_X,F_X=F_X,L_X=L_X)[[2]])
  
  # Update of observations.
  if (any(data_grid == loc)){
    m <- as.vector(c(X,U))
    P <- matrix_P(P_X,P_U)
    y <- unlist(obs[which(obs[, 1] == loc),2:4])
    X <- as.vector(update_of_observations(m,P,y,H,R)[[1]][1:12])
    U <- as.vector(update_of_observations(m,P,y,H,R)[[1]][13:14])
    P_X <- as.matrix(update_of_observations(m,P,y,H,R)[[2]][1:12,1:12])
    P_U <- as.matrix(update_of_observations(m,P,y,H,R)[[2]][13:14,13:14])
  }
  
  # Update of states.
  if (any(ode_grid == loc)){
    P <- matrix_P(P_X,P_U)
    J <- jacobian_h(X,U,pop,gamma,eta)
    m <- c(X,U)
    h_val <- h(X,U,pop,gamma,eta)
    X <- as.vector(update_of_states(m,P,h_val,J)[[1]][1:12])
    U <- as.vector(update_of_states(m,P,h_val,J)[[1]][13:14])
    P_X <- as.matrix(update_of_states(m,P,h_val,J)[[2]][1:12,1:12])
    P_U <- as.matrix(update_of_states(m,P,h_val,J)[[2]][13:14,13:14])
  }
}

# ------------
# Save results
# ------------

# Specify directory for results
directory_res = "~/Documents/GitHub/proboder/Results"
# Save results to the specified directory
save_matrices_as_Rdata(X_values, U_values, P_X_values, P_U_values, directory_res)
  
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
U_value <- processed_data$U_scaled

# Create data frame for real beta values (if available)
real_beta_df <- data.frame(time = time_grid, real_beta = real_beta)

# --------
# Plotting
# --------

# Plot contact rate
plot_contact_rate(type, U_plot, ymin, ymax, U_value, real_beta_df)