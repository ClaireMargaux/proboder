#################################### PROBODER ##################################
################################ Claire Descombes ##############################

# Import functions.
source('~/Documents/GitHub/proboder/functions.R')

# Necessary packages.
library(Matrix)
library(numDeriv)
library(matrixcalc)

# Import data (in any case: date-S-I-D data).
directory <- "~/Documents/GitHub/proboder/Data" # directory of data
type <- 'simulated' # set 'real' for real data, 'simulated' for simulated data
region <- 'BE' # 'BE' or 'GE' available (if 'real' data selected)
daily_or_weekly <- 'weekly' # choose either 'daily' or 'weekly' (if 'real' data selected)

data <- load_data(type,region,daily_or_weekly,directory)
obs <- data$observations
pop <- data$population

# Sanity check.
head(obs)
summary(obs)

#####################################
########## INITIALIZATION ###########
#####################################

# X: initialization of solution of SIRD-ODE and its two first derivatives.
X <- as.vector(c(data= c(obs[1,2],rep(0,11))))
# U: initialization of latent parameter (contact rate) and its first derivative.
U <- as.vector(c(0,0))

# Fixed parameters.
gamma <- 0.06 # recovery_rate
eta <- 0.002 # fatality_rate
l <- 14 # length_scale

# Drift matrices.
F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
F_X <- as.matrix(sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12)))

# Dispersion matrices.
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4))
L_X <- as.matrix(L_X)

# Observation matrix (for observation of S,I and D).
H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,14)))

# Observation noise.
R <- matrix(0.001, nrow = 3, ncol = 3)

# Noise of priors.
P_X <- matrix(0.001, nrow = 12, ncol = 12)
P_U <- matrix(0.001, nrow = 2, ncol = 2)

#####################################
############# ALGORITHM #############
#####################################

#' @param recovery_rate
#' @param fatality_rate
#' @param length_scale
#' @return contact_rate

# Data grid.
data_grid <- obs[,'date']

# ODE grid.
ode_grid <- data_grid # more points could be added

# Overall time grid.
time_grid <- sort(unique(c(data_grid, ode_grid)))

# Arrays to store values of X, P_X and U, P_U.
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
    y <- obs[which(obs[, 1] == loc),2:4]
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

#####################################
########### VISUALIZATION ###########
#####################################

library(ggplot2)

# Extract relevant data
U_scaled <- sigmoid(U_values[1,])
U_plot <- data.frame(time = time_grid, U_value = U_scaled)
P_plot <- data.frame(
  time = time_grid,
  ymin = U_scaled - sqrt(P_U_values[1, 1, ]),
  ymax = U_scaled + sqrt(P_U_values[1, 1, ])
)
real_beta_df <- data.frame(time = time_grid, real_beta = real_beta)

# Plotting
ggplot() +
  geom_line(data = U_plot, aes(x = time, y = U_value, color = "Estimated Contact Rate"), size = 1) +
  geom_ribbon(data = P_plot, aes(x = time, ymin = ymin, ymax = ymax), fill = "lightgreen", alpha = 0.5) +
  geom_line(data = real_beta_df, aes(x = time, y = real_beta, color = "Real Contact Rate"), linetype = "dashed") +
  labs(x = "Time", y = "Contact rate", title = "Contact rate with Error Area",
       color = "Legend") +  
  scale_color_manual(values = c("Estimated Contact Rate" = "darkgreen", "Real Contact Rate" = "lightblue4"),
                     labels = c("Estimated Contact Rate", "Real Contact Rate")) +  # Specify legend labels
  coord_cartesian(ylim = c(-1, 2)) +
  theme_minimal() +
  theme(legend.position = "top")
