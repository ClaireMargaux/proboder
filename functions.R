#################################### PROBODER ##################################
################################ Claire Descombes ##############################
################################### Functions ##################################

# Link function.
sigmoid <- function(z){
  return(1/(1+exp(-z)))
}

# ODE.
f <- function(X0,beta,pop,gamma,eta){
  # Arguments:
  #   arg1: X0, matrix(nrow=4, ncol=1), solution of the ODE
  #   arg2: beta, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: pop, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: sol, numeric, evaluation of the ODE f, vector of size 4
  
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  P <- pop
  
  S_out <- -beta*S*I/P
  I_out <- beta*S*I/P - gamma*I - eta*I
  R_out <- gamma*I
  D_out <- eta*I
  
  sol <- c(S_out,I_out,R_out,D_out)
  
  return(sol)
}

# Measurement model.
h <- function(X,U,pop,gamma,eta){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, matrix(nrow=2, ncol=1), latent parameter of the ODE and it's first derivative
  #   arg3: pop, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: evaluation of the measurement model h, matrix(nrow=4, ncol=1)
  
  beta <- U[1]
  beta <- sigmoid(beta) # rescaling of beta to [0,1]
  
  X0 <- X[1:4]
  X1 <- X[5:8]
  ODE <- f(X0,beta,pop,gamma,eta)
  
  sol <- matrix(X1 - ODE, nrow = 4, ncol = 1)
  
  return(sol)
}

#' Jacobian of measurement model using jacobian function.
jacobian_h <- function(X,U,pop,gamma,eta){
  h1 <- function(X) h(X,U=U,pop=pop,gamma=gamma,eta=eta)
  h2 <- function(U) h(X=X,U,pop=pop,gamma=gamma,eta=eta)
  out <- cbind(jacobian(h1, x = X),jacobian(h2, x = U))
  return(out)
}

# Jacobian of f 'by hand'.
jacobian_U <- function(X,U,pop,gamma,eta){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: pop, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: matrix(nrow=2, ncol=2), Jacobian of ODE f wrt U
  
  X0 <- X[1:4]
  
  P <- pop
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  beta <- U
  
  d_dbeta <- matrix(
    data=c(-S * I / P, S * I / P, 0, 0),
    nrow=1, ncol=4
  )
  
  sol <- d_dbeta
  
  return(sol)
}

jacobian_X <- function(X,U,pop,gamma,eta){
  # Arguments:
  #   arg1: X, matrix(nrow=12, ncol=1), solution of the ODE and its 2 first derivatives
  #   arg2: U, numeric with values in [0,1], latent parameter of the ODE
  #   arg3: pop, integer, total population
  #   arg4: gamma, numeric, recovery rate
  #   arg5: eta, numeric, fatality rate
  #
  # Returns:
  #   output: matrix(nrow=4, ncol=4), Jacobian of ODE f wrt X0
  
  X0 <- X[1:4]
  
  P <- pop
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  beta <- U
  
  d_dS <- matrix(
    data=c(-beta * I / P, -beta * S / P, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dI <- matrix(
    data=c(beta * I / P, beta * S / P - gamma - eta, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dR <- matrix(
    data=c(0, gamma, 0, 0),
    nrow=1, ncol=4
  )
  
  d_dD <- matrix(
    data=c(0, eta, 0, 0),
    nrow=1, ncol=4
  )
  
  sol <- rbind(d_dS,d_dI,d_dR,d_dD)
  
  return(sol)
}

# Prediction step.
prediction_U <- function(m_U,P_U,F_U,L_U){
  # Arguments:
  #   arg1: m_U, numeric, (previous) mean of U
  #   arg2: P_U, numeric, (previous) (co)variance of U
  #   arg3: F_U, matrix(nrow=2, ncol=2), drif matrix of U
  #   arg4: L_U, matrix(nrow=2, ncol=1), dispersion matrix of U
  #
  # Returns:
  #   output: predicted mean and covariance U
  
  exp_F_U <- expm(F_U)
  m_U_out <- exp_F_U %*% m_U
  
  top_row_U <- cbind(F_U,L_U %*% t(L_U))
  bottom_row_U <- cbind(matrix(0,nrow=2,ncol=2),-t(F_U))
  Gamma_U <- rbind(top_row_U,bottom_row_U)
  exp_Gamma_U <- expm(Gamma_U)
  
  M_1_U <- exp_Gamma_U[1:2,1:2]
  M_2_U <- exp_Gamma_U[1:2,3:4]
  
  B_U <- M_2_U %*% t(M_1_U)
  
  P_U_out <- exp_F_U %*% P_U %*% t(exp_F_U) + B_U
  
  out <- list(m_U_out,P_U_out)
  
  return(out)
}

prediction_X <- function(m_X,P_X,F_X,L_X){
  # Arguments:
  #   arg1: m_X, vector(nrow=13), (previous) mean of X
  #   arg2: P_X, matrix(nrow=12, ncol=12), (previous) covariance of X
  #   arg3: F_X, matrix(nrow=12, ncol=12), drif matrix of X
  #   arg4: L_X, matrix(nrow=12, ncol=4), drif matrix of X
  #
  # Returns:
  #   output: predicted mean and covariance of X
  
  exp_F_X <- expm(F_X)
  m_X_out <- exp_F_X %*% m_X
  
  top_row_X <- cbind(F_X,L_X %*% t(L_X))
  bottom_row_X <- cbind(matrix(0,nrow=12,ncol=12),-t(F_X))
  Gamma_X <- rbind(top_row_X,bottom_row_X)
  exp_Gamma_X <- expm(Gamma_X)
  
  M_1_X <- exp_Gamma_X[1:12,1:12]
  M_2_X <- exp_Gamma_X[1:12,13:24]
  
  B_X <- M_2_X %*% t(M_1_X)
  
  P_X_out <- exp_F_X %*% P_X %*% t(exp_F_X) + B_X
  
  out <- list(m_X_out,P_X_out)
  
  return(out)
}

# Update step on tau_OBS.
update_of_observations <- function(m,P,y,H,R){
  # Arguments:
  #   arg1: m, vector(nrow=14), predicted mean of U and X
  #   arg2: P, matrix(nrow=14, ncol=14), predicted covariance of U and X
  #   arg3: y, vector(nrow=3), observations
  #   arg4: H, matrix(nrow=3, ncol=14), observation matrix
  #   arg5: R, matrix(nrow=3, ncol=3), observation noise
  #
  # Returns:
  #   output: vector(ncol=14) and matrix(nrow=14, ncol=14), updated
  #           mean and covariances of U and X
  
  v <- y - H %*% m # residual
  S <- H %*% P %*% t(H) + R # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(H) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

# Update step on tau_ODE.
update_of_states <- function(m,P,h,J){
  # Arguments:
  #   arg1: m, vector(nrow=14), predicted mean of U and X
  #   arg2: P, matrix(nrow=14, ncol=14), predicted covariance of U and X
  #   arg3, h, vector(nrow=4), measurement model
  #   arg4: J, matrix(nrow=4, ncol=14), Jacobian of the measurement model
  #
  # Returns:
  #   output: vector(ncol=14) and matrix(nrow=14, ncol=14), updated
  #           mean and covariances of U and X
  
  v <- -h # residual
  S <- J %*% P %*% t(J) # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(J) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

matrix_P <- function(P_X,P_U){
  # Arguments:
  #   arg1: P_X, matrix(nrow=12, ncol=12), predicted covariance of X
  #   arg2: P_U, matrix(nrow=2, ncol=2), predicted covariance of U
  #
  # Returns:
  #   output: matrix(nrow=14, ncol=14), predicted covariance of X and U
  
  P_X <- as.matrix(P_X)
  P_U <- as.matrix(P_U)
  
  nrow_P <- nrow(P_X) + nrow(P_U)
  ncol_P <- ncol(P_X) + ncol(P_U)
  P <- matrix(0, nrow = nrow_P, ncol = ncol_P)
  P[1:nrow(P_X), 1:ncol(P_X)] <- P_X  # Top left corner
  P[(nrow_P - nrow(P_U) + 1):nrow_P, (ncol_P - ncol(P_U) + 1):ncol_P] <- P_U  # Bottom right corner
  
  return(P)
}

load_data <- function(type, region, daily_or_weekly, directory) {
  # Check if type is valid
  if (!(type %in% c('simulated', 'real'))) {
    stop("Error: Invalid type. Type must be 'simulated' or 'real'.")
  }
  
  # Check if region is valid
  if (!(region %in% c('BE', 'GE'))) {
    stop("Error: Invalid region. Region must be 'BE' or 'GE'.")
  }
  
  # Check if daily_or_weekly is valid
  if (!(daily_or_weekly %in% c('daily', 'weekly'))) {
    stop("Error: Invalid value for daily_or_weekly. Must be 'daily' or 'weekly'.")
  }
  
  # Load data based on type
  if (type == 'simulated') {
    load(file.path(directory, "simulated_data.Rdata"))
    obs <- observations_simulate
    population <- readRDS(file.path(directory, "simulated_pop.Rds"))
  } else {
    region_filename <- paste0("real_data_", region, "_", daily_or_weekly, ".Rdata")
    population_filename <- paste0("real_pop_", region, "_", daily_or_weekly, ".Rds")
    load(file.path(directory, region_filename))
    obs <- observations
    population <- readRDS(file.path(directory, population_filename))
  }
  
  # Return loaded data
  return(list(observations = obs, population = population))
}