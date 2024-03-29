#' Sigmoid function
#'
#' Calculate the sigmoid (logistic) function for the given input.
#'
#' @param z Numeric input.
#' @return Sigmoid value of the input.
#' @examples
#' sigmoid(0) # Returns 0.5
#' sigmoid(1) # Returns 0.7310586
#' @export
sigmoid <- function(z) {
  return(1 / (1 + exp(-z)))
}

#' Logit function
#'
#' Calculate the logit (inverse of sigmoid) function for the given input.
#'
#' @param p Numeric input.
#' @return Logit value of the input.
#' @examples
#' logit(0.5) # Returns 0
#' logit(0.7310586) # Returns 1
#' @export
logit <- function(p) {
  return(log(p / (1 - p)))
}

#' ODE function
#'
#' This function evaluates the ordinary differential equation (ODE) for a given set of parameters.
#'
#' @param X0 Numeric vector of size 4 representing the solution of the ODE.
#' @param beta Numeric value representing the latent parameter of the ODE.
#' @param pop Integer representing the total population.
#' @param gamma Numeric value representing the recovery rate.
#' @param eta Numeric value representing the fatality rate.
#' @return Numeric vector of size 4 representing the evaluation of the ODE.
#' @export
f <- function(X0, beta, pop, gamma, eta) {
  # Arguments:
  #   X0: Numeric vector of size 4, solution of the ODE
  #   beta: Numeric with values in [0,1], latent parameter of the ODE
  #   pop: Integer, total population
  #   gamma: Numeric, recovery rate
  #   eta: Numeric, fatality rate
  #
  # Returns:
  #   sol: Numeric vector of size 4 representing the evaluation of the ODE
  
  S <- X0[1]
  I <- X0[2]
  R <- X0[3]
  D <- X0[4]
  P <- pop
  
  S_out <- - beta * S * I / P
  I_out <- beta * S * I / P - gamma * I - eta * I
  R_out <- gamma * I
  D_out <- eta * I
  
  sol <- c(S_out, I_out, R_out, D_out)
  
  return(sol)
}

#' Jacobian of the ODE function using jacobian function
#'
#' This function calculates the Jacobian matrix of the ODE function using the jacobian function.
#'
#' @param X0 Numeric vector of size 4 representing the solution of the ODE.
#' @param beta Numeric value representing the latent parameter of the ODE.
#' @param pop Integer representing the total population.
#' @param gamma Numeric value representing the recovery rate.
#' @param eta Numeric value representing the fatality rate.
#' @return The Jacobian matrix of the ODE function.
#' @export
jacobian_f <- function(X0, beta, pop, gamma, eta) {
  # Arguments:
  #   X0: Numeric vector of size 4, solution of the ODE
  #   beta: Numeric with values in [0,1], latent parameter of the ODE
  #   pop: Integer, total population
  #   gamma: Numeric, recovery rate
  #   eta: Numeric, fatality rate
  #
  # Returns:
  #   The Jacobian matrix of the ODE function
  
  f1 <- function(X0) f(X0, beta = beta, pop = pop, gamma = gamma, eta = eta)
  f2 <- function(beta) f(X0 = X0, beta, pop = pop, gamma = gamma, eta = eta)
  out <- cbind(jacobian(f1, x = X0), jacobian(f2, x = beta))
  return(out)
}

#' Measurement model
#'
#' This function evaluates the measurement model for a given set of parameters.
#'
#' @param X Numeric vector of size 12 representing the solution of the ODE and its 2 first derivatives.
#' @param U Numeric vector of size 2 representing the latent parameter of the ODE and its first derivative.
#' @param pop Integer representing the total population.
#' @param gamma Numeric representing the recovery rate.
#' @param eta Numeric representing the fatality rate.
#' @return Numeric vector of size 4 representing the evaluation of the measurement model.
#' @export
h <- function(X, U, pop, gamma, eta) {
  # Arguments:
  #   X: Numeric vector of size 12, solution of the ODE and its 2 first derivatives
  #   U: Numeric vector of size 2, latent parameter of the ODE and its first derivative
  #   pop: Integer, total population
  #   gamma: Numeric, recovery rate
  #   eta: Numeric, fatality rate
  #
  # Returns:
  #   sol: Numeric vector of size 4 representing the evaluation of the measurement model
  
  beta <- U[1]
  beta <- sigmoid(beta) # rescaling of beta to [0,1]
  
  X0 <- X[1:4]
  X1 <- X[5:8]
  ODE <- f(X0, beta, pop, gamma, eta)
  
  sol <- c(X1 - ODE)
  
  return(sol)
}

#' Jacobian of measurement model using jacobian function
#'
#' This function calculates the Jacobian matrix of the measurement model using the jacobian function.
#'
#' @param X Vector of length 12 representing the solution of the ODE and its 2 first derivatives.
#' @param U Vector of length 2 representing the latent parameter of the ODE and its first derivative.
#' @param pop Integer representing the total population.
#' @param gamma Numeric representing the recovery rate.
#' @param eta Numeric representing the fatality rate.
#' @return The Jacobian matrix of the measurement model.
#' @export
jacobian_h <- function(X, U, pop, gamma, eta) {
  # Arguments:
  #   X: Vector of length 12 representing the solution of the ODE and its 2 first derivatives
  #   U: Vector of length 2 representing the latent parameter of the ODE and its first derivative
  #   pop: Integer, total population
  #   gamma: Numeric, recovery rate
  #   eta: Numeric, fatality rate
  #
  # Returns:
  #   The Jacobian matrix of the measurement model
  
  h1 <- function(X) h(X, U = U, pop = pop, gamma = gamma, eta = eta)
  h2 <- function(U) h(X = X, U, pop = pop, gamma = gamma, eta = eta)
  out <- cbind(jacobian(h1, x = X), jacobian(h2, x = U))
  return(out)
}

#' Prediction function for latent parameter U or state variable X
#'
#' This function calculates the predicted mean and covariance of either U or X using the given parameters.
#'
#' @param m Numeric vector, the (previous) mean of U or X.
#' @param P Numeric matrix, the (previous) covariance of U or X.
#' @param drift Numeric matrix, drift matrix of U or X.
#' @param dispersion Numeric matrix, dispersion matrix of U or X.
#' @param noise_wiener Numeric matrix, noise of the Wiener process.
#' @return A list containing the predicted mean and covariance.
#' @export
prediction <- function(m, P, drift, dispersion, noise_wiener) {
  # Arguments:
  #   m: Numeric vector, the (previous) mean of U or X.
  #   P: Numeric matrix, the (previous) covariance of U or X.
  #   drift: Numeric matrix, drift matrix of U or X.
  #   dispersion: Numeric matrix, dispersion matrix of U or X.
  #   noise_wiener: Numeric matrix, noise of the Wiener process.
  #
  # Returns:
  #   A list containing the predicted mean and covariance.
  
  n <- nrow(drift)
  
  exp_F <- expm(drift)
  m_out <- as.vector(exp_F %*% m)

  top_row <- cbind(drift, dispersion %*% noise_wiener %*% t(dispersion))
  bottom_row <- cbind(matrix(0, nrow = n, ncol = n), -t(drift))
  Gamma <- rbind(top_row, bottom_row)
  exp_Gamma <- expm(Gamma)
  
  M_1 <- exp_Gamma[1:n, 1:n]
  M_2 <- exp_Gamma[1:n, (n + 1):(2 * n)]
  M_3 <- exp_Gamma[(n + 1):(2 * n), (n + 1):(2 * n)]
  
  # Alternative 1: compute Qd
  Qd <- M_2 %*% t(M_1)
  P_out <- as.matrix(exp_F %*% P %*% t(exp_F) + Qd)

  # Alternative 2: matrix fraction decomposition
  #M_3_inv <- svd.inverse(as.matrix(M_3))
  #P_out <- as.matrix((M_1 %*% P + M_2) %*% M_3_inv)
  
  out <- list(m_out, P_out)
  
  return(out)
}

#' Update step for observations (Kalman Filter)
#'
#' This function updates the mean and covariance of U and X using observations.
#'
#' @param m Numeric vector, predicted mean of U and X.
#' @param P Numeric matrix, predicted covariance of U and X.
#' @param y Numeric vector, observations.
#' @param H Numeric matrix, observation matrix.
#' @param R Numeric matrix, observation noise.
#' @return A list containing the updated mean and covariance of U and X.
#' @export
update_of_observations <- function(m, P, y, H, R) {
  # Arguments:
  #   m: Numeric vector, predicted mean of U and X.
  #   P: Numeric matrix, predicted covariance of U and X.
  #   y: Numeric vector, observations.
  #   H: Numeric matrix, observation matrix.
  #   R: Numeric matrix, observation noise.
  #
  # Returns:
  #   A list containing the updated mean and covariance of U and X.
  
  v <- y - H %*% m # residual
  S <- H %*% P %*% t(H) + R # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(H) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

#' Update step for states (Extended Kalman Filter)
#'
#' This function updates the mean and covariance of U and X using the measurement model and its Jacobian.
#'
#' @param m Numeric vector, predicted mean of U and X.
#' @param P Numeric matrix, predicted covariance of U and X.
#' @param h Numeric vector, measurement model.
#' @param J Numeric matrix, Jacobian of the measurement model.
#' @return A list containing the updated mean and covariance of U and X.
#' @export
update_of_states <- function(m, P, h, J) {
  # Arguments:
  #   m: Numeric vector, predicted mean of U and X.
  #   P: Numeric matrix, predicted covariance of U and X.
  #   h: Numeric vector, measurement model.
  #   J: Numeric matrix, Jacobian of the measurement model.
  #
  # Returns:
  #   A list containing the updated mean and covariance of U and X.
  
  v <- -h # residual
  S <- J %*% P %*% t(J) # innovation covariance
  S_inv <- svd.inverse(S)
  K <- P %*% t(J) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  out <- list(m_out,P_out)
  
  return(out)
}

#' Combine covariance matrices P_U and P_X
#'
#' This function combines the predicted covariance matrices of U and X into a single covariance matrix.
#'
#' @param P_U Numeric matrix, predicted covariance of U.
#' @param P_X Numeric matrix, predicted covariance of X.
#' @return A single matrix representing the combined covariance of X and U.
#' @export
matrix_P <- function(P_U, P_X) {
  # Arguments:
  #   P_U: Numeric matrix, predicted covariance of U.
  #   P_X: Numeric matrix, predicted covariance of X.
  #
  # Returns:
  #   A single matrix representing the combined covariance of U and X.
  
  P_U <- as.matrix(P_U)
  P_X <- as.matrix(P_X)
  
  ncol_P <- ncol(P_X) + ncol(P_U)
  nrow_P <- nrow(P_X) + nrow(P_U)
  P <- matrix(0, nrow = nrow_P, ncol = ncol_P)
  P[1:nrow(P_U), 1:ncol(P_U)] <- P_U  # Top left corner
  P[(nrow_P - nrow(P_X) + 1):nrow_P, (ncol_P - ncol(P_X) + 1):ncol_P] <- P_X  # Bottom right corner
  
  return(P)
}

#' Perform inference using a state-space model
#'
#' This function performs inference over a specified time grid using a state-space model.
#' It iteratively updates the state estimates and covariances based on observations and predictions.
#'
#' @param time_grid Numeric vector, time grid for the inference.
#' @param obs Data frame, contains the dates and the compartment counts.
#' @param initial_params List of initial parameters obtained from the initialization function.
#'
#' @return A list with the inferred values of X, U, P_X, and P_U.
#' @export
inference <- function(time_grid, obs, initial_params){
  # Arguments:
  #   time_grid: Numeric vector, time grid for the inference.
  #   obs: Data frame, contains the dates and the compartment counts.
  #   initial_params: List of initial parameters obtained from the initialization function, has the following entries:
  #
  #   X: Vector of length 12 representing the solution of the ODE and its 2 first derivatives.
  #   U: Vector of length 2 representing the latent parameter of the ODE and its first derivative.
  #   P_U: Numeric matrix, predicted covariance of U.
  #   P_X: Numeric matrix, predicted covariance of X.
  #   F_U: Numeric matrix, drift matrix of U.
  #   F_X: Numeric matrix, drift matrix of X.
  #   L_U: Numeric matrix, dispersion matrix of U.
  #   L_U: Numeric matrix, dispersion matrix of U.
  #   noise_wiener_U: Numeric matrix, noise of the Wiener process driving U.
  #   noise_wiener_X: Numeric matrix, noise of the Wiener process driving X.
  #   H: Numeric matrix, noise of the observations.
  #   pop: Integer, total population.
  #   gamma: Numeric, recovery rate.
  #   eta: Numeric, fatality rate.
  #
  # Returns:
  #   A list with the inferred values of X, U, P_X and P_U.
  
  X <- initial_params$X
  U <- initial_params$U
  P_X <- initial_params$P_X
  P_U <- initial_params$P_U
  F_X <- initial_params$F_X
  F_U <- initial_params$F_U
  L_X <- initial_params$L_X
  L_U <- initial_params$L_U
  noise_wiener_X <- initial_params$noise_wiener_X
  noise_wiener_U <- initial_params$noise_wiener_U
  R <- initial_params$R
  H <- initial_params$H
  pop <- initial_params$pop
  gamma <- initial_params$gamma
  eta <- initial_params$eta
  
  # Initialize lists to store inferred values
  X_values <- matrix(data = NA, nrow = length(X), ncol = length(time_grid))
  U_values <- matrix(data = NA, nrow = length(U), ncol = length(time_grid))
  P_X_values <- array(data = NA, dim = c(length(X), length(X), length(time_grid)))
  P_U_values <- array(data = NA, dim = c(length(U), length(U), length(time_grid)))
  
  # Iterate over each time point in the time grid
  for (i in seq_along(time_grid)) {
    loc <- time_grid[i]
    
    # Store current values and covariances
    X_values[, i] <- X
    U_values[, i] <- U
    P_X_values[, , i] <- P_X
    P_U_values[, , i] <- P_U
    
    # Prediction step
    U <- as.vector(prediction(m = U, P = P_U, drift = F_U, dispersion = L_U, noise_wiener = noise_wiener_U)[[1]])
    P_U <- as.matrix(prediction(m = U, P = P_U, drift = F_U, dispersion = L_U, noise_wiener = noise_wiener_U)[[2]])
    X <- as.vector(prediction(m = X, P = P_X, drift = F_X, dispersion = L_X, noise_wiener = noise_wiener_X)[[1]])
    P_X <- as.matrix(prediction(m = X, P = P_X, drift = F_X, dispersion = L_X, noise_wiener = noise_wiener_X)[[2]])
    
    # Update of observations if available at the current time point
    if (any(obs$date == loc)) {
      m <- c(X, U)
      P <- matrix_P(P_X, P_U)
      y <- unlist(obs[obs$date == loc, 2:4])
      updated_obs <- update_of_observations(m = m, P = P, y = y, H = H, R = R)
      X <- as.vector(updated_obs[[1]][1:12])
      U <- as.vector(updated_obs[[1]][13:14])
      P_X <- as.matrix(updated_obs[[2]][1:12, 1:12])
      P_U <- as.matrix(updated_obs[[2]][13:14, 13:14])
    }
    
    # Update of states if required at the current time point
    if (any(ode_grid == loc)) {
      P <- matrix_P(P_X, P_U)
      J <- jacobian_h(X, U, pop, gamma, eta)
      m <- c(X, U)
      h_val <- h(X, U, pop, gamma, eta)
      updated_states <- update_of_states(m = m, P = P, h = h_val, J = J)
      X <- as.vector(updated_states[[1]][1:12])
      U <- as.vector(updated_states[[1]][13:14])
      P_X <- as.matrix(updated_states[[2]][1:12, 1:12])
      P_U <- as.matrix(updated_states[[2]][13:14, 13:14])
    }
  }
  
  # Return inferred values
  list(X_values = X_values, U_values = U_values, P_X_values = P_X_values, P_U_values = P_U_values)
}