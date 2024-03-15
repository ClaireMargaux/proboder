#################################### PROBODER ##################################
################################ Claire Descombes ##############################
############################ Functions for inference ###########################

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

#' Prediction function for latent parameter U or state variable X
#'
#' This function calculates the predicted mean and covariance of either U or X using the given parameters.
#'
#' @param m Numeric vector, the (previous) mean of U or X.
#' @param P Numeric matrix, the (previous) covariance of U or X.
#' @param F Numeric matrix, drift matrix of U or X.
#' @param L Numeric matrix, dispersion matrix of U or X.
#' @return A list containing the predicted mean and covariance.
#' @export
prediction <- function(m, P, drift, dispersion) {
  # Arguments:
  #   m: Numeric vector, the (previous) mean of U or X.
  #   P: Numeric matrix, the (previous) covariance of U or X.
  #   drift: Numeric matrix, drift matrix of U or X.
  #   dispersion: Numeric matrix, dispersion matrix of U or X.
  #
  # Returns:
  #   A list containing the predicted mean and covariance.
  
  exp_F <- expm(drift)
  m_out <- exp_F %*% m
  m_out <- as.vector(m_out)
  
  n <- nrow(drift)
  j <- ncol(dispersion)
  Q <- diag(0.1, nrow = j, ncol = j)
  top_row <- cbind(drift, dispersion %*% Q %*% t(dispersion))
  #top_row <- cbind(drift, dispersion %*% t(dispersion))
  bottom_row <- cbind(matrix(0, nrow = n, ncol = n), -t(drift))
  Gamma <- rbind(top_row, bottom_row)
  exp_Gamma <- expm(Gamma)
  
  M_1 <- exp_Gamma[1:n, 1:n]
  M_2 <- exp_Gamma[1:n, (n + 1):(2 * n)]
  
  G <- M_2 %*% t(M_1)
  
  P_out <- exp_F %*% P %*% t(exp_F) + G
  P_out <- as.matrix(P_out)
  
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
