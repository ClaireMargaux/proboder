#' (Scaled) sigmoid function
#'
#' Calculate the sigmoid (logistic) function for the given input, up to a scaling factor of 5.
#'
#' @param z Numeric input.
#' @return Sigmoid value of the input.
#' @export
sigmoid <- function(z) {
  return(5 / (1 + exp(-z)))
}

#' (Scaled) logit function
#'
#' Calculate the logit (inverse of sigmoid) function for the given input, up to a scaling factor of 5.
#'
#' @param p Numeric input.
#' @return Logit value of the input.
#' @export
logit <- function(p) {
  return(-log(5/p - 1))
}

#' ODE function
#'
#' This function evaluates the ordinary differential equation (ODE) for a given set of parameters.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param X0 Numeric vector representing the solution of the ODE.
#' @param beta Numeric value representing the latent parameter of the ODE.
#' @param pop Integer representing the total population.
#' @param lambda Numeric value representing the latency rate.
#' @param gamma Numeric value representing the recovery rate.
#' @param eta Numeric value representing the fatality rate.
#' @return Numeric vector representing the evaluation of the ODE.
#' @export
f <- function(model, X0, beta, pop, lambda, gamma, eta = 0) {
  # Arguments:
  #   model: String, type of model to be used ('SEIRD' or 'SEIR' available)
  #   X0: Numeric vector, solution of the ODE
  #   beta: Numeric with value in R+, latent parameter of the ODE
  #   pop: Integer, total population
  #   lambda: Numeric with value in R+, latency rate
  #   gamma: Numeric with value in R+, recovery rate
  #   eta: Numeric with value in R+, fatality rate
  #
  # Returns:
  #   sol: Numeric vector representing the evaluation of the ODE

  if (model == 'SEIRD') {
    dS <- - beta * X0[1] * X0[3] / pop
    dE <- beta * X0[1] * X0[3] / pop - lambda * X0[2]
    dI <- lambda * X0[2] - gamma * X0[3] - eta * X0[3]
    dR <- gamma * X0[3]
    dD <- eta * X0[3]
    res <- c(dS, dE, dI, dR, dD) 
  }
  
  if (model == 'SEIR') {
    dS <- - beta * X0[1] * X0[3] / pop
    dE <- beta * X0[1] * X0[3] / pop - lambda * X0[2]
    dI <- lambda * X0[2] - gamma * X0[3] - eta * X0[3]
    dR <- gamma * X0[3]
    res <- c(dS, dE, dI, dR)
  }
  
  return(res)
}

#' Jacobian of the ODE function using jacobian function
#'
#' This function calculates the Jacobian matrix of the ODE function using the jacobian function.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param X0 Numeric vector representing the solution of the ODE.
#' @param beta Numeric value representing the latent parameter of the ODE.
#' @param pop Integer representing the total population.
#' @param lambda Numeric value representing the latency rate.
#' @param gamma Numeric value representing the recovery rate.
#' @param eta Numeric value representing the fatality rate.
#' @return The Jacobian matrix of the ODE function.
#' @export
jacobian_f <- function(model, X0, beta, pop, lambda, gamma, eta = 0) {
  # Arguments:
  #   model: String, type of model to be used ('SEIRD' or 'SEIR' available)
  #   X0: Numeric vector, solution of the ODE
  #   beta: Numeric with values in R+, latent parameter of the ODE
  #   pop: Integer, total population
  #   lambda: Numeric with value in R+, latency rate
  #   gamma: Numeric with value in R+, recovery rate
  #   eta: Numeric with value in R+, fatality rate
  #
  # Returns:
  #   The Jacobian matrix of the ODE function
  
  f1 <- function(X0) f(model = model, X0, beta = beta, pop = pop, lambda = lambda, gamma = gamma, eta = eta)
  f2 <- function(beta) f(model = model, X0 = X0, beta, pop = pop, lambda = lambda, gamma = gamma, eta = eta)
  out <- cbind(jacobian(f1, x = X0), jacobian(f2, x = beta))
  return(out)
}

#' Measurement model
#'
#' This function evaluates the measurement model for a given set of parameters.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param X Numeric vector representing the solution of the ODE and its 2 first derivatives.
#' @param U Numeric vector of size 2 representing the latent parameter of the ODE and its first derivative.
#' @param pop Integer representing the total population.
#' @param lambda Numeric value representing the latency rate.
#' @param gamma Numeric representing the recovery rate.
#' @param eta Numeric representing the fatality rate.
#' @return Numeric vector representing the evaluation of the measurement model.
#' @export
h <- function(model, X, U, pop, lambda, gamma, eta = 0) {
  # Arguments:
  #   model: String, type of model to be used ('SEIRD' or 'SEIR' available)
  #   X: Numeric vector, solution of the ODE and its 2 first derivatives
  #   U: Numeric vector of size 2, latent parameter of the ODE and its first derivative
  #   pop: Integer, total population
  #   lambda: Numeric with value in R+, latency rate
  #   gamma: Numeric with value in R+, recovery rate
  #   eta: Numeric with value in R+, fatality rate
  #
  # Returns:
  #   sol: Numeric vector representing the evaluation of the measurement model
  
  beta <- U[1]
  beta <- sigmoid(beta) # rescaling of beta to [0,5]
  
  if (model == 'SEIRD') {
    X0 <- X[1:5]
    X1 <- X[6:10]
  }
  
  if (model == 'SEIR') {
    X0 <- X[1:4]
    X1 <- X[5:8]
  }

  ODE <- f(model = model,
           X0 = X0, 
           beta = beta, 
           pop = pop, 
           lambda = lambda, 
           gamma = gamma, 
           eta = eta)
  
  sol <- c(X1 - ODE)
  
  return(sol)
}

#' Jacobian of measurement model using jacobian function
#'
#' This function calculates the Jacobian matrix of the measurement model using the jacobian function.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param X Vector representing the solution of the ODE and its 2 first derivatives.
#' @param U Vector of length 2 representing the latent parameter of the ODE and its first derivative.
#' @param pop Integer representing the total population.
#' @param lambda Numeric value representing the latency rate.
#' @param gamma Numeric representing the recovery rate.
#' @param eta Numeric representing the fatality rate.
#' @return The Jacobian matrix of the measurement model.
#' @export
jacobian_h <- function(model, X, U, pop, lambda, gamma, eta = 0) {
  # Arguments:
  #   model: String, type of model to be used ('SEIRD' or 'SEIR' available)
  #   X: Vector representing the solution of the ODE and its 2 first derivatives
  #   U: Vector of length 2 representing the latent parameter of the ODE and its first derivative
  #   pop: Integer, total population
  #   lambda: Numeric, latency rate
  #   gamma: Numeric, recovery rate
  #   eta: Numeric, fatality rate
  #
  # Returns:
  #   The Jacobian matrix of the measurement model
  
  h1 <- function(X) h(model = model, X = X, U = U, pop = pop, lambda = lambda, gamma = gamma, eta = eta)
  h2 <- function(U) h(model = model, X = X, U, pop = pop, lambda = lambda, gamma = gamma, eta = eta)
  out <- cbind(jacobian(h1, x = X), jacobian(h2, x = U))
  return(out)
}

#' Transition matrices calculation
#'
#' Calculates the transition matrices for a given drift, dispersion, Wiener noise, and time steps.
#'
#' @param drift Numeric matrix representing the drift.
#' @param dispersion Numeric matrix representing the dispersion.
#' @param noise_wiener Numeric matrix representing the noise of the Wiener process.
#' @param steps Numeric, time steps.
#' @return A list containing the transition matrices A and Q.
#' @export
transition_matrices <- function(drift, dispersion, noise_wiener, steps){
  # Arguments:
  #   drift: Numeric matrix representing the drift.
  #   dispersion: Numeric matrix representing the dispersion.
  #   noise_wiener: Numeric matrix representing the noise of the Wiener process.
  #   steps: Numeric, time steps.
  #
  # Returns:
  #   A list containing the transition matrices A and Q.
  
  l <- nrow(drift)
  
  top_row <- cbind(steps*drift, steps*dispersion %*% noise_wiener %*% t(dispersion))
  bottom_row <- cbind(matrix(0, nrow = l, ncol = l), -steps*t(drift))
  Gamma <- rbind(top_row, bottom_row)
  exp_Gamma <- expm(Gamma)
  M_1 <- exp_Gamma[1:l, 1:l]
  M_2 <- exp_Gamma[1:l, (l + 1):(2 * l)]
  M_3 <- exp_Gamma[(l + 1):(2 * l), (l + 1):(2 * l)]
  
  A <- expm(steps*drift)
  Q <- M_2 %*% t(M_1)
  
  res <- list(A,Q)
  
  return(res)
}

#' Prediction function for latent parameter U or state variable X
#'
#' This function calculates the predicted mean and covariance of either U or X using the given parameters.
#'
#' @param m Numeric vector, the (previous) mean of U or X.
#' @param P Numeric matrix, the (previous) covariance of U or X.
#' @param A Numeric matrix, transition matrix of U or X.
#' @param Q Numeric matrix, transition noise of U or X.
#' @return A list containing the predicted mean and covariance.
#' @export
prediction <- function(m, P, A, Q) {
  # Arguments:
  #   m: Numeric vector, the (previous) mean of U or X.
  #   P: Numeric matrix, the (previous) covariance of U or X.
  #   A: Numeric matrix, transition matrix of U or X.
  #   Q: Numeric matrix, transition noise of U or X.
  #
  # Returns:
  #   A list containing the predicted mean and covariance.

  m_out <- as.vector(A %*% m)
  P_out <- as.matrix(A %*% P %*% t(A) + Q)
  
  out <- list(m_out, P_out)
  
  return(out)
}

#' Update step for observations (Kalman Filter)
#'
#' This function updates the mean and covariance of (X,U) using observations.
#'
#' @param m Numeric vector, predicted mean of (X,U).
#' @param P Numeric matrix, predicted covariance of (X,U).
#' @param y Numeric vector, observations.
#' @param H Numeric matrix, observation matrix.
#' @param R Numeric matrix, observation noise.
#' @return A list containing the updated mean and covariance of U and X.
#' @export
update_of_observations <- function(m, P, y, H, R) {
  # Arguments:
  #   m: Numeric vector, predicted mean of (X,U).
  #   P: Numeric matrix, predicted covariance of (X,U).
  #   y: Numeric vector, observations.
  #   H: Numeric matrix, observation matrix.
  #   R: Numeric matrix, observation noise.
  #
  # Returns:
  #   A list containing the updated mean and covariance of U and X.
  
  v <- y - H %*% m # residual
  S <- H %*% P %*% t(H) + R # innovation covariance
  S_inv <- svd.inverse(S)
  # S_inv <- solve(S)
  K <- P %*% t(H) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  cond_S_obs <- kappa(S) # condition number
  
  out <- list(m_out = m_out, P_out = P_out, cond_S_obs = cond_S_obs)
  
  return(out)
}

#' Update step for states (Extended Kalman Filter)
#'
#' This function updates the mean and covariance of (X,U) using the measurement model and its Jacobian.
#'
#' @param m Numeric vector, predicted mean of (X,U).
#' @param P Numeric matrix, predicted covariance of (X,U).
#' @param h Numeric vector, measurement model.
#' @param J Numeric matrix, Jacobian of the measurement model.
#' @param jit Numeric, value to add as jitter to the innovation covariance.
#' @return A list containing the updated mean and covariance of U and X.
#' @export
update_of_states <- function(m, P, h, J, jit) {
  # Arguments:
  #   m: Numeric vector, predicted mean of (X,U).
  #   P: Numeric matrix, predicted covariance of (X,U).
  #   h: Numeric vector, measurement model.
  #   J: Numeric matrix, Jacobian of the measurement model.
  #   jit: Numeric, value to add as jitter to the innovation covariance.
  #
  # Returns:
  #   A list containing the updated mean and covariance of U and X.
  
  if (is.null(jit)) {
    jit_mat <- diag(0, nrow = nrow(J))
  } else {
    jit_mat <- diag(jit, nrow = nrow(J))
  }
  
  v <- - h # residual
  S <- J %*% P %*% t(J) + jit_mat # innovation covariance
  S_inv <- svd.inverse(S)
  # S_inv <- solve(S)
  K <- P %*% t(J) %*% S_inv # Kalman gain
  m_out <- m + K %*% v # updated mean
  P_out <- P - K %*% S %*% t(K) # updated covariance
  
  cond_S_ode <- kappa(S) # condition number
  
  return(list(m_out = m_out, P_out = P_out, cond_S_ode = cond_S_ode))
}

#' Combine covariance matrices P_X and P_U
#'
#' This function combines the predicted covariance matrices of X and U into a single covariance matrix.
#'
#' @param P_X Numeric matrix, predicted covariance of X.
#' @param P_U Numeric matrix, predicted covariance of U.
#' @return A single matrix representing the combined covariance of X and U.
#' @export
matrix_P <- function(P_X, P_U) {
  # Arguments:
  #   P_X: Numeric matrix, predicted covariance of X.
  #   P_U: Numeric matrix, predicted covariance of U.
  #
  # Returns:
  #   A single matrix representing the combined covariance of X and U.
  
  P_U <- as.matrix(P_U)
  P_X <- as.matrix(P_X)
  
  ncol_P <- ncol(P_X) + ncol(P_U)
  nrow_P <- nrow(P_X) + nrow(P_U)
  P <- matrix(0, nrow = nrow_P, ncol = ncol_P)
  P[1:nrow(P_X), 1:ncol(P_X)] <- P_X  # Top left corner
  P[(nrow_P - nrow(P_U) + 1):nrow_P, (ncol_P - ncol(P_U) + 1):ncol_P] <- P_U  # Bottom right corner
  
  return(P)
}

#' Perform inference using a state-space model
#'
#' This function performs inference over a specified time grid using a state-space model.
#' It iteratively updates the state estimates and covariances based on observations and predictions.
#'
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
#' @param grids List of time grids for inference.
#' @param obs Data frame, contains the dates and the compartment counts.
#' @param jit Numeric, value to add as jitter to the innovation covariance (optional, default 0.001).
#' @param jit_X Numeric, value to add as jitter to the covariance matrix of X (optional, default 1e-9).
#' @param jit_U Numeric, value to add as jitter to the covariance matrix of U (optional, default 1e-9).
#' @param initial_params List of initial parameters obtained from the initialization function.
#'
#' @return A list with the inferred values of X, U, P_X, and P_U.
#' @export
inference <- function(model, grids, obs, 
                      jit = 0.001, 
                      jit_X = 1e-9, jit_U = 1e-9,
                      initial_params){
  
  # Arguments:
  #   model: String, type of model to be used ('SEIRD' or 'SEIR' available).
  #   grids: List of time grids for inference.
  #   obs: Data frame, contains the dates and the compartment counts.
  #   jit: Numeric, value to add as jitter to the innovation covariance (optional, default 0.001).
  #   jit_X: Numeric, value to add as jitter to the covariance matrix of X (optional, default 1e-9).
  #   jit_U: Numeric, value to add as jitter to the covariance matrix of U (optional, default 1e-9).
  #   initial_params: List of initial parameters obtained from the initialization function.
  # 
  # Returns:
  #   A list with the inferred values of X, U, P_X and P_U.
  
  time_grid <- grids$time_grid
  data_grid <- grids$data_grid
  ode_grid <- grids$ode_grid
  steps <- grids$steps
  
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
  lambda <- initial_params$lambda
  gamma <- initial_params$gamma
  eta <- initial_params$eta
  
  # Transition matrices and noise
  trans_X <- transition_matrices(drift = F_X, dispersion = L_X, noise_wiener = noise_wiener_X, steps = steps) 
  A_X <- as.matrix(trans_X[[1]])
  Q_X <- as.matrix(trans_X[[2]])
  trans_U <- transition_matrices(drift = F_U, dispersion = L_U, noise_wiener = noise_wiener_U, steps = steps) 
  A_U <- as.matrix(trans_U[[1]])
  Q_U <- as.matrix(trans_U[[2]])
  
  # Create lists to store inferred values and condition numbers
  X_values <- matrix(data = NA, nrow = length(X), ncol = length(time_grid))
  U_values <- matrix(data = NA, nrow = length(U), ncol = length(time_grid))
  P_X_values <- array(data = NA, dim = c(length(X), length(X), length(time_grid)))
  P_U_values <- array(data = NA, dim = c(length(U), length(U), length(time_grid)))
  cond_numbers <- matrix(data = NA, nrow = 4, ncol = (length(time_grid)-1))
  
  # Initialize these lists
  X_values[, 1] <- X
  P_X_values[, , 1] <- P_X
  U_values[, 1] <- U
  P_U_values[, , 1] <- P_U

  # Get vector indicating which compartments have been observed
  get_observation_vector <- function(obs_data) {
    ind <- c()
    if ("S" %in% colnames(obs_data)) {
      ind <- c(ind,1)
    }
    if ("E" %in% colnames(obs_data)) {
      ind <- c(ind,2)
    }
    if ("I" %in% colnames(obs_data)) {
      ind <- c(ind,3)
    }
    if ("R" %in% colnames(obs_data)) {
      ind <- c(ind,4)
    }
    if ("D" %in% colnames(obs_data)) {
      ind <- c(ind,5)
    }
    return(ind)
  }
  ind <- get_observation_vector(obs)
  
  # Observations without time
  obs_without_time <- obs[-1]
  
  # Iterate over the time grid
  for (i in 2:length(time_grid)) {
    
    # Get previous values and covariances
    X_prev <- as.matrix(X_values[,i-1])
    P_X_prev <- as.matrix(P_X_values[,,i-1])
    U_prev <- as.matrix(U_values[,i-1])
    P_U_prev <- as.matrix(P_U_values[,,i-1])
    
    # Prediction step
    prediction_X <- prediction(m = X_prev, P = P_X_prev, A = A_X, Q = Q_X)
    X_pred <- as.vector(prediction_X[[1]])
    P_X_pred <- as.matrix(prediction_X[[2]])
    prediction_U <- prediction(m = U_prev, P = P_U_prev, A = A_U, Q = Q_U)
    U_pred <- as.vector(prediction_U[[1]])
    P_U_pred <- as.matrix(prediction_U[[2]])
    
    # Update of observations if available at the current time point
    if (any(data_grid == time_grid[i])) {
      m <- c(X_pred, U_pred)
      P <- matrix_P(P_X_pred, P_U_pred)
      y <- unlist(obs_without_time[i,])
      updated_obs <- update_of_observations(m = m, P = P, y = y, H = H, R = R)
      
      if (model == 'SEIRD') {
        X_pred <- as.vector(updated_obs[[1]][1:15])
        P_X_pred <- as.matrix(updated_obs[[2]][1:15, 1:15]) 
        U_pred <- as.vector(updated_obs[[1]][16:17])
        P_U_pred <- as.matrix(updated_obs[[2]][16:17, 16:17])
      }
      
      if (model == 'SEIR') {
        X_pred <- as.vector(updated_obs[[1]][1:12])
        P_X_pred <- as.matrix(updated_obs[[2]][1:12, 1:12])
        U_pred <- as.vector(updated_obs[[1]][13:14])
        P_U_pred <- as.matrix(updated_obs[[2]][13:14, 13:14])
      }
    }
    
    # Condition number
    cond_S_obs <- updated_obs[[3]]
    
    # Update of states if required at the current time point
    if (any(ode_grid == time_grid[i])) {
      m <- c(X_pred, U_pred)
      P <- matrix_P(P_X_pred, P_U_pred)
      
      h_val <- h(model, X_pred, U_pred, pop, lambda, gamma, eta)
      J_val <- jacobian_h(model, X_pred, U_pred, pop, lambda, gamma, eta)
      
      updated_states <- update_of_states(m = m, P = P, h = h_val, J = J_val, jit = jit)
      
      if (model == 'SEIRD') {
        X_pred <- as.vector(updated_states[[1]][1:15])
        P_X_pred <- as.matrix(updated_states[[2]][1:15, 1:15]) 
        U_pred <- as.vector(updated_states[[1]][16:17])
        P_U_pred <- as.matrix(updated_states[[2]][16:17, 16:17])
      }
      
      if (model == 'SEIR') {
        X_pred <- as.vector(updated_states[[1]][1:12])
        P_X_pred <- as.matrix(updated_states[[2]][1:12, 1:12])
        U_pred <- as.vector(updated_states[[1]][13:14])
        P_U_pred <- as.matrix(updated_states[[2]][13:14, 13:14])
      }
    }

    # Condition number
    cond_S_ode <- updated_states[[3]]
    
    # Function to check if a matrix is positive semidefinite
    is_positive_semidefinite <- function(matrix, tol = 1e-10) {
      eigenvalues <- eigen(matrix)$values
      eigenvalues_real <- Re(eigenvalues)
      all(eigenvalues_real >= -tol)
    }

    # Initial jitter values
    jitter_for_X <- jit_X  # Initial jitter value for P_X_pred
    jiter_for_U <- jit_U  # Initial jitter value for P_U_pred
    
    # Ensure positive semidefiniteness of P_X_pred
    while (!is_positive_semidefinite(P_X_pred)) {
      P_X_pred <- P_X_pred + diag(jitter_for_X, nrow = nrow(P_X_pred))
      jitter_for_X <- jitter_for_X * 10  # Increase jitter value
    }
    
    # Ensure positive semidefiniteness of P_U_pred
    while (!is_positive_semidefinite(P_U_pred)) {
      P_U_pred <- P_U_pred + diag(jiter_for_U, nrow = nrow(P_U_pred))
      jiter_for_U <- jiter_for_U * 10  # Increase jitter value
    }
    
    # Compute condition numbers
    cond_P_X <- kappa(P_X_pred)
    cond_P_U <- kappa(P_U_pred)
    
    # Store current values and covariance matrices
    X_values[, i] <- X_pred
    U_values[, i] <- U_pred
    P_X_values[, , i] <- P_X_pred
    P_U_values[, , i] <- P_U_pred
    cond_numbers[, (i-1)] <- c(cond_S_obs, cond_S_ode, cond_P_X, cond_P_U)
    
  }
  
  # Return inferred values
  list(X_values = X_values, U_values = U_values, 
       P_X_values = P_X_values, P_U_values = P_U_values, 
       cond_numbers = cond_numbers)
}
