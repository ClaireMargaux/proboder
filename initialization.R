#' Initialization Function for ProbODER
#' 
#' Initializes parameters and matrices for the given epidemiological model.
#' 
#' @param obs Data frame with dates and observations of SEIRD counts.
#' @param beta0 Initial value for contact rate in (0,1).
#' @param beta0prime Initial value for 1st derivative of contact rate.
#' @param lambda Latency rate.
#' @param gamma Recovery rate.
#' @param eta Fatality rate.
#' @param l Length scale.
#' @param noise_obs Noise of the observations.
#' @param noise_wiener_X Noise of the Wiener process modelling X.
#' @param noise_wiener_U Noise of the Wiener process modelling U.
#' @param pop Population.
#' @return A list containing initialized parameters and matrices.
#' @export
initialization <- function(obs, beta0, beta0prime, 
                           lambda = 0.6, gamma = 0.4, eta = 0.2, 
                           l = 1, scale = 1, noise_obs = 0.1,
                           noise_wiener_X = 1, noise_wiener_U = 1,
                           pop){
  
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
  
  # Initialize compartment counts and their two first derivatives
  X0 <- numeric(5)
  X0[1] <- ifelse(1 %in% ind, obs[1, which(colnames(obs) == "S")], pop)
  X0[2] <- ifelse(2 %in% ind, obs[1, which(colnames(obs) == "E")], 0)
  X0[3] <- ifelse(3 %in% ind, obs[1, which(colnames(obs) == "I")], 0)
  X0[4] <- ifelse(4 %in% ind, obs[1, which(colnames(obs) == "R")], 0)
  X0[5] <- ifelse(5 %in% ind, obs[1, which(colnames(obs) == "D")], 0)
  X1 <- f(X0, beta0, pop, lambda, gamma, eta) # Initial values for 1st derivatives
  X2 <- diag(jacobian_f(X0, beta0, pop, lambda, gamma, eta)) # Initial values for 2nd derivatives
  X <- as.vector(c(X0,X1,X2))
  
  # Initialize latent parameter (contact rate) and its first derivative
  U <- as.vector(c(
    logit(beta0), # Initial value for contact rate rescaled to the real line
    beta0prime) # Initial value for 1st derivative
  )
  
  # Drift matrices
  F_X <-  as.matrix(sparseMatrix(i = 1:10, j = 6:15, x = 1, dims = c(15,15)))
  F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
  
  # Dispersion matrices
  L_X <-  as.matrix(sparseMatrix(i = 11:15, j = 1:5, x = 1, dims = c(15,5)))
  L_U <- matrix(c(0,scale), nrow = 2, ncol = 1)
  
  # Observation matrix
  H <- as.matrix(sparseMatrix(i = 1:length(ind), j = ind, x = 1, dims = c(length(ind),17)))
  
  # Observation noise
  R <- diag(noise_obs, nrow = length(ind), ncol = length(ind))
  #R <- cov(obs[,-1])
  
  # Noise of priors
  P_X <- matrix(0, nrow = 15, ncol = 15)
  P_U <- matrix(0, nrow = 2, ncol = 2)

  # Noise of Wiener process
  noise_wiener_X <- diag(noise_wiener_X, nrow = ncol(L_X), ncol = ncol(L_X))
  noise_wiener_U <- diag(noise_wiener_U, nrow = ncol(L_U), ncol = ncol(L_U))
  
  out <- list(X = X, U = U, P_X = P_X, P_U = P_U, 
              F_X = F_X, F_U = F_U, L_X = L_X, L_U = L_U, 
              noise_wiener_X = noise_wiener_X, noise_wiener_U = noise_wiener_U,
              R = R, H = H, pop = pop, 
              lambda = lambda, gamma = gamma, eta = eta, l = l)

  return(out)
}
