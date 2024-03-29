#' Initialization Function for ProbODER
#' 
#' Initializes parameters and matrices for the given epidemiological model.
#' 
#' @param model Model type ('SIR' or 'SID').
#' @param obs Data frame with dates and observations of SIR or SID counts.
#' @param beta0 Initial value for contact rate in (0,1).
#' @param beta0prime Initial value for 1st derivative of contact rate.
#' @param gamma Recovery rate.
#' @param eta Fatality rate.
#' @param l Length scale.
#' @param noise_wiener Noise of the Wiener processes.
#' @param pop Population.
#' @return A list containing initialized parameters and matrices.
#' @export
initialization <- function(model, obs, beta0, beta0prime, gamma = 0, eta = 0, l = 1, noise_wiener = 0.1, pop){

  # Initialize latent parameter (contact rate) and its first derivative
  U <- as.vector(c(
    logit(beta0), # Initial value for contact rate rescaled to the real line
    beta0prime) # Initial value for 1st derivative
  )
  
  # Initialize solution of SIRD-ODE and its two first derivatives
  if(model == 'SIR'){
    X0 <- c(obs[1,2], obs[1,3], obs[1,4], 0) # Initial values for S, I, R and D
  }else if(model == 'SID'){
    X0 <- c(obs[1,2], obs[1,3], 0, obs[1,4]) # Initial values for S, I, R and D
  }else{
    print('Wrong model type!')
  }
  X1 <- f(X0, beta0, pop, gamma, eta) # Initial values for 1st derivatives
  X2 <- diag(jacobian_f(X0, beta0, pop, gamma, eta)) # Initial values for 2nd derivatives
  X <- as.vector(c(X0,X1,X2))
  
  # Drift matrices
  F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
  F_X <- as.matrix(sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12)))
  
  # Dispersion matrices
  L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
  L_X <- as.matrix(sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4)))
  
  # Observation matrix (for observation of S,I, and R or D)
  if(model == 'SID'){
    H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,4), x = 1, dims = c(3,14)))
  }else if(model == 'SIR'){
    H <- as.matrix(sparseMatrix(i = c(1,2,3), j = c(1,2,3), x = 1, dims = c(3,14)))
  }else{
    print('Wrong model type!')
  }
  
  # Observation noise
  R <- cov(obs[,-1])
  
  # Noise of priors
  P_U <- matrix(0, nrow = 2, ncol = 2)
  P_X <- matrix(0, nrow = 12, ncol = 12)
  
  # Noise of Wiener process
  noise_wiener_U <- diag(noise_wiener, nrow = ncol(L_U), ncol = ncol(L_U))
  noise_wiener_X <- diag(noise_wiener, nrow = ncol(L_X), ncol = ncol(L_X))
  
  out <- list(X = X, U = U, P_X = P_X, P_U = P_U, 
              F_X = F_X, F_U = F_U, L_X = L_X, L_U = L_U, 
              noise_wiener_X = noise_wiener_X, noise_wiener_U = noise_wiener_U,
              R = R, H = H, pop = pop, gamma = gamma, eta = eta, l = l)
  
  return(out)
}
