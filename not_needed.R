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