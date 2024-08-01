#' Load Required Packages
#'
#' This function checks if a set of required R packages is installed and loads them if available. 
#' If a package is not found, it attempts to install it from CRAN before loading.
#'
#' @return NULL
#' @export
load_required_packages <- function() {
  # List of packages to load
  required_packages <- c(
    "Matrix",
    "numDeriv",
    "matrixcalc",
    "knitr",
    "kableExtra",
    "ggplot2",
    "gridExtra",
    "tictoc",
    "progress",
    "dplyr",
    "lubridate"
  )
  
  # Load each package if not already loaded
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Initialization Function for ProbODER
#' 
#' Initializes parameters and matrices for the given epidemiological model.
#' 
#' @param model String, type of model to be used ('SEIRD' or 'SEIR' available).
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
#' @param num_initial_values Number of values over which we take the mean to initialize X.
#' @param start_S Starting count of S compartment (optional).
#' @param start_E Starting count of E compartment (optional).
#' @param start_I Starting count of I compartment (optional).
#' @param start_R Starting count of R compartment (optional).
#' @param start_D Starting count of D compartment (optional).
#' @return A list containing initialized parameters and matrices.
#' @export
initialization <- function(model, obs, 
                           beta0, beta0prime, 
                           lambda, gamma, eta = 0, 
                           l, scale = 1, noise_obs,
                           noise_X, noise_U,
                           noise_wiener_X, noise_wiener_U,
                           pop,
                           num_initial_values = 1,
                           start_S = NULL, start_E = NULL, 
                           start_I = NULL, start_R = NULL, 
                           start_D = NULL){
  
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
  
  if (!model %in% c("SEIRD", "SEIR")) {
    warning("The model should be either 'SEIRD' or 'SEIR'.")
  }
  
  if (model == 'SEIRD') {
  
    # Initialize compartment counts and their two first derivatives
    X0 <- numeric(5)
    X0[1] <- if (!is.null(start_S)) {
      start_S
    } else {
      ifelse("S" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "S")]), pop)
    }
    X0[2] <- if (!is.null(start_E)) {
      start_E
    } else {
      ifelse("E" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "E")]), 1)
    }
    X0[3] <- if (!is.null(start_I)) {
      start_I
    } else {
      ifelse("I" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "I")]), 0)
    }
    X0[4] <- if (!is.null(start_R)) {
      start_R
    } else {
      ifelse("R" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "R")]), 0)
    }
    X0[5] <- if (!is.null(start_D)) {
      start_D
    } else {
      ifelse("D" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "D")]), 0)
    }
    X0 <- unlist(X0)
    X1 <- f(model, X0, beta0, pop, lambda, gamma, eta) # Initial values for 1st derivatives
    X2 <- diag(jacobian_f(model, X0, beta0, pop, lambda, gamma, eta)) # Initial values for 2nd derivatives
    X <- as.vector(c(X0,X1,X2))
  
    # Drift and dispersion matrices
    F_X <-  as.matrix(sparseMatrix(i = 1:10, j = 6:15, x = 1, dims = c(15,15)))
    L_X <-  as.matrix(sparseMatrix(i = 11:15, j = 1:5, x = 1, dims = c(15,5)))
    
    # Observation matrix
    H <- as.matrix(sparseMatrix(i = 1:length(ind), j = ind, x = 1, dims = c(length(ind),17)))
    
    # Noise of prior
    P_X <- matrix(noise_X, nrow = 15, ncol = 15)
  }
  
  if (model == 'SEIR') {
    
    # Initialize compartment counts and their two first derivatives
    X0 <- numeric(4)
    X0[1] <- if (!is.null(start_S)) {
      start_S
    } else {
      ifelse("S" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "S")]), pop)
    }
    X0[2] <- if (!is.null(start_E)) {
      start_E
    } else {
      ifelse("E" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "E")]), 1)
    }
    X0[3] <- if (!is.null(start_I)) {
      start_I
    } else {
      ifelse("I" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "I")]), 0)
    }
    X0[4] <- if (!is.null(start_R)) {
      start_R
    } else {
      ifelse("R" %in% colnames(obs), mean(obs[1:num_initial_values, which(colnames(obs) == "R")]), 0)
    }
    X0 <- unlist(X0)
    X1 <- f(model, X0, beta0, pop, lambda, gamma, eta) # Initial values for 1st derivatives
    X2 <- diag(jacobian_f(model, X0, beta0, pop, lambda, gamma, eta)) # Initial values for 2nd derivatives
    X <- as.vector(c(X0,X1,X2))
    
    # Drift and dispersion matrices
    F_X <-  as.matrix(sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12)))
    L_X <-  as.matrix(sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4)))
    
    # Observation matrix
    H <- as.matrix(sparseMatrix(i = 1:length(ind), j = ind, x = 1, dims = c(length(ind),14)))
    
    # Noise of prior
    P_X <- matrix(noise_X, nrow = 12, ncol = 12)
    
    # Fatality rate
    eta <- 0
  }
  
  # Initialize latent parameter (contact rate) and its first derivative
  U <- as.vector(c(
    logit(beta0), # Initial value for contact rate rescaled to the real line
    beta0prime) # Initial value for 1st derivative
  )
  
  # Drift and dispersion matrices
  F_U <- matrix(c(0,-(sqrt(3)/l)^2,1,-2*sqrt(3)/l), nrow = 2, ncol = 2)
  L_U <- matrix(c(0,scale), nrow = 2, ncol = 1)
  
  # Observation noise
  if(length(noise_obs) == 1){
    R <- diag(noise_obs, nrow = length(ind), ncol = length(ind))
  }else{
    R <- noise_obs
  }
  
  # Noise of prior
  P_U <- matrix(noise_U, nrow = 2, ncol = 2)

  # Noise of Wiener processes
  noise_wiener_X <- diag(noise_wiener_X, nrow = ncol(L_X), ncol = ncol(L_X))
  noise_wiener_U <- diag(noise_wiener_U, nrow = ncol(L_U), ncol = ncol(L_U))
  
  out <- list(X = X, U = U, P_X = P_X, P_U = P_U, 
              F_X = F_X, F_U = F_U, L_X = L_X, L_U = L_U, 
              noise_wiener_X = noise_wiener_X, noise_wiener_U = noise_wiener_U,
              R = R, H = H, pop = pop, 
              lambda = lambda, gamma = gamma, eta = eta, l = l)

  return(out)
}

#' Generate Data Grid for inference
#' 
#' Generate a data grid with optional additional points between observations for ODE evaluation.
#' 
#' @param obs Data frame containing the observed compartments.
#' @param num_points_between Integer specifying the number of points to add between observations.
#' @param daily_or_weekly To define the time steps (optional).
#' 
#' @return A list containing the overall time grid, data grid, ODE grid, and time steps.
#' @export
generate_grid <- function(obs, num_points_between, daily_or_weekly = NULL) {
  
  # Extract the observed time points
  data_grid <- obs$t
  
  # Initialize the ODE grid with observed time points
  ode_grid <- data_grid
  
  if (num_points_between > 0) {
    # Generate equidistant points between observed time points
    for (i in 1:(length(data_grid) - 1)) {
      equidistant_points <- seq(
        from = data_grid[i], 
        to = data_grid[i + 1], 
        length.out = num_points_between + 2
      )[-c(1, num_points_between + 2)]
      
      # Append the equidistant points to the ODE grid and sort
      ode_grid <- sort(c(ode_grid, equidistant_points))
    }
  }
  
  # Combine data and ODE grids, ensuring unique and sorted time points
  time_grid <- sort(unique(c(data_grid, ode_grid)))
  
  # Calculate time steps
  if(is.null(daily_or_weekly)){
    steps <- ode_grid[2] - ode_grid[1]
  } else {
    if(daily_or_weekly == 'daily'){
      steps = 1
    } else {
      steps = 7
    }
  }
  
  # Return the grids and time steps as a list
  list(
    time_grid = time_grid, 
    data_grid = data_grid, 
    ode_grid = ode_grid, 
    steps = steps
  )
}
