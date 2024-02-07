# PROBODER - Initialization file #
# Claire Descombes               #

library(Matrix)

#' @param lengthscale
#' @return 
#' 

# Drift matrices.
F_U <- matrix(c(0,-(sqrt(3)/lengthscale)^2,1,-2*sqrt(3)/lengthscale), nrow = 2, ncol = 2)
F_X <- sparseMatrix(i = 1:8, j = 5:12, x = 1, dims = c(12,12))
  
# Dispersion matrices.
L_U <- matrix(c(0,1), nrow = 2, ncol = 1)
L_X <- sparseMatrix(i = 9:12, j = 1:4, x = 1, dims = c(12,4))

# Prediction step.


# Update step on tau_OBS.


# Update step on tau_ODE.
