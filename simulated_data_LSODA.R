############################################
####### SIMULATED DATA USING LSODA #########
############################################

directory <- "~/Documents/GitHub/proboder/Data/LSODA" # directory for data
noise <- 1 # noise to add on the data

# Time grid
steps <- 0.2
grid <- seq(0,10,by=steps)
n <- length(grid)

# Parameters
lambda <- 0.6
gamma <- 0.4
eta <- 0.2
pop <- 1000
parms  <- c(lambda = lambda, gamma = gamma, eta = eta, pop = pop)

# Beta function: sine
# beta <- function(t){0.3*sin(t)+0.5}
# beta_val <- beta(grid)
# beta0 <- beta_val[1]

# Beta function: log
beta <- function(t){0.3*log(t+1)+0.1}
beta_val <- beta(grid)
beta0 <- beta_val[1]

# Start values for steady state
y <- xstart <- c(S = 984, E = 10, I = 2, R = 2, D = 2)

# ODE solver
library(deSolve)
SPCmod <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    beta <- beta(t)
    dS <- - beta * S * I / pop
    dE <- beta * S * I / pop - lambda * E
    dI <- lambda * E - gamma * I - eta * I
    dR <- gamma * I
    dD <- eta * I
    res <- c(dS, dE, dI, dR, dD)
    list(res)
  })
}

# Solving
out <-  lsoda(xstart, grid, SPCmod, parms)

# Data frame
df <- data.frame(t <- grid, S <- out[,"S"], E <- out[,"E"], I <- out[,"I"], R <- out[,"R"], D <- out[,"D"])
colnames(df) <- c("t","S","E","I","R","D")
df_beta <- data.frame(t <- grid, beta_val <- beta_val)

# Adding noise
y <- df[,c("S","E","I","R","D")]
y_with_noise <- y + matrix(rnorm(n*n, mean = 0, sd = sqrt(noise)), nrow = n)
df_with_noise <- data.frame(t <- grid, y_with_noise)
colnames(df_with_noise) <- c("t","S","E","I","R","D")

# Plotting
countsplot <- ggplot(df, aes(x = t, y = S)) +
  geom_line(aes(color = 'S')) +
  geom_line(aes(x = t, y = E, color = 'E')) +
  geom_line(aes(x = t, y = I, color = 'I')) +
  geom_line(aes(x = t, y = R, color = 'R')) +
  geom_line(aes(x = t, y = D, color = 'D')) +
  theme_minimal() +
  scale_y_continuous(trans='log10') + 
  ylab("log(count)") +
  xlab("time") + 
  ggtitle("LSODA") 

countsplot

betaplot <- ggplot(df_beta, aes(x = t, y = beta_val)) +
  geom_line(aes(color = 'beta')) +
  theme_minimal() +
  ylab("beta") +
  xlab("time") +
  ylim(c(0,1)) +
  ggtitle("Simulated beta") 

betaplot

df_params <- data.frame(pop = pop, lambda = lambda, gamma = gamma, eta = eta, obs_noise = noise)

################################
######## SAVE THE DATA #########
################################

save(df, file = file.path(directory, "simulated_data_LSODA.Rdata"))
save(df_with_noise, file = file.path(directory, "simulated_noisy_data_LSODA.Rdata"))
save(df_params, file = file.path(directory, "simulated_params_LSODA.Rdata"))
saveRDS(beta_val, file = file.path(directory, "simulated_beta_LSODA.Rds"))
save(df_beta, file = file.path(directory, "simulated_beta_LSODA.Rdata"))