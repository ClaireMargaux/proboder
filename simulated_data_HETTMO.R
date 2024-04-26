#############################################
####### SIMULATED DATA USING HETTMO #########
#############################################

directory <- "~/Documents/GitHub/proboder/Data/HETTMO" # directory for data

library(HETTMO)

params_unstratified = set_parameters()
params_unstratified$p_detect1 <- 1
params_unstratified$p_detect2 <- 1

population_simulate <- params_unstratified$popsize

# Time grid
steps <- 1
grid <- seq(1,45,by=steps)
n <- length(grid)

simulate <- simulate_data(params = params_unstratified, ts = grid)
PCR <- simulate
reproduction_number <- simulate[[2]]
real_beta <- simulate[[3]]
I <- simulate[[4]]
E <- simulate[[5]]
R <- simulate[[6]]
sero <- simulate[[7]]
observations_simulate <- data.frame(t = grid, E = E, I = I, R = R)
HETTMO_data_simulated <- list(PCR,reproduction_number,real_beta,I,E,R,sero)

df_params <- data.frame(pop = population_simulate, lambda = 0.3703704, gamma = 0.3703704, eta = 0)

################################
######## SAVE THE DATA #########
################################

save(observations_simulate, file = file.path(directory, "simulated_data_HETTMO.Rdata"))
save(population_simulate, file = file.path(directory, "simulated_params_HETTMO.Rdata"))
saveRDS(real_beta, file = file.path(directory, "simulated_beta_HETTMO.Rds"))
save(HETTMO_data_simulated, file = file.path(directory, "HETTMO_save.Rdata"))