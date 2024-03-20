#################################### PROBODER ##################################
################################ Claire Descombes ##############################
###################################### Data ####################################

directory <- "~/Documents/GitHub/proboder/Data" # directory for data

################################
####### SIMULATED DATA #########
################################

library(HETTMO)
params_unstratified = set_parameters()
params_unstratified$p_detect1 <- 1
params_unstratified$p_detect2 <- 1
population_simulate <- params_unstratified$popsize
simulate <- simulate_data(params = params_unstratified, ts = 1:45)
I <- simulate[[1]]
S <- rep(0,length(I))
for (i in 1:length(I)){
  S[i] <- population_simulate - sum(I[1:i])
}
R <- simulate[[6]]
real_beta <- simulate[[3]]
date <- 1:length(S)
observations_simulate <- data.frame(date,S,I,R)

################################
######## SAVE THE DATA #########
################################

save(observations_simulate, file = file.path(directory, "simulated_data.Rdata"))
saveRDS(population_simulate, file = file.path(directory, "simulated_pop.Rds"))
saveRDS(real_beta, file = file.path(directory, "simulated_real_beta.Rds"))