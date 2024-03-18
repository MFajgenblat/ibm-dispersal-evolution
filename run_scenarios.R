################################################################################
# This file launches the replicate runs for all scenarios
################################################################################

#-------------------------------------------------------------------------------
# Setting up furrr for parallel running of the IBMs
#-------------------------------------------------------------------------------

library(furrr)
plan("multisession", workers = 36)

source("model.R")

#-------------------------------------------------------------------------------
# Generating a parameter grid
#-------------------------------------------------------------------------------

parameter_grid <- expand.grid(Species_pool_scenario = 1:16,
                              Lambda = c(0.1,1,10),
                              Chi = c(0,0.1),
                              Psi = c(0,0.01),
                              Mu = c(10e-4,10e-5),
                              replicate_id = 1:25)

#-------------------------------------------------------------------------------
# Running the the IBMs
#-------------------------------------------------------------------------------

all_histories <- do.call(rbind,
                         (future_map(1:nrow(parameter_grid),
                                     function(run) run_ibm(Species_pool_scenario = parameter_grid$Species_pool_scenario[run],
                                                           Lambda = parameter_grid$Lambda[run],
                                                           Chi = parameter_grid$Chi[run],
                                                           Psi = parameter_grid$Psi[run],
                                                           Mu = parameter_grid$Mu[run],
                                                           replicate_id = parameter_grid$replicate_id[run]),
                                     .progress = T)))

#-------------------------------------------------------------------------------
# Writing the results
#-------------------------------------------------------------------------------

write.table(all_histories, paste0("results_scenarios.csv"), col.names = T, row.names = F, quote = F, sep = ";")
