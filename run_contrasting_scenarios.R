#-------------------------------------------------------------------------------
# Setting up furrr for parallel running of the IBMs
#-------------------------------------------------------------------------------

library(furrr)
plan(multisession,workers = availableCores())

#-------------------------------------------------------------------------------
# Loading the IBM
#-------------------------------------------------------------------------------

source("model.R")

#-------------------------------------------------------------------------------
# Generating a parameter grid
#-------------------------------------------------------------------------------

parameter_grid <- expand.grid(species_pool = 3,
                              evodisp_scenario = c(1,2,4,5),
                              delta_fix = c(0.25/5, 0.25*5),
                              lambda = c(0.1,1,10),
                              chi = 0,
                              psi = c(0,0.01),
                              replicate_id = 1:20)

#-------------------------------------------------------------------------------
# Running the the IBMs
#-------------------------------------------------------------------------------

all_histories <- do.call(rbind,
                         (future_map(1:nrow(parameter_grid),
                                     function(run) run_ibm(species_pool = parameter_grid$species_pool[run],
                                                           evodisp_scenario = parameter_grid$evodisp_scenario[run],
                                                           delta_fix = parameter_grid$delta_fix[run],
                                                           lambda = parameter_grid$lambda[run],
                                                           chi = parameter_grid$chi[run],
                                                           psi = parameter_grid$psi[run],
                                                           replicate_id = parameter_grid$replicate_id[run]),
                                     .progress = T)))

#-------------------------------------------------------------------------------
# Writing the results
#-------------------------------------------------------------------------------

write.table(all_histories, "results_contrasting_scenarios.csv", col.names = T, row.names = F, quote = F, sep = ";")
