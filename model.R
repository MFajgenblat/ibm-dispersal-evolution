#-------------------------------------------------------------------------------
# Landscape, metacommunity & genetics settings
#-------------------------------------------------------------------------------

# Number of patches, environmental niche values and species (identical by design)
N_patches <- N_environments <- N_max_species <- 5

# Environmental values of the patches
environment <- seq(0.2, 0.8, length.out = N_environments)

# Available species pool configurations
species_pools <- list(c(3), c(1,5), c(1,2,3,4,5))

# Species that cannot evolve dispersal
species_fixed_dispphen <- list(c(), c(1,2,3,4,5), c(3), c(1,2,4,5), c(2,3,4,5))

# Size of the haploid set of loci coding the environmental phenotype
# Identitcal to the size of the haploid set of loci coding the dispersal phenotype
N_loci <- 50

# Number of unique mainland haplotypes per species
N_mainland_haplotypes <- 100

# Number of unique mainland genotypes per species
N_mainland_genotypes <- 10000

# Duration of a single simulation
N_generations <- 100

#-------------------------------------------------------------------------------
# Fixed model parameters
#-------------------------------------------------------------------------------

# Local carrying capacity of patch
K <- 250

# Per capita maximum growth rate
r_max <- 10

# Environmental selection strength
sigma <- 0.002

# Genetic diversity in mainland population (max = 0.5)
gamma <- 0.1

# Mutation rate
mu = 10e-4

# Initial dispersal propensity
delta_init <- 0.2

#-------------------------------------------------------------------------------
# Variable model parameters: default settings
#-------------------------------------------------------------------------------

# Selected species pool
species_pool <- 3

# Species that cannot evolve dispersal scenario
evodisp_scenario <- 1

# Fixed dispersal phenotype for species that cannot evolve dispersal
delta_fix <- 0.25/5

# Expected amount of mainland individuals seeded per patch per generation
lambda <- 1

# Cost of dispersal (probability of additional mortality for dispersing phenotypes)
chi <- 1

# Probability of complete patch extinction probability per generation
psi <- 0

#-------------------------------------------------------------------------------
# Utility functions
#-------------------------------------------------------------------------------

# Robust resample function without quirks if only one item is to be sampled
resample <- function(x, ...) x[sample.int(length(x), ...)]

# Round to either lower or upper integer value, with the decimal part of the number as probability
probabilistic_round <- function(x) {
  return(floor(x) + resample(c(floor(x - floor(x)), ceiling(x - floor(x))), 1, replace = F, prob = c(x - floor(x), 1 - (x - floor(x)))))
}

# Create standing genetic variation through Robertsonian chromosomal translocations
# Alleles are reshuffeled across loci, i.e. keeping the phenotype identical
reshuffle_alleles <- function(genotype, rate) {
  number <- min(rpois(1, rate), probabilistic_round(N_loci/2))
  locations <- resample(1:N_loci, number, replace = F)
  antilocations <- resample(which(!(1:N_loci %in% locations)), number)
  locations_values_stored <- genotype[locations]
  antilocations_values_stored <- genotype[antilocations]
  genotype[locations] <- antilocations_values_stored
  genotype[antilocations] <- locations_values_stored
  return(genotype)
}

#-------------------------------------------------------------------------------
# Column legend of the output dataframe
#-------------------------------------------------------------------------------

#   1 --> Generation
#   2 --> Patch
#   3 --> Species
#   4:(3+N_loci*2) --> Diploid set of environmental loci
#   (4+N_loci*2):(3+N_loci*4) --> Diploid set of dispersal loci

#-------------------------------------------------------------------------------
# Defining the IBM function
#-------------------------------------------------------------------------------

run_ibm <- function(species_pool = species_pool,
                    evodisp_scenario = evodisp_scenario,
                    delta_fix = delta_fix,
                    lambda = lambda,
                    chi = chi,
                    psi = psi,
                    replicate_id = 1) {
  
  # MAINLAND SETUP
  
  # Create modal environmental genotypes for each species
  modal_env_genotypes <- matrix(0, nrow = N_max_species, ncol = N_loci)
  for (i in 1:N_max_species) {modal_env_genotypes[i,resample(1:N_loci, probabilistic_round(environment[i]*N_loci))] <- 1}
  
  # Create modal dispersal genotypes for each species
  modal_disp_genotypes <- matrix(0, nrow = N_max_species, ncol = N_loci)
  for (i in 1:N_max_species) {modal_disp_genotypes[i,resample(1:N_loci, probabilistic_round(delta_init*N_loci))] <- 1}
  
  # Create haplotype pools for each species
  mainland_env_haplotypes <- lapply(1:N_max_species, function(i) t(replicate(N_mainland_haplotypes, reshuffle_alleles(modal_env_genotypes[i,], probabilistic_round(gamma*N_loci)))))
  mainland_disp_haplotypes <- lapply(1:N_max_species, function(i) t(replicate(N_mainland_haplotypes, reshuffle_alleles(modal_disp_genotypes[i,], probabilistic_round(gamma*N_loci)))))
  
  # Create the mainland combination by randomly combining two haplotypes per species
  mainland_individuals <- cbind(Species = rep(1:N_max_species, each = N_mainland_genotypes),
                                do.call(rbind, lapply(rep(1:N_max_species, each = N_mainland_genotypes), function(i) c(t(mainland_env_haplotypes[[i]][resample(1:N_mainland_haplotypes, 2, replace = T),])))),
                                do.call(rbind, lapply(rep(1:N_max_species, each = N_mainland_genotypes), function(i) c(t(mainland_disp_haplotypes[[i]][resample(1:N_mainland_haplotypes, 2, replace = T),])))))
  
  # POPULATION DYNAMICS
  
  for (generation in 1:N_generations) {
    
    # Subsetting the current generation, augmented by any mainland individuals
    if (generation == 1) {
      history <- cbind(Generation = NA, Patch = NA, Species = NA, Environmental_phenotype = NA, dispersal_phenotype = NA)[-1,]
      current_generation <- cbind(Generation = NA, Patch = NA, Species = NA, matrix(NA, ncol = 4*N_loci, nrow = 1))[-1,]
    } else {
      current_generation <- upcoming_generation
    }
    
    # Mainland rain
    current_generation <- rbind(current_generation,
                                do.call(rbind, lapply(1:N_patches, function(i) cbind(Generation = generation, Patch = i, mainland_individuals[resample(which(mainland_individuals[,1] %in% species_pools[[species_pool]]), rpois(1, lambda), replace = T),,drop = F]))))
    
    
    # If the patches have at least one individual:
    if (nrow(current_generation) > 0) {
      
      # Fitness assessment
      fitness <- exp(-(environment[current_generation[,2]] - rowMeans(current_generation[,4:(4+N_loci*2),drop = F]))^2/sigma)
      
      # Logistic coefficient
      popsizes <- sapply(1:N_patches, function(i) sum(current_generation[,2] == i))
      logistic_coefficient <- r_max * (1 - popsizes / K) + popsizes / K
      logistic_coefficient[popsizes > K] <- (K / popsizes)[popsizes > K]
      
      # Computing reproductive output
      reproductive_output <- rpois(nrow(current_generation), fitness * logistic_coefficient[current_generation[,2]])
      
      # If the reproductive output is positive:
      if (sum(reproductive_output) > 0) {
        
        # Mate selection
        partner_1 <- c(unlist(sapply(1:nrow(current_generation), function(i) rep(i, reproductive_output[i]))))
        partner_2 <- c(unlist(sapply(1:nrow(current_generation), function(i) resample(which((current_generation[,2] == current_generation[i,2]) & (reproductive_output > 0)), reproductive_output[i], replace = T))))
        
        # Recombination
        recombination_matrix <- replicate(2, matrix(rbinom(sum(reproductive_output)*N_loci, 1, prob = 0.5), nrow = sum(reproductive_output), ncol = N_loci), simplify = F)
        recombined_genotypes <- cbind(recombination_matrix[[1]] * current_generation[partner_1,4:(3+N_loci)] + (1 - recombination_matrix[[1]]) * current_generation[partner_1,(4+N_loci):(3+2*N_loci)],
                                      recombination_matrix[[1]] * current_generation[partner_2,4:(3+N_loci)] + (1 - recombination_matrix[[1]]) * current_generation[partner_2,(4+N_loci):(3+2*N_loci)],
                                      recombination_matrix[[2]] * current_generation[partner_1,(4+2*N_loci):(3+3*N_loci)] + (1 - recombination_matrix[[2]]) * current_generation[partner_1,(4+3*N_loci):(3+4*N_loci)],
                                      recombination_matrix[[2]] * current_generation[partner_2,(4+2*N_loci):(3+3*N_loci)] + (1 - recombination_matrix[[2]]) * current_generation[partner_2,(4+3*N_loci):(3+4*N_loci)])
        
        # Mutation
        where2mutate <- matrix(rbinom(sum(reproductive_output)*N_loci*4, 1, mu), nrow = sum(reproductive_output), ncol = N_loci*4)
        mutated_genotypes <- (1 - where2mutate) * recombined_genotypes + where2mutate * (1 - recombined_genotypes)
        
        # Dispersal
        dispersal_phenotypes <- rowMeans(mutated_genotypes[,(1+N_loci*2):(N_loci*4),drop = F])/5
        dispersal_phenotypes[which(unlist(sapply(1:nrow(current_generation), function(i) rep(current_generation[i,3], reproductive_output[i]), simplify = F)) %in% species_fixed_dispphen[[evodisp_scenario]])] <- rep(delta_fix/5, sum(reproductive_output))
        
        # Compiling the upcomping generation
        natal_patches <- unlist(sapply(1:nrow(current_generation), function(i) rep(current_generation[i,2], reproductive_output[i])))
        upcoming_generation <- cbind(Generation = generation,
                                     Patch = sapply(1:sum(reproductive_output), function(i) resample(c(natal_patches[i], 1:N_patches), 1, replace = F, prob = c(1 - dispersal_phenotypes[i], dispersal_phenotypes[i]*c(rep(1/N_patches, N_patches))))),
                                     Species = unlist(sapply(1:nrow(current_generation), function(i) rep(current_generation[i,3], reproductive_output[i]), simplify = F)),
                                     mutated_genotypes)[which(rbinom(sum(reproductive_output), 1, chi*dispersal_phenotypes) == 0),,drop = F]
        
        # Patch extinction
        upcoming_generation <- upcoming_generation[upcoming_generation[,2] %in% which(rbinom(5, 1, rep(psi, 5)) == 0),,drop = F]
        
        # Add new generation to metacommunity history
        history <- rbind(history,
                         data.frame(upcoming_generation[,1:3,drop = F],
                                    Environmental_phenotype = rowMeans(upcoming_generation[,4:(3+N_loci*2), drop = F]),
                                    Dispersal_phenotype = rowMeans(upcoming_generation[,(4+N_loci*2):(3+N_loci*4), drop = F])))
        
      
      # If there is no reproductive output:
      } else {
        
        upcoming_generation <- cbind(Generation = NA, Patch = NA, Species = NA, matrix(NA, ncol = 4*N_loci, nrow = 1))[-1,]
        
      }
    
    # If all patches are empty:
    } else {
      
      upcoming_generation <- cbind(Generation = NA, Patch = NA, Species = NA, matrix(NA, ncol = 4*N_loci, nrow = 1))[-1,]
      
    }
    
  }
  
  # Compile, randomly sample and return history
  return(data.frame(species_pool = species_pool,
                    evodisp_scenario = evodisp_scenario,
                    delta_fix = delta_fix,
                    lambda = lambda,
                    chi = chi,
                    psi = psi,
                    replicate_id = replicate_id,
                    history)[sample(1:nrow(history), round(nrow(history)/100)),])
  
}
