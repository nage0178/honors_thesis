source("general_functions/common_functions.R")
source("general_functions/generalized_loci_number_functions.R")


# The goal is to implement an underdominance model with assortative mating determined by two separate loci with a third locus that represents metamorphose time.
# This locus, E,  has directional selection. 
# The A locus determines which pool the individual mates in and the I locus determines the species identity. 
# Fitnesses are affected by I, which is underdominant, and E. 
# Fitnesses are multiplicative. 
# Loci are in the order AIE.
# Questions: When does linkage develop? When does linkage degrade? Is it stable to have two species? Does metamorphose rate fix? 

# Gives messages and stops when inputs do not make biological sense
# The messages describe what is being checked
# Args: 
#   f: Matrix of initial genotype frequencies
#   d.check: Parameter d, govern non-randomness of mating 
#   h.I: the decrease in fitness of heterozygotes at the identity locus
#   h.E: Together with s.E, governs the fitness of being a E_b E_m in a shallow pond
#   s.E: Governs the fitness cost of being E_b E_b in shallow ponds
#   r.check: vector of recombination rates
check.initial <- function(f, d.check, h.I, r.check, h.E, s.E) {
  failure <- FALSE

  if (abs(1 - sum(f)) > 10^-12) {
    message("Initial frequencies do not sum to one.")
    failure <- TRUE
  }
  
  if (any(f < 0)) {
    message("Initial frequencies must be positive.")
    failure <- TRUE
  }
  
  if (d.check < 0 | d.check > 1) {
    message("d is not between zero and one. ")
    failure <- TRUE
  }
  
  if (h.I <= 0 | h.I > 1) {
    message("h.I must be greater than zero and less than or equal to one. ")
    failure <- TRUE
  }
  
  if (h.E <= 0 | h.E > 1) {
    message("h.E must be greater than zero and less than or equal to one. ")
    failure <- TRUE
  }
  
  if (s.E <= 0 | s.E > 1) {
    message("s.E must be greater than zero and less than or equal to one. ")
    failure <- TRUE
  }
  
  
  if (any(r.check < 0) | any(r.check > .5)) {
    message("Recombination rate must be between 0 and .5")
    failure <- TRUE
  }
  
  stopifnot(failure == FALSE)
}

# Generates the vector of fitnesses
# Args:
#   geno.mat: a matrix of all the phased genotypes
#   loci.position.I: the position on the chromosome of the identity locus
#   loci.position.E: the position on the chromosome of the metamorphose rate locus
#   h.I: the decrease in fitness of heterozygotes at the identity locus
#   h.E: Together with s.E, governs the fitness of being a E_b E_m in a shallow pond
#   s.E: Governs the fitness cost of being E_b E_b in shallow ponds
# Returns:
#   A vector of fitnesses for each phased genotype
selection.assign <- function(geno.mat, loci.position.I, loci.position.E, h.I, h.E, s.E) {
  # Determines the number of copies of Ib for each phased genotype
  i_count <- locus_count(geno.mat, loci.position.I) * -1
  # Determines selection based off the number of copies of B2
  for (position in 1:length(i_count)) {
    if (i_count[position] == 0) {
      i_count[position] <- 1
    }
    if (i_count[position] == -1) {
      i_count[position] <- 1 - h.I
    }
    if (i_count[position] == -2) {
      i_count[position] <- 1
    }
  }
  
  # Determines the number of copies of E.b for each phased genotype
  E_count <- locus_count(geno.mat, loci.position.E) * -1
  # Determines selection based off the number of copies of E.b
  for (position in 1:length(i_count)) {
    if (E_count[position] == 0) {
      E_count[position] <- 1
    }
    if (E_count[position] == -1) {
      E_count[position] <- 1 - h.E * s.E
    }
    if (E_count[position] == -2) {
      E_count[position] <- 1 - s.E
    }
  }
  # Multiplies fitness based on I and E loci and returns the total fitness
  fitness <- i_count * E_count
  fitness
}

# Runs the model for a given set of parameters and initial genotype frequencies
# Args:
#   f.selection: 36 x 1 matrix of initial genotype frequencies
#   h.I: the decrease in fitness of heterozygotes at the identity locus
#   h.E: Together with s.E, governs the fitness of being a E_b E_m in a shallow pond
#   s.E: Governs the fitness cost of being E_b E_b in shallow ponds
#   d: Parameter d, govern non-randomness of mating 
#   recomb: Recombination rates in a vector
# Returns: 
#   A list of: Matrix of genotype frequencies at equilibrium or after the maximum number of generations, whichever occurs first
#              Number of generations run
#              Where or not convergence was reached (0 if reached, 1 if not reached)
run.Model <- function(f.selection, h.I, h.E, s.E, d, recomb) {
  # Checks that the input values of number of initial frequencies, d (assortative mating), fitnesses, and recombination make biological sense
  # Gives an message if not and stops

  check.initial(f.selection, d, h.I, recomb, h.E, s.E)
  
  # Stores the current generation. The value of generation after exiting the while loop is 1 greater than the number of generations run
  generation <- 1
  # Might stop before this time if convergence is reached
  generation.max <- 10000
  
  # Generates the matrix of phased genotypes
  geno.store <- geno.all(haplo.all(haplo.1(), 3))
  # Generates the array that determines which phased genotyps change during recombination
  rec.array <- find.rec.mat(geno.store)
  
  num.geno <- dim(geno.store)[1] / 2 # Number of genotypes
  num.haplo <- dim(haplo.all(haplo.1(), 2))[1] # Number of haplotypes
  
  # Stores the frequencies immediately after selection for each generation. This happens after mating each generation.  
  f.selection.store <- matrix(NA, ncol = num.geno, nrow = generation.max)
  # Stores the frequencies immediately after mating for each generation. This happens after first each generation. 
  f.mating.store <- matrix(NA, ncol = num.geno, nrow = generation.max) 
  
  # Used in the mating.pool function. 
  # Determines the proportion of each genotype in each mating pool according to d
  pool.1 <- pool.assign.1(geno.store, 1, d)
  pool.2 <- rep(1, num.geno) - pool.1
  # Used to convert genotype frequencies to haplotype frequencies
  haplo.mat <- find.haplo.conv(haplo.all(haplo.1(), 3), geno.store)
  
  # Governs fitnesses, based off of I locus. 
  # Heterozygotes at the I locus are h less fit than homozygotes.
  # geno.mat, loci.position.I, loci.position.E, h.I, h.E, s.E
  w <- selection.assign(geno.store, 2, 3, h.I, h.E, s.E)
  
  # Finds the equilibrium frequencies of each genotype
  # First the individuals mate in two pools, which is determined by d. Then selection occurs according to w.
  # This is repeated for generation.max generations or until the difference between each genotype frequency is less than 10^-14, whichever comes first
  delta.f.mating <- rep(1, num.geno) # stores the difference in frequencies between generations, altered in the loop as the stopping condition
  while (generation < generation.max + 1 && any(delta.f.mating > 10^-14)) {
    # Calculates the genotype frequencies after mating
    f.mating.1 <- mating.pool(f.selection, pool.1, haplo.mat, num.geno, recomb, rec.array)
    f.mating.2 <- mating.pool(f.selection, pool.2, haplo.mat, num.geno, recomb, rec.array)
    f.mating.store[generation, ] <- f.mating.1 + f.mating.2
    
    
    # Calculates the genotype frequencies after selection (and density regulation)
    f.selection <- selection(f.mating.store[generation, ], w)
    f.selection.store[generation, ] <- t(f.selection)
    
    # Calculates the difference in frequencies between generations
    if (generation > 1) {
      delta.f.mating <- abs(f.mating.store[generation, ] - f.mating.store[generation - 1, ])
    }
    generation <- generation + 1
  }
  
  # Gives a message if convergence isn't reached.
  # Sets convergence to be 1 if convergences isn't reached and 0 if it is reached. 
  if (generation > generation.max) {
    message(paste("Did not reach equilibrium." ))
    convergence <- 1
  } else {
    convergence <- 0
  }

  
  # Returns a list with frequencies after selection, number of generations run, and whether convergence was reached
  list(f.selection = f.selection.store[generation - 1, ], generation = generation - 1, convergence = convergence)
}