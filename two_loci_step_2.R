source("general_functions/common_functions.R")
source("general_functions/generalized_loci_number_functions.R")

# The goal is to implement an underdominance model with assortative mating determined by two separate loci. 
# The A locus determines which pool the individual mates in and the I locus determines the species identity. 
# Fitnesses are affected only by I and are underdominant. 
# Questions: When does linkage develop? When does linkage degrade? Is it stable to have two species? 

# All vectors and matrices are in the order: aa, ab, ac, ad, bb, bc, bd, cc, cd, dd
# Where a = A1.Im, b = A1.Ib, c = A2.Im, d = A2.Ib

# Gives messages and stops when inputs do not make biological sense
# The messages describe what is being checked
# Args: 
#   f: Matrix of initial genotype frequencies
#   d.check: Parameter d, govern non-randomness of mating 
#   h: Parameter than governs the fitness cost of being heterozygous at the I locus
#   r.check: Recombination rate
check.initial <- function(f, d.check, h, r.check) {
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
  
  if (h <= 0 | h > 1) {
     message("h must be greater than zero and less than or equal to one. ")
     failure <- TRUE
  }
  
  if (r.check < 0 | r.check > .5) {
    message("Recombination rate must be between 0 and .5")
    failure <- TRUE
  }
  
  stopifnot(failure == FALSE)
}

# Generates the vector of fitnesses
# Args:
#   geno.mat: a matrix of all the phased genotypes
#   loci.position: the position on the chromosome of the identity locus
#   h: the decrease in fitness of heterozygotes at the identity locus
# Returns:
#   A vector of fitnesses for each phased genotype
selection.assign <- function(geno.mat, loci.position, h) {
  # Determines the number of copies of Ib for each phased genotype
  i_count <- locus_count(geno.mat, loci.position) * -1
  # Determines selection based off the number of copies of B2
  for (position in 1:length(i_count)) {
    if (i_count[position] == 0) {
      i_count[position] <- 1
    }
    if (i_count[position] == -1) {
      i_count[position] <- 1 - h
    }
    if (i_count[position] == -2) {
      i_count[position] <- 1
    }
  }
  i_count
}

# Runs the model for a given set of parameters and initial genotype frequencies
# Args:
#   f.selection: 10 x 1 matrix of initial genotype frequencies
#   h: Parameter than governs the fitness cost of being heterozygous at the I locus
#   d: Parameter d, govern non-randomness of mating 
#   recomb: Recombination rate
# Returns: 
#   A list of: Matrix of genotype frequencies at equilibrium or after the maximum number of generations, whichever occurs first
#              Number of generations run
#              Where or not convergence was reached (0 if reached, 1 if not reached)
run.Model <- function(f.selection, h, d, recomb) {
  # Checks that the input values of number of initial frequencies, d (assortative mating), fitnesses, and recombination make biological sense
  # Gives an message if not and stops
  check.initial(f.selection, d, h, recomb)
  
  # Stores the current generation. The value of generation after exiting the while loop is 1 greater than the number of generations run
  generation <- 1
  # Might stop before this time if convergence is reached
  generation.max <- 10000
  
  # Generates the matrix of phased genotypes
  geno.store <- geno.all(haplo.all(haplo.1(), 2))
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
  haplo.mat <- find.haplo.conv(haplo.all(haplo.1(), 2), geno.store)
  
  # Governs fitnesses, based off of I locus. 
  # Heterozygotes at the I locus are h less fit than homozygotes.
  w <- selection.assign(geno.store, 2, h)
  
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
    message(paste("Did not reach equilibrium. d is ", d, ". h is ", h ))
    convergence <- 1
  } else {
    convergence <- 0
  }
  
  # Returns a list with frequencies after selection, number of generations run, and whether convergence was reached
  list(f.selection = f.selection.store[generation - 1, ], generation = generation - 1, convergence = convergence)

}