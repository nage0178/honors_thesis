source("general_functions/common_functions.R")

# The goal is to program an underdominance model with assortative mating with a single locus. 
# This locus determines both species identity and which mating pool the individual joins.
# Questions:  Is having two species ever stable? If so, what are the conditions? 

# All genotypes are stored in the order IM.IM, IM.IB, IB.IB

# Gives messages and stops when inputs do not make biological sense
# The messages describe what is being checked
# Args: 
#   geno: Number of phased genotypes
#   haplo: Number of haplotypes
#   f: Matrix of initial genotype frequencies
#   d.check: Parameter d, govern non-randomness of mating 
#   fitness: Vector of fitnesses for each genotype
check.initial <- function(geno, haplo, f, d.check, fitness) {
  failure <- FALSE
  
  if (floor(geno) != geno | geno < 1) {
    message("The number of phased genotypes must be a positive integer")
    failure <- TRUE
  }
  
  if (floor(haplo) != haplo | haplo < 1) {
    message("The number of haplotypes must be a positive integer")
    failure <- TRUE
  }
  
  if (geno != sum(c(1:haplo))) {
    message("The number of phased genotypes does not match the number of haplotypes.")
    failure <- TRUE
  }
  
  if (abs(1 - sum(f)) > 10^-12) {
    message("Initial frequencies do not sum to one.")
    failure <- TRUE
  }
  
  if (any(f < 0)) {
    message("Initial frequencies must be positive.")
    failure <- TRUE
  }
  
  if (d.check < 0 | d.check > 1) {
    message("d must be between zero and one. ")
    failure <- TRUE
  }
  
  if (any(fitness < 0)) {
    message("Fitnesses must be nonnegative.")
    failure <- TRUE
  }
  
  if (geno != 3) {
    message("Number of phased genotypes is not equal to 3. The genotype to haplotype converstion matrix (haplo.mat) needs to be changed in the run.Model function.")
    failure <- TRUE
  }
 
  if (length(fitness) != geno) {
    message("The number of fitnesses must equal the number of genotypes. ")
    failure <- TRUE
  }
  
  stopifnot(failure == FALSE)
}

# Runs the model for set of parameters and initial genotype frequencies
# Args: 
#   num.geno: Number of phased genotypes
#   num.haplo: Number of haplotypes
#   f.selection: Matrix of initial genotype frequencies
#   w: Vector of fitnesses for each genotypes
#   d: Parameter d, govern non-randomness of mating 
#   show.plots: Determines if graphs of genotype frequencies over time should be graphed (TRUE/FALSE). Default is TRUE
#   multiple.plots: Determines if the final frequencies (FALSE) or the frequencies for each generation (TRUE) should be retuned. Default is FALSE
# Returns: 
#   If multiple.plots is false, matrix of genotype frequencies at equilibrium or after the maximum number of generations, whichever occurs first
#   if multiple.plots is true, list with matrix of genotype frequncies after selection for each generation and the number of generations run
run.Model <- function(num.geno, num.haplo, f.selection, w, d, show.plots=TRUE, multiple.plots = FALSE) {
  # Checks that the input values of number of phased genotypes, number of haplotypes, 
  # initial frequencies, d (assortative mating) and fitnesses make biological sense
  # Gives an message if not and stops
  check.initial(num.geno, num.haplo, f.selection, d, w)
  
  # Stores the current generation. The value of generation after exiting the while loop is 1 greater than the number of generations run
  generation <- 1
  # Might stop before this time if convergence is reached
  generation.max <- 10000
  
  # Stores the frequencies immediately after mating for each generation. This is the first event that happens each generation. 
  f.mating.store <- matrix(NA, ncol = num.geno, nrow = generation.max) 
  # Stores the frequencies immediately after selection for each generation. This happens after mating in each generation. 
  f.selection.store <- matrix(NA, ncol = num.geno, nrow = generation.max)
  
  # Used in the mating.pool function. 
  # Determines the proportion of each genotype in each mating pool according to d
  pool.1 <- c((1 + d) / 2, 1/2, (1 - d) / 2)
  pool.2 <- rep(1, num.geno) - pool.1
  # Used to convert genotype frequencies to haplotype frequencies
  haplo.mat<- matrix(c(1, 1/2, 0,
                       0, 1/2, 1), 
                     ncol = num.geno, nrow = num.haplo, byrow = TRUE)

  # Finds the equilibrium frequencies of each genotype
  # First the individuals mate in two pools, which is determined by d. Then selection occurs according to w
  # This is repeated for generation.max generations or until the difference between each genotype frequency is less than 10^-12, whichever comes first
  delta.f.mating <- rep(1, num.geno) # altered in the loop as the stopping condition
  while (generation < generation.max + 1 && any(delta.f.mating > 10^-12)) {
    # Calculates the genotype frequencies after mating
    f.mating <- mating.pool(f.selection, pool.1, haplo.mat, num.geno) + 
                mating.pool(f.selection, pool.2, haplo.mat, num.geno)
    f.mating.store[generation, ] <- f.mating
    
    # Calculates the genotype frequencies after selection
    f.selection <- selection(f.mating, w)
    f.selection.store[generation, ] <- t(f.selection)
    
    # Calculates the difference in frequencies between generations
    if (generation > 1) {
      delta.f.mating <- abs(f.mating.store[generation, ] - f.mating.store[generation-1, ])
    }
    generation <- generation + 1
  }
  
  if (generation > generation.max) {
    message("Did not reach equilibrium.")
  }
  
  if (show.plots) {
    par(mfrow = c(1,2))
    
      # Plots the frequencies after mating
      matplot(c(1:(generation-1)), f.mating.store[1:generation - 1, ], type = "l", lwd = 3,
              xlab = "Generation", ylab = "Genotype frequency",
              ylim = c(0, max(f.selection.store, f.mating.store, na.rm = T)), main = "After Mating")
      legend(legend = c("MM", "MB", "BB"), "right", col = 1:3, lty = 1:3, lwd = 3)

      # Plots the frequencies after selection
      matplot(c(1:(generation-1)), f.selection.store[1:generation - 1, ], type = "l", lwd = 3,
              xlab = "Generation", ylab = "Genotype frequency",
              ylim = c(0, max(f.selection.store, f.mating.store, na.rm = T)), main = "After Selection")
      legend(legend = c("MM", "MB", "BB"), "right", col=1:3, lty = 1:3, lwd = 3)
  }
  
  if(multiple.plots){
    # Returns a list of the frequencies after selection and the number of generations run
    list(f.selection = f.selection.store, generation = generation - 1)
  } else{
    # Returns final frequencies
    f.selection
  }

}

