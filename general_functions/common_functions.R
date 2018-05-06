#This is a file of functions that are used in multiple models/files. 

# Calculates the genotype frequencies after random mating within the pool. 
# First the genotype frequencies in the pool are calculated. Then recombination occurs. These frequencies 
# are converted to haploid frequencies. Offspring genotype frequencies are calculated by finding the 
# outer product of the haploid frequencies. Frequencies are then 
# scaled so they add to the same proportion as entered the pool. 
# Args: 
#   f: Matrix of the genotype frequencies of the entire population
#   pool.matrix: A matrix that determines the proportion of each genotype that mates in a given pool
#   hap.mat: A matrix to convert phased genotype frequencies to haplotype frequencies
#   num.geno.m: Number of phased genotypes
#   rec: Recombination rate
#   rec.array: array that determines which phased genotypes change with recombination
# Returns: 
#   A matrix of the gneotype frequencies after mating random mating within the pool, scaled so the sum of the frequencies before and after mating are the same
mating.pool <- function(f, pool.matrix, hap.mat, num.geno.m, rec, rec.array) {
  # Determines the frequencies of each phased genotype that mate in the pool
  f.pool.mating <- pool.matrix * f
  
  # Determines phased genotype frequencies after recombination
  # If there is only one locus (two genotypes), recombination doesn't occur
  if (length(f) > 3) {
    f.pool.mating <- recombination(rec.array, rec, f.pool.mating)
  }
  
  # Calculates the haplotype frequencies and the frequency of each haplotype cross (which
  # is the genotype frequency for the offspring)
  f.hap <- c(hap.mat %*% f.pool.mating)

  f.unscale <- haplo.comb(f.hap, num.geno.m)
  
  # Frequencies are scaled because initial frequencies do not add to one if there are multiple mating pools.
  # Scaling causes the sum of the frequencies after mating to be the same as the sum of frequencies before mating.
  f.scale <- f.unscale / sum(f.unscale) * sum(f.pool.mating)
  f.scale
}

# Randomly combines haplotypes, used in mating pool function and to find inital frequencies
# Args: 
#   f.hap: Vector of haplotype frequencies
#   num.geno.m: Number of phased genotypes
# Returns: 
#   Vector of phased genotype frequencies
haplo.comb <- function(f.hap, num.geno.m) {
  outprod <- f.hap %o% f.hap
  
  f.position <- 1    # Used to determine which position of the genotype frequencies in f.unscale
  f.unscale <- matrix(NA, ncol = 1, nrow = num.geno.m)
  
  # Puts the frequency of each haploid cross into a vector. If the two parental haplotypes are the same, 
  # the frequency of the cross is equal to the frequency of the offspring genotype. 
  # If the parental haplotypes are different, the frequencies of the two crosses are added together. (eg Axa + axA)
  for (row in 1:length(f.hap)) {
    for (column in row:length(f.hap)) {
      if (row == column) {
        f.unscale[f.position] <- outprod[row, column]
      } else {
        f.unscale[f.position] <- outprod[row, column] + outprod[column, row]
      }
      f.position <- f.position + 1
    }
  }
  f.unscale
}

# Calculates the genotype frequencies after selection 
# Args: 
#   f: Matrix of genotype frequencies of the entire population
#   selection.v: A vector that determines selection on each genotype
#   density.reg: Determines if density regulation occurs. Default is FALSE
# Returns: 
#   Matrix of genotype frequencies after selection
selection <- function(f, selection.v, density.reg = FALSE) {
  after.selection <- f * selection.v
  w.average <- sum(after.selection)
  f.selection <- after.selection / w.average
  # Density regulates for models with multiple selection pools that density regulate separately
  if (density.reg) {
    f.selection * sum(f)
  } else {
    f.selection
  }
}

# Calculates the phased genotype frequencies after recombination. 
# Called from the mating.pool function
# Args: 
#   rec.array: An array that determines which genotypes change with recombination
#              rec.array[,,n] is the matrix of how recombination between the n and n+1 loci 
#              change genotype frequencies.
#   f: The matrix of genotype frequencies for the entire population
#   recomb: Recombination rate, with more than two loci this is a vector
# Returns: 
#   The frequency of all genotypes after recombination
recombination <- function(rec.array, recomb, f) {
  for (position in 1: (dim(rec.array)[3])) {
    f <- f * (1 - recomb[position]) + rec.array[, , position] %*% (f * recomb[position])
  }
  f
}