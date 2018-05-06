# This is a series of functions that are needed to make common_functions work with a generalized number of loci. 
# The general idea is that all possible haplotypes for a number of loci are stored in a matrix. 
# This is then used to make a matrix of all possible phased genotypes where every set of two consecutive rows 
# is a single phased genotype. These are then used to find vectors and matricies needed in other functions. 
# Because this is generalized, changing the number of loci and the order of the loci on the chromosome requires 
# very few modifications. Alleles are represented as 0 or 1. It is assumed all loci have only two alleles.
# Note that none of these functions require frequencies of genotypes. 


# Returns the matrix of haplotypes with only 1 locus
# Is used in the haplo.all function to pass as the first matrix used in recursion
haplo.1 <- function () {
  matrix(c(0, 1), nrow = 2, ncol = 1)
}

# Enumerates all of the haplotypes with any number of loci assuming 2 alleles per locus
# Alleles are represented as 0 or 1 
# Each row is a distinct haplotype and each column is an allele
# This works by recursion
# Args:
#   start.mat: A matrix of haplotypes for a number of loci (less than or equal to number of loci desired in the end)
#   loci: the number of loci that need to be added to start.mat
# Returns:
#   a matrix of haplotypes with the specified number of loci if no more loci need to be added
haplo.all <- function(start.mat, loci) {
  # If more loci need to be added, adds one more loci and calls this function again with the new matrix
  # as the starting matrix and decreasing the number of loci left to add by one. 
  if (loci - 1 > 0) {
    row <- 1
    new.mat <- matrix(NA, ncol = dim(start.mat)[2] + 1, nrow = dim(start.mat)[1] * 2)
    for (col in 0:1) {
      for (row in 1:dim(start.mat)[1]) {
        new.mat[dim(start.mat)[1] * col + row, ] <- c(col, start.mat[row, ])
      }
    }
    haplo.all(new.mat, loci - 1 )
  }
  # If no more loci need to be added to the matrix, returns the matrix of haplotypes
  else{
    start.mat
  }
}

# This function generates all possible phased genotypes, which are stored in a matrix
# Every 2 consecutive, nonoverlapping lines is a phased genotype. 
# Args:
#   haplo.mat: matrix will all haplotypes for a given number of loci
# Returns:
#   matrix with all phased genotypes 
geno.all <- function(haplo.mat) {
  # Used to store the phased genotype. Rows (2n - 1) and 2n are nth phased genotype 
  geno.mat <- matrix(NA, ncol = dim(haplo.mat)[2], nrow = sum(1:dim(haplo.mat)[1]) * 2)
  # Keeps track of which row in geno.mat is next
  row <- 1
  
  # Puts all combinations of 2 haplotypes into a matrix of phased genotypes
  # Repeats with the order of genotypes flipped are not included (eg Ab x aB and aBx Ab is only included once)
  
  # Loops through haplotypes to determines the 1st haplotype in a phased genotype
   for (haplo_1 in 1: dim(haplo.mat)[1]) {
     # loops through haplotypes to determines the 2nd haplotype in a phased genotype
     for (haplo_2 in haplo_1:dim(haplo.mat)[1]) {
       geno.mat[row, ] <- haplo.mat[haplo_1, ]
       geno.mat[row + 1, ] <- haplo.mat[haplo_2, ]
       row = row + 2
     }
   }
  geno.mat
}

# This function find the recombination array. This determines which 
# genotypes change after recombination between any two loci. 
# Args: 
#   geno.mat: a matrix of all the phased genotypes
# Returns: 
#   an array that is used to determine frequencies after recombination
find.rec.mat <- function(geno.mat) {
  # A vector used to order the haplotypes in a phased genotype so that AB and BA are always recorded in the same order
  # Is equal to (10^(n-1), 10^(n-2), ... , 10) where n is the number of loci
  ten.v <- rep(0, dim(geno.mat)[2])
  for (ten.place in 1:length(ten.v)) {
    ten.v[ten.place] <- 10 ^ (length(ten.v) - ten.place)
  }
  
  # An array that stores the which genotypes change after recombination. It is initially set so that no genotypes change
  # The for loops modify this. Its dimensions are [# of phased genotypes, # phased genotypes, # loci - 1]
  # rec.arry[, , n] determines how phased genotypes change with recombination between the n and n + 1 loci
  rec.matrix <- diag(rep(1,(dim(geno.mat)[1] / 2)), (dim(geno.mat)[1] / 2))
  rec.array <- array(rep(rec.matrix, dim(geno.mat)[2] - 1), dim = c(dim(geno.mat)[1] / 2, dim(geno.mat)[1] / 2, dim(geno.mat)[2] - 1))
  
  # The outer loop is for each of the possible locations of recombination. It is not necessary for only two loci.
  for (position in 1: (dim(geno.mat)[2] - 1)) {
    # Loop changes which phased genotype is being recombined 
    for (row in 1:(dim(geno.mat)[1] / 2)) {
      # Determines what the new phased genotypes are after recombination
      recombine_1 <- geno.mat[row * 2 - 1, ]
      recombine_2 <- geno.mat[row * 2, ]
      recombine <-  matrix(c(recombine_1[1:position], recombine_2[(position+1): dim(geno.mat)[2]], 
                             recombine_2[1:position], recombine_1[(position+1): dim(geno.mat)[2]]),
                           nrow = 2, ncol = dim(geno.mat)[2], byrow = TRUE)
      
      # Flips the order of the haplotypes if they are not in the right order
      # In each haplotype pair, the right order is determined by the number that is obtained by reading 
      # the values of 0 and 1 across a row in the matrix. The row that has a lower number is first in the pair of haplotypes. 
      # (eg [0, 0, 1] becomes 001 is before [0, 1, 0] 010). 
      if (recombine[1, ] %*% ten.v > recombine[2, ] %*% ten.v) {
        recombine <-rbind(recombine[2, ], recombine[1, ])
      }
       
      # If phased genotype is different after recombination, this finds which phased genotype the new phased genotype 
      # corresponds to in the matrix of phased genotypes. Its index (eg the 8th phased genotype in geno.mat) is then used to 
      # change the correct value in the rec.array.
      if (any(recombine != geno.mat[(row * 2 - 1):(row * 2), ])) {
        rec.array[row, row, position] <- 0
        
        # Loops through the phased genotype matrix
        for (find.same in 1:(dim(geno.mat)[1] / 2)) {
          # If the phased genotype matches the recombined genotype, changes the appropriate entry in rec.array
          if (all(recombine == geno.mat[(find.same * 2 - 1):(find.same * 2), ])) {
            rec.array[find.same, row, position] <- 1
          }
        }
      }
    }
  }
  rec.array
}
# Find the matrix used to convert from phased genotype frequencies to haplotypes frequencies. 
# Needed for the mating.pool function
# Args:
#   geno.mat: matrix of all the phased genotypes
#   haplo.mat: matrix of all haplotypes
# Returns:
#   matrix used to convert to haplotype frequencies given phased genotype frequencies. 
find.haplo.conv <- function (haplo.mat, geno.mat) {
  # Matrix to convert to haplotype frequencies, initially all filled with 0s. 
  # Dimensions are [# haplotypes, # phased genotypes]
  haplo.conv <- matrix(rep(0, (dim(haplo.mat)[1]) * dim(geno.mat)[1] / 2), nrow = dim(haplo.mat)[1], ncol = dim(geno.mat)[1] / 2, byrow = TRUE)
  
  # Loops through all of the phased genotypes
  for (row_g in 1:(dim(geno.mat)[1] / 2)) {
    # Loops through all of the haplotypes
    for (row_h in 1:dim(haplo.mat)[1]) {
      # Adds a 1/2 to the [haplotype, genotype] position in haplo.conv for every copy of the a given haplotype in a given genotype 
      if (all(haplo.mat[row_h,] == geno.mat[row_g * 2 - 1,])) {
        haplo.conv[row_h, row_g] <- 1/2
      }
      if (all(haplo.mat[row_h,] == geno.mat[row_g*2,])) {
        haplo.conv[row_h, row_g] <- haplo.conv[row_h, row_g] + 1/2
      }
    }
  }
  haplo.conv
}

# This functions determines the the number of copies of a specific allele in each phased genotype
# Args:
#   geno.mat: a matrix of all the phased genotypes
#   loci.position: the locus at which alleles are being counted
# Returns:
#   A vector with the number of copies of the allele represented by "1" for each phased genotype
locus_count <- function(geno.mat, loci.position) {
  # used to store the number of copies of the "1" allele
  locus.v <- rep(-1, (dim(geno.mat)[1] / 2))
  # Counts the number of copies of the "1" allele for each phased genotype
  for (row in 1:(dim(geno.mat)[1] / 2)) {
    locus.v[row] <- geno.mat[(2 * row) - 1, loci.position] + geno.mat[2 * row, loci.position] 
  }
  locus.v
}

# This function determines the proportion of inidividuals that mate in pool 1
# Args:
#   geno.mat: matrix of all the phased genotypes
#   loci.position: the position of loci that determines assortative mating
#   d: Assortative mating parameter
pool.assign.1 <- function(geno.mat, loci.position, d) {
  # Counts the number of copies of the assorative mating allele
  a_count <- locus_count(geno.mat, loci.position)
  a_return <- rep(NA, length(a_count))
  # Determines which pool each phased genotype mates on based off of the number of copies of the "A2" allele
  for (position in 1:length(a_count)) {
    if (a_count[position] == 0) {
      a_return[position] <- (1 + d) / 2
    }
    if (a_count[position] == 1) {
      a_return[position] <- 1/2
    }
    if (a_count[position] == 2) {
      a_return[position] <- (1 - d) / 2
    }
  }
  a_return
}

# Calculates the allele frequencies for a matrix of genotype frequencies
# Args: 
#   geno.mat: 
#   f.mat: Matrix of genotype frequencies. Each row is a the genotypes at a given time for a given model run. Each column is the phased genotype
# Return:
#   Matrix with the allele frequencies of each allele for each row in f.mat

allele.freq <- function(geno.mat, f.mat) {
  # Gives the number of loci
  num.locus <- dim(geno.mat)[2]
  
  # Matrix used to store the allele frequencies 
  f <- matrix(NA, ncol = num.locus, nrow = dim(f.mat)[1])
  
  # Calculates the allele frequencies for each allele for each generation before and after mating
  for(i in 1:num.locus){
    f.locus <- locus_count(geno.mat, i)
    f[ , i] <- f.mat %*% f.locus/2
  }
  f
}