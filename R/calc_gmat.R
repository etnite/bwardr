#' Calculate Realized Genomic Relationship Matrices
#'
#' @param data Dataframe or matrix with individuals in rows and markers in columns.
#'   Markers should be encoded in minor-allele dosage format
#'   (i.e. 0 = homozygous for major allele, 1 = heterozygous, 2 = homozygous
#'   for minor allele).
#' @param method String to select the method to use for calculating the 
#' relationship matrix. Available methods are:
#' \enumerate{
#'   \item "GOF" - observed allele frequencies (VanRaden, 2008)
#'   \item "GD" - weighted markers by recipricals of expected variance (Forni et al., 2011)
#'   \item "G05" - allele frequencies fixed at 0.5 (Forni et al., 2011)
#'   \item "GM" - allele frequencies fixed at mean for each locus (Forni et al., 2011)
#' }
#' @return The full realized genetic relationship matrix (G)
#' @details This function is adapted from an example from Chapter 11 of
#'   Genetic Data Analysis for Plant and Animal Breeding, by Isik, Holland and
#'   Maltecca, 2017 (DOI 10.1007/978-3-319-55177-7). However, their code is
#'   more complex, as it includes one additional computation method that uses
#'   observed pedigree data, and it outputs both realized genetic and additive
#'   relationships. This function only returns the realized genetic relationship
#'   matrix.
#' @export
calc_gmat <- function(data, method = 1) {

  ## Sanity check on supplied method integer
  if (!method %in% c("GOF", "GD", "G05", "GM")) {
    stop("Please select 'GOF', 'GD', 'G05', or 'GM' for G matrix calculation method")
  }
  
  ## Create marker matrix encoded (-1, 0, 1)
  ## For homozgous, het, other homozygous
  M <- as.matrix(data) - 1
  
  ## Calculate or estimate p for each SNP, where p = (2 * minor allele freq)
  if (method %in% c("GOF", "GD")) {
      p1 <- round((apply(M, 2, sum) + nrow(M))/(nrow(M) * 2), 3)
      p <- 2 * (p1 - 0.5)
  } else if (method == "G05") {
      p1 <- array(0.5, ncol(M))
      p <- p1
  } else {
      p1 <- round((apply(M, 2, mean)), 3)
      p <- 2 * (p1 - 0.5)
  }
  
  ## Scale M matrix to create Z matrix
  P <- matrix(p, byrow = T, nrow = nrow(M), ncol = ncol(M))
  Z <- as.matrix(M - P)
  
  ## Estimate G matrix from Z matrix
  if (method %in% c("GOF", "G05", "GM")) {
      b <- 1 - p1
      c <- p1 * b
      d <- 2 * (sum(c))
      ZZt <- Z %*% t(Z)
      G <- (ZZt/d)
  } else {
      D <- 1/(ncol(M) * (2 * p1 * (1 - p1)))
      G <- Z %*% (D * t(Z))
  }
  
  return(G)
}