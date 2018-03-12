#' Perform Eigenstrat marker PCA
#' 
#' @param geno A matrix of SNP data with SNPs in rows, genotypes in columns.
#'   SNPs should be encoded in minor-allele dosage format, so that:
#'   2 = homozygous for minor allele
#'   1 = heterozygous
#'   0 = homozygous for major allele
#'   NA = missing data
#' @return The Eigen decomposition of the input genotypic matrix
#' @details This function performs principle component analysis on a matrix of
#'   SNP data in the method of EIGENSTRAT (Price et al, 2006; doi:10.1038/ng1847)
#' @importFrom stats cov
#' @export
eigstrat <- function(geno) {
  col_names <- colnames(geno)
  
  ## Center SNPs, set NA to 0
  geno <- apply(geno, 1, function(x) {x - mean(x, na.rm = TRUE)})
  geno[is.na(geno)] <- 0
  
  ## apply() will transpose matrix if nrows < ncols
  ## Usually won't happen with SNP data
  if(identical(rownames(geno), col_names)) {
    geno <- t(geno)
  }
  
  ## Normalize SNPs
  p <- as.vector(1 + rowSums(geno))/(2 + 2*ncol(geno))
  denom <- sqrt(p * (1 - p))
  geno <- t(t(geno) / denom)
  
  ## Calculate covariance matrix
  cov_mat <- cov(geno)
  
  ## Perform Eigen decomposition
  return(eigen(cov_mat))
}
