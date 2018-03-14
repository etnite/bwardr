#' Perform Eigenstrat marker PCA
#' 
#' @param geno A matrix of SNP data encoded in minor-allele dosage format:
#'   2 = homozygous for minor allele
#'   1 = heterozygous
#'   0 = homozygous for major allele
#'   NA = missing data
#' @param snps Either "rows" to indicate that SNPs are in rows and individuals 
#'   are in columns, or "cols" to indicate the SNPs are in columns and
#'   individuals are in rows
#' @return The Eigen decomposition of the input genotypic matrix
#' @details This function performs principle component analysis on a matrix of
#'   SNP data in the method of EIGENSTRAT (Price et al, 2006; doi:10.1038/ng1847)
#' @importFrom stats cov
#' @export
eigstrat <- function(geno, snps = "rows") {
  
  ## Check matrix orientation
  if (! snps %in% c("rows", "cols")) {
    stop("Please specify whether SNPs are in 'rows' or 'cols'")
  }
  if (snps == "rows") { geno <- t(geno)}
  
  ## Center SNPs
  geno <- scale(geno, scale = FALSE)
  
  ## Normalize SNPs
  p_est <- as.vector(1 + colSums(geno, na.rm = TRUE))/(2 + 2*nrow(geno))
  denom <- sqrt(p_est * (1 - p_est))
  geno <- geno/denom
  
  ## Set NAs to 0; Calculate covariance matrix
  geno[is.na(geno)] <- 0
  cov_mat <- cov(t(geno))
  
  ## Perform Eigen decomposition
  return(eigen(cov_mat))
}
