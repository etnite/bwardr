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
#' @param nvecs Either "all" to return all eigenvectors, or an integer to return
#'   the first n eigenvectors
#' @return A list with components
#'   \item{eigvecs}{Data frame containing the eigenvectors of the input
#'   genotypic matrix. The number of output eigenvectors is specified by the
#'   nvecs argument}
#'   \item{eigvals}{Vector of eigenvalues of the input genotypic matrix}
#' @details This function performs principle component analysis on a matrix of
#'   SNP data in the method of EIGENSTRAT (Price et al, 2006; doi:10.1038/ng1847).
#'   Note that while the user can select the number of eigenvectors to return,
#'   the function always returns all the eigenvalues.
#' @importFrom stats cov
#' @export
eigstrat <- function(geno, snps = "rows", nvecs = 25) {
  
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
  eig <- eigen(cov_mat)
  
  ## If nvecs is set to "all", give it a numeric val. equal to the number of
  ## eigenvectors
  if (nvecs == "all") {
    nvecs <- ncol(eig$vectors)
  }
  
  ## Set nvecs to number of eigenvectors if necessary
  if (nvecs > ncol(eig$vectors)) {
    nvecs = ncol(eig$vectors)
    warning(paste("nvecs has been reduced to the number of eigenvectors:", nvecs))
  }
  
  ## Create eigenvector data frame
  eigvecs <- as.data.frame(eig$vectors[, 1:nvecs])
  colnames(eigvecs) <- paste("PC", c(1:ncol(eigvecs)), sep = "")
  rownames(eigvecs) <- colnames(cov_mat)

  ## Assign eigenvalues to eigvals
  eigvals <- eig$values
  
  return(list("eigvecs" = eigvecs, "eigvals" = eigvals))  
}
