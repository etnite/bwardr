#' Perform PCA on marker data
#' 
#' @param geno A matrix of SNP data with SNPs in rows, genotypes in columns.
#'   SNPs should be encoded in minor-allele dosage format, so that:
#'   2 = homozygous for minor allele
#'   1 = heterozygous
#'   0 = homozygous for major allele
#'   NA = missing data
#' @return The Eigen decomposition of the input genotypic matrix
#' @details This function performs principle component analysis on a matrix of
#'   SNP data. Note that there are a few transpositions required, as R would
#'   typically not perform the centering/scaling/PCA in the correct orientation
#' @importFrom stats prcomp
#' @export
marker_pca <- function(geno) {
  col_names <- colnames(geno)
  
  ## Center/scale SNPs, set NA to 0
  geno <- t(scale(t(geno)))
  geno[is.na(geno)] <- 0
  
  ## Perform the PCA
  pca <- prcomp(geno, center = FALSE, scale. = FALSE)
  eigvecs <- pca$x
  eigvals <- pca$sdev^2
  
  return(list("vectors" = eigvecs, "values" = eigvals))
}