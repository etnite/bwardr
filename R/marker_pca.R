#' Perform PCA on marker data
#' 
#' @param geno A matrix of SNP data with SNPs in rows, genotypes in columns.
#'   SNPs should be encoded in minor-allele dosage format, so that:
#'   2 = homozygous for minor allele
#'   1 = heterozygous
#'   0 = homozygous for major allele
#'   NA = missing data
#' @param snps Either "rows" to indicate that SNPs are in rows and individuals 
#'   are in columns, or "cols" to indicate the SNPs are in columns and
#'   individuals are in rows
#' @return The Eigen decomposition of the input genotypic matrix
#' @details This function performs principle component analysis on a matrix of
#'   SNP data. Note that there are a few transpositions required, as R would
#'   typically not perform the centering/scaling/PCA in the correct orientation
#' @importFrom stats prcomp
#' @export
marker_pca <- function(geno, snps = "rows") {

  ## Check matrix orientation
  if (! snps %in% c("rows", "cols")) {
    stop("Please specify whether SNPs are in 'rows' or 'cols'")
  }
  if (snps == "cols") { geno <- t(geno)}
  
  ## Center/scale SNPs, set NA to 0
  R <- (1/nrow(geno)) * (t(geno) %*% geno)
  
  ## Perform eigen decomposition and return - that's it
  return(eigen(R))
}