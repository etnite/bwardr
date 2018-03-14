#' Perform PCA on marker data
#' 
#' @param geno A matrix of SNP data with SNPs in minor-allele dosage format:
#'   2 = homozygous for minor allele
#'   1 = heterozygous
#'   0 = homozygous for major allele
#'   NA = missing data
#' @param snps Either "rows" to indicate that SNPs are in rows and individuals 
#'   are in columns, or "cols" to indicate the SNPs are in columns and
#'   individuals are in rows
#' @return The Eigen decomposition of the input genotypic matrix
#' @details This function performs a "quick and dirty" principle component 
#'   analysis on the input matrix of SNP data. A genetic relationship matrix
#'   is calculated first, and then an Eigen decomposition is performed in this
#'   matrix
#' @export
marker_pca <- function(geno, snps = "rows") {

  ## Check matrix orientation
  if (! snps %in% c("rows", "cols")) {
    stop("Please specify whether SNPs are in 'rows' or 'cols'")
  }
  if (snps == "cols") { geno <- t(geno)}
  
  ## Fast genetic relationship matrix generation
  R <- (1/nrow(geno)) * (t(geno) %*% geno)
  
  ## Perform eigen decomposition and return - that's it
  return(eigen(R))
}