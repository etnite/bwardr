#' Flip alleles in a minor-allele dosage matrix
#' 
#' @param genmat A numeric matrix of SNPs x individuals, with SNPs encoded
#'   as {0, 1, 2}. Can be in either orientation (SNPs x individuals or
#'   individuals x SNPs)
#' @param snps Either "rows" to indicate that SNPs are in rows and individuals 
#'   are in columns, or "cols" to indicate the SNPs are in columns and
#'   individuals are in rows
#' @return Numeric matrix of same dimensions as input matrix, but with all SNPs
#'   encoded in proper minor-allele dosage format
#' @details SNP data is often stored in numeric matrices in minor-allele dosage
#'   format, where 0 indicates no copies of the minor allele (i.e. homozygous
#'   for the major allele), 1 indicates 1 copy of the minor allele (i.e.
#'   heterozygous), and 2 indicates 2 copies of the minor allele (i.e.
#'   homozygous for the minor allele). However, there are situations in which
#'   the alleles get "flipped" between major vs. minor (e.g. when subsetting
#'   or filtering genotypic data). flip_alleles() is intended to selectively
#'   reset SNPs which have become flipped to the proper minor-allele dosage
#'   orientation.
#' @export
flip_alleles <- function(genmat, snps = "rows") {
  
  ## Input sanity checks
  if (! snps %in% c("rows", "cols")) {
    stop("Please specify either 'rows' or 'cols' for snps argument")
  }
  if (!all(genmat %in% c(0, 1, 2))) {
    stop("Input genotypic matrix should be encoded {0, 1, 2}")
  }
  
  ## These lines transpose matrix if necessary
  ## Then, for any row (i.e. SNP) with mean > 1, 2 is subtracted, and
  ## the absolute value is taken
  if (snps == "cols") {genmat <- t(genmat)}
  snp_means <- rowMeans(genmat, na.rm = TRUE)
  genmat[snp_means > 1, ] <- abs(genmat[snp_means > 1, ] - 2)
  if (snps == "cols") {genmat <- t(genmat)}
  
  return(genmat)
}
