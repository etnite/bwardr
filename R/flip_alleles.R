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
  if (min(genmat, na.rm = TRUE) < 0 || max(genmat, na.rm = TRUE) > 2) {
    stop("Input genotypic matrix should be encoded {0, 1, 2}")
  }
  
  ## These lines first shift SNPs to (-1, 0, 1) to make use of rowSums(),
  ## Then shift SNPs to (0, 1, 2) or (0, -1, -2)
  ## depending on which allele is minor
  ## Then the absolute value is taken to normalize to 0 = major, 2 = minor
  genmat <- genmat - 1
  if (snps == "cols") {genmat <- t(genmat)}
  allele_count <- rowSums(genmat, na.rm = TRUE)
  genmat[allele_count <= 0, ] <- genmat[allele_count <= 0, ] + 1
  genmat[allele_count > 0, ] <- genmat[allele_count > 0, ] - 1
  genmat <- abs(genmat)
  if (snps == "cols") {genmat <- t(genmat)}
  
  return(genmat)
}