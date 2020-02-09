#' Flip alleles in a minor-allele dosage matrix
#' 
#' @param genmat A numeric matrix of SNPs x individuals, with SNPs encoded
#'   as {0, 1, 2}. Can be in either orientation (SNPs x individuals or
#'   individuals x SNPs)
#' @param snps Either "rows" to indicate that SNPs are in rows and individuals 
#'   are in columns, or "cols" to indicate the SNPs are in columns and
#'   individuals are in rows
#' @return List containing the following elements:
#' * dosmat - NumericSNP  matrix in minor-allele dosage format
#' * flipped_loci - Vector listing markers for which alleles were flipped
#' @details SNP data is often stored in numeric matrices in minor-allele dosage
#'   format, where 0 indicates no copies of the minor allele (i.e. homozygous
#'   for the major allele), 1 indicates 1 copy of the minor allele (i.e.
#'   heterozygous), and 2 indicates 2 copies of the minor allele (i.e.
#'   homozygous for the minor allele). However, there are situations in which
#'   the alleles get "flipped" between major vs. minor (e.g. when subsetting
#'   or filtering genotypic data). In addition, SNPs may be numerically encoded
#'   in various ways, e.g. by encoding the allele in the reference genome as
#'   0, and the alternate allele as 2. flip_alleles() is intended to selectively 
#'   "flip" SNPs into the proper minor-allele dosage orientation.
#' @export
#' @seealso gt2num
#' @md
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
  flip <- snp_means > 1
  
  if (sum(flip) == 0) {
    flipped_loci <- NULL
  } else {
    flipped_loci <- rownames(genmat)[flip]
    genmat[flip, ] <- abs(genmat[flip, ] - 2)
  }
  
  ## Transpose back if necessary
  if (snps == "cols") {genmat <- t(genmat)}
  
  ## Construct output list and return
  out_list <- list("dosmat" = genmat, "flipped_loci" = flipped_loci)
  return(out_list)
}
