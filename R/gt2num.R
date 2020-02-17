#' Convert VCF GT fields to numeric format
#' 
#' @param genomat A matrix of VCF-style genotype calls, with SNPs in rows, and
#'   samples/individuals in columns
#' @return A list containing the following elements:
#' * genomat Converted numeric genotype matrix
#' * removed_loci Vector of loci removed due to having greater than two alleles
#' @details This function will convert a matrix of genotype calls from a VCF file
#'   into numeric format. It will convert homozygous reference allele calls to 0,
#'   heterozygous calls to 1, and homozygous alternate allele calls to 2.
#'   Note that the output IS NOT necessarily in minor-allele dosage format. 
#'   See ?flip_alleles for more details. Loci with more than two alleles will 
#'   be removed and reported. An input matrix of the required format can be 
#'   generated with vcfR::extract.gt(), with bcftools query, or manually in R.
#' @seealso flip_alleles
#' @md
#' @export
gt2num <- function(genomat) {
  
  if (!is.matrix(genomat)) {
    stop("genomat parameter must be a matrix")
  }

  ## Define genotypes
  refs <- c("0/0", "0|0")
  alts <- c("1/1", "1|1")
  hets <- c("0/1", "1/0", "0|1", "1|0")
  misses <- c("./.", ".|.")
  
  ## Perform numeric substitution
  genomat[genomat %in% refs] <- 0
  genomat[genomat %in% alts] <- 2
  genomat[genomat %in% hets] <- 1
  genomat[genomat %in% misses] <- NA
  
  ## Record loci with > 2 alleles; remove from matrix
  biallelic <- apply(genomat, 1, function(x) { all(x %in% c(NA, "0", "1", "2")) })
  if (all(biallelic)) {
    removed_loci <- NULL
  } else {
    removed_loci <- rownames(genomat)[!biallelic]
    genomat <- genomat[biallelic, ]
  }
  
  ## Construct and return output list
  class(genomat) <- "numeric"
  return(list("genomat" = genomat, "removed_loci" = removed_loci))
}