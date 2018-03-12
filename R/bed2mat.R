#' Subset a BED file and convert to SNP matrix
#'
#' @param bed A bed.matrix object created with package 'gaston'
#' @param genos A character vector listing genotypes (i.e. individuals) to
#'   include in the final genotypic matrix
#' @param min_maf Number between 0 and 0.5 - the minimum minor-allele-frequency
#'   threshold to apply after subsetting
#' @param max_het Number between 0 and 1 - the maximum SNPwise heterozygosity to
#'   allow after subsetting
#' @return A subsetted, possibly filtered SNP matrix created from the input
#'   BED file
#' @details This function takes a BED (PLINK binary pedigree/map) file,
#'   optionally subsets individuals from it, and optionally filters SNPs by minor-
#'   allele frequency and/or percent heterozygous calls. It then outputs a SNP 
#'   matrix. It is advised to re-filter SNPs whenever subsetting individuals, 
#'   as allele frequencies and the proportion of heterozygous calls may change.
#' @export
bed2mat <- function(bed, genos=NULL, min_maf=0, max_het=1) {
  
  ## Perform optional subsetting/filtering
  if (!is.null(genos)) { bed <- gaston::select.inds(bed, id %in% genos) }
  bed <- gaston::select.snps(bed, maf >= min_maf)
  bed <- gaston::select.snps(bed, hz <= max_het)
  bed <- gaston::set.stats(bed, verbose = FALSE)
  
  ## "Flip" SNPs for which 2 encodes the major allele after subset/filter
  gen_mat <- gaston::as.matrix(bed)
  gen_mat[, bed@p > 0.5] <- gen_mat[, bed@p > 0.5] - 2
  gen_mat <- abs(gen_mat)
  
  return(gen_mat)
}