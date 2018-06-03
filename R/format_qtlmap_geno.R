#' Generate files for QTL mapping
#' 
#' @param bed A binary ped (bed) object created with the package 'gaston',
#'   likely either using the functions gaston::read.bed.matrix() or
#'   gaston::read.vcf()
#' @param par_a String identifying the column in the SNP matrix containing
#'   data for the first parent.
#' @param par_b String identifying the column in the SNP matrix containing 
#'   data for the second parent.
#' @param rm_mono Logical indicating whether to remove SNPs that are monomorphic
#'   between the two parents.
#' @param rm_het Logical indicating whether to remove SNPs that are heterozygous
#'   in either parent.
#' @param rm_miss Logical indicating whether to remove SNPs that are missing
#'   in either parent.
#' @param include_pars Logical indicating whether to include the parents in the
#'   output data
#' @param out_fmt String indicating the format of the output table - one of
#'   either "icimapping" or "rqtl"
#' @details If "rqtl" output format is selected, a table suitable for use as the
#'   genotypic portion of the separate, rotated (i.e. "csvsr") R/qtl input
#'   format is produced.
#' @export
format_qtlmap_geno <- function(bed, 
                               par_a, par_b, 
                               rm_mono = TRUE, rm_het = TRUE, rm_miss = TRUE,
                               include_pars = TRUE,
                               out_fmt = "rqtl") {
  
  ## Perform input sanity checks
  if (!(par_a %in% bed@ped$id)) {
    stop("Parent A is not present in supplied bed matrix")
  }
  
  if (!(par_b %in% bed@ped$id)) {
    stop("Parent B is not present in supplied bed matrix")
  }
    
  if (!(is.logical(rm_mono) & 
        is.logical(rm_het) & 
        is.logical(rm_miss) & 
        is.logical(include_pars))) {
    stop("rm_mono, rm_het, rm_miss, and include_pars must be TRUE or FALSE")
  }
  
  if (!(out_fmt %in% c("icimapping", "rqtl"))) {
    stop("Please select either 'icimapping' or 'rqtl' for out_fmt")
  }
  
  ## The integer 0 is assigned as the parent A allele, and 2 as the parent B allele
  ## This requires "flipping" some SNPs
  ## After this, parent B should only have the 0 allele when a SNP is monomorphic
  genomat <- as.data.frame(t(gaston::as.matrix(bed)))
  flip <- rownames(genomat)[genomat[[par_a]] == 2]
  genomat[flip, ] <- genomat[flip, ] - 2
  genomat <- abs(genomat)
  
  ## Flip the alleles in the bed object's snp info for good measure
  a2 <- data.frame("id" = bed@snps$id, "A2" = bed@snps$A2, stringsAsFactors = FALSE)
  bed@snps$A2[bed@snps$id %in% flip] <- bed@snps$A1[bed@snps$id %in% flip]
  bed@snps$A1[bed@snps$id %in% flip] <- a2$A2[a2$id %in% flip]
  
  ## Reorder genotypes so that the two parents come first
  colorder <- setdiff(colnames(genomat), c(par_a, par_b))
  colorder <- c(par_a, par_b, colorder)
  genomat <- genomat[colorder]
  
  ## Isolate the parents
  parents <- data.frame("snp" = rownames(genomat),
                        "a" = genomat[[par_a]],
                        "b" = genomat[[par_b]],
                        stringsAsFactors = FALSE)
  
  ## Find monomorphic SNPs if specified
  drop_snps <- NULL
  if (rm_mono) {
    drop_snps <- union(drop_snps, parents$snp[parents$a == parents$b])
  }
  
  ## Find SNPs with either parent het if specified
  if (rm_het) {
    drop_snps <- union(drop_snps, parents$snp[parents$a == 1 | parents$b == 1])
  }
  
  ## Find SNPs with either parent missing if specified
  if (rm_miss) {
    drop_snps <- union(drop_snps, parents$snp[is.na(parents$a) | is.na(parents$b)])
  }
  
  ## Remove SNPs
  genomat <- genomat[!rownames(genomat) %in% drop_snps, ]
  bed@snps <- bed@snps[!bed@snps$id %in% drop_snps, ]
  
  ## Format for either ICIMapping or R/qtl output
  genomat$id <- rownames(genomat)
  if (out_fmt == "icimapping") {
      genomat[is.na(genomat)] <- -1
      id_df <- bed@snps[c("id", "chr", "dist", "pos", "A1", "A2")]
      out_df <- merge(id_df, genomat, by = "id")
  } else {
      genomat[genomat == 0] <- "A"
      genomat[genomat == 2] <- "B"
      genomat[genomat == 1] <- "H"
      genomat[is.na(genomat)] <- "-"
      id_df <- bed@snps[c("id", "chr", "dist")]
      out_df <- merge(id_df, genomat, by = "id")
      names(out_df)[2:3] <- ""
  }
  
  ## Remove parents if specified
  if (!include_pars) {
    out_df[c(par_a, par_b)] <- NULL
  }
  
  return(out_df)
}
