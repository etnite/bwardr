#' Generate files for QTL mapping
#' 
#' @param bed A binary ped (bed) object created with the package 'gaston',
#'   likely either using the functions gaston::read.bed.matrix() or
#'   gaston::read.vcf()
#' @param par_a String identifying the column in the SNP matrix containing
#'   data for the first parent.
#' @param par_b String identifying the column in the SNP matrix containing 
#'   data for the second parent.
#' @param rm_het Logical indicating whether to remove SNPs that are heterozygous
#'   in either parent.
#' @param rm_miss Logical indicating whether to remove SNPs that are missing
#'   in either parent.
#' @param include_pars Logical indicating whether to include the parents in the
#'   output data
#' @param out_fmt String indicating the format of the output table - one of
#'   either "icimapping" or "rqtl"
#' @return list with components
#'   \item{abh}{Dataframe in A/B/H format (see details)}
#'   \item{flipped}{Character vector consisting of IDs of SNPs with flipped alleles}
#' @details This function will always remove SNPs that are monomorphic between
#'   the two parents, SNPs with missing data in both parents, and SNPs with 
#'   with missing data in one parent and a heterozygous call in the other. The
#'   user can optionally specify whether to remove SNPs with missing data for
#'   either parent, or heterozygous calls in either parent. If "rqtl" output 
#'   format is selected, a table suitable for use as the genotypic portion of 
#'   the separate, rotated (i.e. "csvsr") R/qtl input format is produced.
#' @export
format_qtlmap_geno <- function(bed, 
                               par_a, par_b, 
                               rm_het = TRUE, rm_miss = TRUE,
                               include_pars = TRUE,
                               out_fmt = "rqtl") {
  
  ## Perform input sanity checks
  if (!(par_a %in% bed@ped$id)) {
    stop("Parent A is not present in supplied bed matrix")
  }
  
  if (!(par_b %in% bed@ped$id)) {
    stop("Parent B is not present in supplied bed matrix")
  }
    
  if (!(is.logical(rm_het) & 
        is.logical(rm_miss) & 
        is.logical(include_pars))) {
    stop("rm_het, rm_miss, and include_pars must be TRUE or FALSE")
  }
  
  if (!(out_fmt %in% c("icimapping", "rqtl"))) {
    stop("Please select either 'icimapping' or 'rqtl' for out_fmt")
  }
  
  ## Convert bed object to matrix
  genomat <- as.data.frame(t(gaston::as.matrix(bed)))
  
  ## Isolate the parents
  parents <- data.frame("snp" = rownames(genomat),
                        "a" = genomat[[par_a]],
                        "b" = genomat[[par_b]],
                        stringsAsFactors = FALSE)
  parents$sum <- parents$a + parents$b
  
  ### Select SNPs for removal ###
  ## Always select SNPs that are monomorphic between parents
  drop_snps <- NULL
  drop_snps <- union(drop_snps, parents$snp[parents$a == parents$b])
  
  ## Find SNPs with either parent het if specified; always remove SNPs with one
  ## parent het, the other missing
  if (rm_het) {
    drop_snps <- union(drop_snps, parents$snp[parents$a == 1 | parents$b == 1])
  } else {
    drop_snps <- union(drop_snps, parents$snp[is.na(parents$a) & parents$b == 1])
    drop_snps <- union(drop_snps, parents$snp[parents$a == 1 & is.na(parents$b)])
  }
  
  ## Find SNPs with either parent missing if specified; always remove SNPs with
  ## both parents missing
  if (rm_miss) {
    drop_snps <- union(drop_snps, parents$snp[is.na(parents$a) | is.na(parents$b)])
  } else {
    drop_snps <- union(drop_snps, parents$snp[is.na(parents$a) & is.na(parents$b)])
  }
  
  ## Perform the SNP removal
  genomat <- genomat[!rownames(genomat) %in% drop_snps, ]
  bed@snps <- bed@snps[!bed@snps$id %in% drop_snps, ]
  
  ### Flip SNPs ###
  
  ## The integer 0 is assigned as the parent A allele, and 2 as the parent B allele
  ## This requires "flipping" some SNPs (on average half of them)
  flip <- parents$snp[parents$a == 2]
  flip <- union(flip, parents$snp[is.na(parents$sum) & parents$b == 0])
  flip <- union(flip, parents$snp[parents$sum == 1 & parents$b == 0])
  flip <- flip[!is.na(flip)]
  genomat[flip, ] <- genomat[flip, ] - 2
  genomat <- abs(genomat)
  
  ## Remove parents if specified; otherwise reorder cols so parents appear first
  if (include_pars) {
    colorder <- setdiff(colnames(genomat), c(par_a, par_b))
    colorder <- c(par_a, par_b, colorder)
    genomat <- genomat[colorder]
  } else {
    genomat[c(par_a, par_b)] <- NULL
  }
  
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
  
  out_list <- list("abh" = out_df, "flipped" = flip)
  return(out_list)
}
