#' Generate files for QTL mapping
#' 
#' @param genomat A matrix or dataframe of SNP data with SNPs in rows and genotypes in
#'   columns, encoded {0,1,2,NA}, where 0 indicates the homozygous state for one 
#'   allele, 2 indicates the homozygous state for the alternate allele, 1
#'   indicates the heterozygous state, and NA indicates missing data. The alleles
#'   encoded by 0 and 2 for each SNP can be arbitrary. Row names should be
#'   SNP identifiers, and column names should be genotype identifiers.
#' @param map Dataframe containing map information for SNPs in genomat. This 
#'   should contain the columns "SNP", "CHR", and "DIST" for the SNP ID, 
#'   chromosome, and genetic distance (if known, or else all 0), respectively
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
#' @param out_fmt String indicating the format of the output table - one of
#'   either "icimapping" or "rqtl"
#' @details It is strongly recommended that users utilize the read.vcf() function
#'   in package 'gaston' to create a binary ped object from an input VCF genotypic
#'   file. This binary ped object can then easily be converted into the required
#'   input matrix using the command t(as.matrix(bed)). If "rqtl" is selected for
#'   the output format, the function will output a table suitable for import to
#'   R/qtl as the genotypic portion of the "csvs" format (which uses separate
#'   genotypic and phenotypic tables)
format_qtlmap_geno <- function(genomat, map, 
                               par_a, par_b, 
                               rm_mono = TRUE, rm_het = TRUE, rm_miss = TRUE, 
                               out_fmt = "rqtl") {
  
  ## Sanity checks
  if (!out_fmt %in% c("icimapping", "rqtl")) {
    stop("Please select either 'icimapping' or 'rqtl' for out_fmt")
  }
  
  if(colnames(map) != c("SNP", "CHR", "DIST")) {
    stop("map dataframe should have column names 'SNP', 'CHR' and 'DIST'")
  }
  
  ## First, 2 is assigned as the parent A allele, and 0 as the parent B allele
  ## This requires "flipping" some SNPs
  ## After this, parent B should only have the 2 allele when SNP is monomorphic
  genomat <- as.data.frame(genomat)
  para_vec <- genomat[[par_a]]
  genomat[para_vec == 0, ] <- genomat[para_vec == 0, ] - 2
  genomat <- abs(genomat)
  
  ## Reorder genotypes so that the two parents come first
  colorder <- setdiff(colnames(genomat), c(par_a, par_b))
  colorder <- c(par_a, par_b, colorder)
  genomat <- genomat[colorder]
  
  ## Remove monomorphic SNPs if specified
  if (rm_mono) {
    genomat <- genomat[genomat[[par_a]] != genomat[[par_b]], ]
  }
  
  ## Remove SNPs with either parent het if specified
  if (rm_het) {
    genomat <- genomat[genomat[[par_a]] != 1 & genomat[[par_b]] != 1, ]
  }
  
  ## Remove SNPs with either parent missing if specified
  if (rm_miss) {
    genomat <- genomat[!is.na(genomat[[par_a]]) & !is.na(genomat[[par_b]]), ]
  }
  
  ## Formatting for either ICIMapping or R/qtl output
  if (out_fmt == "icimapping") {
      genomat[is.na(genomat)] <- -1
  } else {
      genomat[genomat == 2] <- "A"
      genomat[genomat == 0] <- "B"
      genomat[genomat == 1] <- "H"
  }
  
  ## Merge in map data
  genomat$SNP <- rownames(genomat)
  genomat <- merge(map, genomat, by = "SNP")
  
  ## Additional formatting required for R/qtl output
  if (out_fmt == "rqtl") {
    rownames(genomat) <- genomat$SNP
    drops <- c("SNP", "POS", "A1", "A2")
    genomat <- genomat[, !(names(genomat) %in% drops)]
    genomat <- as.data.frame(t(genomat), stringsAsFactors = FALSE)
    genomat <- data.frame("id" = rownames(genomat), genomat, stringsAsFactors = FALSE)
    genomat[1:2, 1] <- ""
  }
  
  return(genomat)
}