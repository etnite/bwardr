#' Convert a VCF file to minor-allele dosage format
#'
#' @param vcf_dat - A vcfR object created by the read.vcfR() function of package
#'   vcfR
#' @return dataframe in the same format as the data rows of the input VCF data, 
#'   except with genotype calls transformed into minor allele dosage format 
#'   (i.e. homozygous for minor allele = 2, het = 1, homo. for major allele = 0)
#'   NOTE: VCF 'FORMAT' column is not present in output, since only the GT field
#'   is extracted from the input VCF
#' @details vcf2dos() uses vectorized operations and is therefore quite fast.
#'   However, since R must read an entire file into memory, large VCF files
#'   will cause problems, especially if the size of the file is similar to
#'   the amount of physical memory on the machine. In the case of large VCF 
#'   files, it may be best to split them apart by chromosome, or else use 
#'   command-line tools such as bcftools, vcftools, or PLINK 1.9
#' @export
vcf2dos <- function(vcf_dat) {
  gt <- vcfR::extract.gt(vcf_dat)
  
  ## Replace VCF geno calls with numeric values
  gt <- plyr::mapvalues(gt,
                        from = c("1/1", "1|1", "0/0", "0|0", "0/1", "1/0", "0|1", "1|0", "./.", ".|."),
                        to = c(1, 1, -1, -1, 0, 0, 0, 0, NA, NA),
                        warn_missing = FALSE)
  
  class(gt) <- "numeric"
  
  ## These lines shift SNPs to (0, 1, 2) or (0, -1, -2)
  ## depending on which allele is minor
  ## Then the absolute value is taken to normalize to 0 = major, 2 = minor
  allele_count <- rowSums(gt, na.rm = TRUE)
  gt[allele_count <= 0, ] <- gt[allele_count <= 0, ] + 1
  gt[allele_count > 0, ] <- gt[allele_count > 0, ] - 1
  gt <- abs(gt)
  
  return(data.frame(vcf_dat@fix, gt,
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    row.names = NULL))
}
