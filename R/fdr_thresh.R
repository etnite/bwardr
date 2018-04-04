#' FDR threshold and SNP highlighting
#'
#' @param gwas_results A dataframe containing at least the following columns:
#'   1) SNP identifiers
#'   2) GWAS p-values.  These should NOT be -log transformed
#' @param id String specifying the name of the column in gwas_results
#'   containing SNP identifiers
#' @param p Character specifying the name of the column in gwas_results
#'   containing GWAS output p-values
#' @param meth String specifying the method of p-value adjustment to use - 
#'   see output of help(p.adjust) for details
#' @param thresh Number between 0 and 1 - the desired FDR threshold
#' @return list consisting of
#'   \item{p.thresh}{approximation of the FDR threshold for p-value distribution}
#'   \item{sig.snps}{vector of names of SNPs passing the FDR threshold}
#' @details Adjusted p-values are calculated from a vector of p-values output
#'   from a GWAS analysis. If any adjusted p-values fall below the specified
#'   threshold, an approximate threshold for adjusted p-values is returned, with
#'   a list of the IDs of significant SNPs
#' @importFrom stats lm p.adjust
#' @export
fdr_thresh <- function(gwas_results, id, p, meth = "fdr", thresh = 0.05) {
  
  ## Stop if user-defined FDR threshold is outside valid range
  if(thresh <= 0 || thresh >= 1) {
    stop("Please enter a value between 0 and 1 for the FDR threshold")
  }
  
  ## Calculate adjusted p-values.  If no adjusted p-values fall below the
  ## FDR threshold, then set p.thresh and sig.snps to NULL
  p_adj <- p.adjust(gwas_results[[p]], method = meth)
  if(min(p_adj, na.rm = TRUE) > thresh) {
      cat("No SNP adjusted p-values exceed the FDR threshold of", thresh, "\n")
      p_thresh <- NA
      sig_snps <- NULL
  } else {
      gwas_results$p_adj <- p_adj
      gwas_results <- gwas_results[order(gwas_results[[p]]), ]
      
      ## Isolate SNPs with "significant" adjusted p-values
      sig_subset <- gwas_results[gwas_results$p_adj <= thresh, ]
      sig_subset <- droplevels(sig_subset)
      sig_snps <- sig_subset[[id]]
      
      ## Isolate SNPs with "non-significant" adjusted p-values
      nonsig_subset <- gwas_results[gwas_results$p_adj > thresh, ]
      nonsig_subset <- droplevels(nonsig_subset)
      
      ## Isolate the significant SNP with the highest p-value
      ## and the non-significant SNP with the lowest p-value
      bracket <- rbind(sig_subset[nrow(sig_subset), ],
                       nonsig_subset[1, ])
      
      ## Draw a straight line between the two SNPs in an [adjusted p-value, p-value] 2D space
      ## Then plug the FDR threshold into the regression to find the unadjusted p-value
      ## corresponding to the threshold
      conn_points <- lm(bracket[[p]] ~ bracket$p_adj)
      p_thresh <- conn_points$coefficients[2]*thresh +
        conn_points$coefficients[1]
  }  
  
  fdr_out <- list("p_thresh" = as.numeric(p_thresh), "sig_snps" = sig_snps)
  return(fdr_out)  
}