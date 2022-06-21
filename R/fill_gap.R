#' Fill in a gap in recombination distance
#' 
#' @param df - A dataframe containing info on SNPs for a single gap region. Must
#'   contain the column "dist" for each SNP's recombination distance (in cM), 
#'   but otherwise can contain any additional number of columns
#' @param start Numeric - the start position of the gap in cM (only used when method
#'   is set to `value`)
#' @param end Numeric - the end position of the gap in cM (only used when method
#'   is set to `value`)
#' @param n_snps Integer - the ideal number of SNPs to select from the region.
#'   The actual number of SNPs selected may be lower than this value
#'   when setting method to `value`
#' @param method String - either one of `value` or `percentile`
#' @details This function will attempt to find SNPs to fill in gaps (in terms
#'   of recombination distance) on chromosomes. To do so, it requires a dataframe
#'   containing the cM coordinates of SNPs within the gap region, 
#'   and some number supplied for `n_snps`, which essentially specifies the maximum number of SNPs
#'   to return. Note that the gap region should contain more SNPs than `n_snps`. 
#'   The function then generates "virtual SNPs" which are evenly spaced across
#'   the gap region. These can be spaced according to actual recombination distance
#'   (cM) or else percentile. If using the former, the function will typically
#'   return fewer SNPs than `n_snps` because the actual SNPs present in the region
#'   are likely not distributed uniformly. Therefore multiple virtual SNPs may
#'   correspond to a single actual SNP. In the case of using percentiles, the
#'   function will return n_snps unless the region actually contains fewer SNPs
#'   than `n_snps`. 
#'
#'   This function can be used to find SNPs that are approximately equally spaced
#'   across a chromosome by considering the entire chromosome as a gap. If using
#'   the `value` method, simply set the start value to 0 and the end value to the
#'   expected length of the chromosome in cM.
#' @return The input data.table with added weights column and rows filtered to 
#'   only include selected SNPs. The weights column represents clustering of SNPs - 
#'   i.e. the number of SNPs that are "respresented" by the selected SNP.
#' @importFrom stats ecdf
#' @importFrom utils head tail
#' @export
fill_gap <- function(df, start, end, n_snps, method = "value") {
  
  ## Set up the vector of ideally spaced virtual SNPs, and the vector of what
  ## we actually have
  if (method == "value") {
    ideal_vec <- seq(start, end, length.out = n_snps + 2)
    actual_vec <- df$dist
  } else if (method == "percentile") {
    ideal_vec <- seq(0, 1, length.out = n_snps + 2)
    actual_vec <- ecdf(df$dist)(df$dist)
  } else {
    stop("Please select either 'value' or 'percentile' for method")
  }
  ideal_vec <- head(ideal_vec, -1)
  ideal_vec <- tail(ideal_vec, -1)
  
  ## Find closest actual SNP for each ideal virtual SNP
  ind <- vector("integer", length = length(ideal_vec))
  for (i in 1:length(ind)) {
    ind[i] <- which.min(abs(actual_vec - ideal_vec[i]))
  }
  ind_uniq <- unique(ind)
  weights <- as.vector(table(ind))
  
  ## Subset input data.table and add weights column
  df <- df[ind_uniq, ]
  df$weight <- weights
  return(df)
}