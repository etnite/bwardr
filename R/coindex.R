#' Calculate Coincidence Index
#' 
#' @param obs A named vector of observed phenotypic data
#' @param pred A named vector of predicted phenotypic data
#' @param si Selection index - value from 0 to 1
#' @param best String or number signifying what phenotypic values are considered
#'   "best." Must be one of the following:
#'  \enumerate{
#'    \item "high" to perform directional selection of genotypes with highest trait values
#'    \item "low" to perform directional selection of genotypes with lowest trait values
#'    \item "mean" to perform stabilizing selection around the mean values of the observed and predicted data
#'    \item "median" to perform stabilizing selection around the median values of the observed and predicted data
#'    \item numeric, in which case it is assumed that the ideal phenotypic value has been supplied
#'  }
#' @param naive Logical - if TRUE then a simple proportion matching coefficient is
#'   returned, without correcting for random chance
#' @return The calculated coincidence index
#' @details This function calculates the coincidence index as defined by Hamblin
#'   and Zimmerman, 1986. (https://doi.org/10.1002/9781118061015.ch8). In the context
#'   of this function, its intended use is to compare direct phenotypic selection
#'   with indirect selection based upon genomic data. However, it can be used to
#'   compare selections made using any two different methods. The formula for the
#'   coincidence index is
#'   \deqn{CI = (C - R) / (T - R)}
#'   where T is the number of observed genotypes selected using the selection index,
#'   R is the expected number of genotypes among these correctly selected due to 
#'   chance, and C is the number of coincident genotypes selected using both methods. 
#'   
#'   For instance, if we have a population of 100 genotypes and a selection index of
#'   0.1, T = 100 * 0.1 = 10 and R = 10 * 0.1 = 1. The value of C will depend upon
#'   how many genotypes are selected based upon both observed and predicted data.
#'   
#'   A value of 1 indicates identical selections made using the observed and
#'   predicted phenotypic values, while a value of 0 indicates that selections
#'   based on predicted values were no better than random chance. Note that CI values
#'   can be negative. Small samples sizes will lead to high CI variance.
#' @export
coindex <- function(obs, pred, si, best = "high", naive = FALSE) {
  
  ## Sanity check on inputs
  if (si >= 1 | si <= 0) {
    stop("Selection index must be between 0 and 1")
  }
  
  ## Set the target value for selection, based on best being specified as either
  ## "high", "low", "mean", "median" or a numeric value
  if (best == "high") {
    obs_targ <- max(obs)
    pred_targ <- max(pred)
  } else if (best == "low") {
    obs_targ <- min(obs)
    pred_targ <- min(pred)
  } else if (best == "mean") {
    obs_targ <- mean(obs)
    pred_targ <- mean(pred)
  } else if (best == "median") {
    obs_targ <- median(obs)
    pred_targ <- median(pred)
  } else if (is.numeric(best)) {
    obs_targ <- best
    pred_targ <- best
  } else {
    stop("Incorrect specification of 'best' parameter")
  }
  
  ## Check the intersection between the two input vectors
  ## Throw warning if mismatched names are detected
  inter_names = sort(intersect(names(obs), names(pred)))
  union_names = sort(union(names(obs), names(pred)))
  if (!identical(inter_names, union_names)) {
    warning("Mismatches between the two input vectors. The following elements are being discarded:")
    print(setdiff(union_names, inter_names))
  }
  obs <- obs[inter_names]
  pred <- pred[inter_names]
  
  ## Calculate t and r
  t <- round(length(obs) * si)
  r <- round(t * si)
  
  ## Make selections based on observed and predicted vals; calculate c
  obs <- sort(abs(obs_targ - obs))[1:t]
  pred <- sort(abs(pred_targ - pred))[1:t]
  c <- length(intersect(names(obs), names(pred)))
  
  ## Calculate and return coincidence index
  if (naive) {
    CI <- c / t
  } else {
    CI <- (c - r) / (t - r)
  }
  return(CI)
}