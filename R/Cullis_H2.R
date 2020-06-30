#' Calculate Generalized Heritability from lme4 Model
#' 
#' @param model A lme4 model object
#' @param geno_label A string denoting the label of the random genotypic effect in the
#'   supplied lme4 model object
#' @return A list containing the following elements:
#' * avsed The average standard error of differences between adjusted means estimates
#' * H2 The generalized heritability estimate
#' @details This function calculates generalized heritability using the method of
#'   Cullis et al., 2006 (\url{https://doi.org/10.1198/108571106X154443}).
#'   Specifically, their formula is H2 = 1 - (vblup / (2 * var_g)). Where the
#'   generalized heritability (H2) is a function of the reliability of the BLUPs
#'   (vblup - the average standard error of differences between BLUPs squared),
#'   and the genotypic variance (var_g). This method can be used in unbalanced 
#'   applications where the traditional entry-mean heritability calculation will 
#'   give biased estimates. The method of doing this using lme4 is detailed by 
#'   Ben Bolker at \url{https://stackoverflow.com/questions/38697477/mean-variance-of-a-difference-of-blues-or-blups-in-lme4}.
#'   This method yields values that are slightly different (I have observed up to
#'   0.75%) from ASReml-R's results. Another solution I came across at 
#'   \url{https://shantel-martinez.github.io/resources.html} seems to produce
#'   results that are more divergent from ASReml-R's.
#' @export
#' @md
Cullis_H2 <- function(model, geno_label = "GENO") {
  
  ## Extract genotypic variance
  var_g <- lme4::VarCorr(model, comp = "Variance")[[geno_label]][1]
  
  ## Extract genotypic conditional variances
  convars <- lme4::ranef(model, condVar = TRUE)
  g_convar <- attr(convars[[geno_label]], "postVar")
  
  ## Calculate VBLUP and avsed
  vblup <- 2 * mean(g_convar)
  avsed <- sqrt(vblup)
  
  ## Calculate generalized heritability
  H2 <- 1 - (vblup / (2 * var_g))
  
  ## Create and return output list
  out_list <- list("avsed" = avsed, "H2" = H2)
  return(out_list)
}