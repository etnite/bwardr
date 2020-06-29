#' Calculate Generalized Heritability from lme4 Model
#' 
#' @param model A lme4 model object
#' @param geno_name A string denoting the label of the random genotypic effect in the
#'   supplied lme4 model object
#' @return A list containing the following elements:
#' * SED A vector containing the standard errors of the genotypic random effects
#'   (the length of this vector is the number of genotypes supplied to the model)
#' * H2 The generalized heritability estimate
#' @details This function calculates generalized heritability using the method of
#'   Cullis et al., 2006 (\url{https://doi.org/10.1198/108571106X154443}). This
#'   method can be used in unbalanced applications where the traditional entry
#'   mean heritability calculation will give biased estimates. The method of doing
#'   this using lme4 is detailed at \url{https://shantel-martinez.github.io/resources.html}
#' @export
#' @md
Cullis_H2 <- function(model, geno_name = "GENO") {
  
  ses <- arm::se.ranef(model)[[geno_name]]
  v_BLUP <- ses^2
  var_g <- lme4::VarCorr(model, comp="Variance")[[geno_name]][1]
  reliability <- 1 - (v_BLUP / (2 * var_g))
  H2 <- round(mean(reliability), 3)
  
  out_list <- list("SE" = ses, "H2" = H2)
  return(out_list)
}