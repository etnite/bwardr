#' Monte Carlo training/validation splits for multi-environment data
#' 
#' @param pheno Dataframe containing phenotypic data. The dataframe must 
#'   contain genotype IDs in a column named "IID", and environment designators 
#'   in a column named "ENV". The data can be replicated within environment, and
#'   can have any additional number of columns.
#' @param prop_val Number between 0 and 1 indicating the proportion of phenotypic
#'   data that should be set to the validation set.
#' @param cv_scheme Character string consisting of either "CV1" or "CV2". "CV1"
#'   assigns validation data by genotype (i.e. simulates introducing new genotypes
#'   which have not been tested in any environment). CV2 assigns validation data
#'   by genotype-environment combination (i.e. simulates introducing genotypes into
#'   new environments).
#' @return A dataframe identical to the input phenotypic dataframe, except with
#'   a TRAIN_VAL column added (if not previously present) or else re-randomized 
#'   (if the columns was already present) to indicate training/validation sets.
#' @details This function only performs Monte Carlo (i.e. random subsampling)
#'   cross-validation training/validation assignment. Note that k-fold cross
#'   validation becomes difficult to perform for a CV2 scheme. 
#'   If the CV2 cross-validation scheme is selected, the function will
#'   attempt to choose genotype-environment combinations to assign to the
#'   validation set such that any given genotype will only be assigned to the
#'   validation set in at most one environment. However, if the user sets a
#'   prop_val value higher than (1 / # envs), then some genotypes will be assigned to
#'   the validation set in multiple environments. Note that the value of
#'   prop_val at which this occurs will be lower for data that is unbalanced
#'   across environments. As prop_val is increased, the CV2 cross-validation
#'   scheme will begin to more closely resemble CV1, as a higher proportion of 
#'   genotypes will be assigned to the validation set across all environments.
#' @export
mc_train_val_multienv <- function(pheno, prop_val = 0.2, cv_scheme = "CV1") {
  
  cv_scheme   <- toupper(cv_scheme)
  pheno$TRAIN_VAL <- "train"
  pheno$CONCAT <- paste(pheno$ENV, pheno$IID, sep = "_")

  ## for CV1 cross-validation, define training/validation split by genotype (IID)
  if (cv_scheme == "CV1") {
      rand_iid <- sample(unique(pheno$IID), 
                         size = round(length(unique(pheno$IID)) * prop_val))
      pheno$TRAIN_VAL[pheno$IID %in% rand_iid] <- "val"
    
  ## CV2 is quite a bit more involved
  } else if (cv_scheme == "CV2") {
    
      ## First find unique geno-env combinations
      ph_uniq <- pheno[!duplicated(pheno$CONCAT), ]
      ph_uniq <- ph_uniq[c("IID", "ENV", "CONCAT")]
      
      ## Now calculate weighting of environments
      tot_nval <- round(prop_val * nrow(ph_uniq))
      env_weights <- as.data.frame(table(ph_uniq$ENV, dnn = c("ENV")), stringsAsFactors = FALSE)
      env_weights$prop <- env_weights$Freq / sum(env_weights$Freq)
      env_weights$nval <- round(env_weights$prop * tot_nval)
      env_weights <- env_weights[sample(nrow(env_weights)), ]
      
      ## Initialize two vectors - one a pool of all genotypes
      ## Second for holding output genotype-env combinations
      geno_pool <- unique(ph_uniq$IID)
      ge_vec <- vector()
      
      ## Cycle through envs
      for (e in env_weights$ENV) {
        
        ## number of geno-env combinations to draw for environment
        npull <- env_weights[env_weights == e, "nval"]
        
        ## draw and shuffle geno-env combinations matching target environment
        ## and contained within the pool of available genos
        one_env <- ph_uniq[ph_uniq$ENV == e & ph_uniq$IID %in% geno_pool, ]
        one_env <- one_env[sample(nrow(one_env)), ]
        
        ## If The pool of available genos becomes exhausted, use the remainder
        ## then refresh the genotype pool (minus the genotypes that were just used
        ## to avoid duplication)
        ## NOTE: Even if the genotype pool isn't technically depleted, refreshing
        ## can occur due to missing env-genotype combinations in the data
        ## (e.g. data that is unbalanced across environments)
        if (nrow(one_env) < npull) {
          top <- one_env
          geno_pool <- setdiff(unique(ph_uniq$IID), top$IID)
          bottom <- ph_uniq[ph_uniq$ENV == e & ph_uniq$IID %in% geno_pool, ]
          bottom <- bottom[sample(nrow(bottom)), ]
          one_env <- rbind(top, bottom)
        }
        
        ge_vec <- c(ge_vec, one_env$CONCAT[1:npull])
        geno_pool <- setdiff(geno_pool, one_env$IID[1:npull])
      } ## End of ENV loop
      
    pheno$TRAIN_VAL[pheno$CONCAT %in% ge_vec] <- "val"
  } ## End of CV-scheme if-else
  
  pheno$CONCAT <- NULL
  return(pheno)
}