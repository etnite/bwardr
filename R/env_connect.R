#' Evaluate connectivity between environments
#' 
#' @param pheno Dataframe containing one column for environment identifier, and
#' one column for genotype identifier.
#' @param env String containing the name of the column in pheno containing
#' environment identifiers.
#' @param geno String containing the name of the column in pheno containing
#' genotype identifiers.
#' @param proportion Logical indicating whether to return proportion of lines shared
#' between each pair of environments (i.e. the intersection of lines / union of lines) or
#' else just return the raw number of lines shared between environments
#' @param output String containing either "dataframe" to return a three-column
#' dataframe, or "matrix" to return a sparse matrix.
#' @return Either a dataframe or else a sparse matrix containing the number of
#' genotypes within each environment, and the number of genotypes shared between
#' each pair of environments.
#' @details This function takes as input a dataframe containing one column 
#' containing identifiers for environments, and one column containing identifiers
#' for genotypes, and returns a structure that lists the number of genotypes
#' evaluated within each environment, and the number of genotypes shared between
#' each pair of environments. This information can either be returned as a
#' dataframe with three columns, or else as a matrix.
#' @export
env_connect <- function(pheno, env, geno, proportion = FALSE, output = "dataframe") {
  
  ## Haven't yet figured out an elegant way to handle both dataframes and data.tables
  pheno <- as.data.frame(pheno)
  
  ## Sanity checks
  if (!output %in% c("dataframe", "matrix")) {
    stop("Please select either 'dataframe' or 'matrix' for output")
  }
  if (!env %in% colnames(pheno)) {
    stop(paste("env column", env, "is not present in pheno dataframe"))
  }
  if (!geno %in% colnames(pheno)) {
    stop(paste("geno column", geno, "is not present in pheno dataframe"))
  }
  if (!is.logical(proportion)) {
    stop("Please supply TRUE or FALSE for proportion parameter")
  }
  
  ## Find all pairs of envs
  env_vec <- unique(pheno[[env]])
  env_combos <- expand.grid(env_vec, env_vec)
  colnames(env_combos) <- c("env1", "env2")
  env_combos$n_genos <- 0

  ## Find number of genos in each env and between each pair of envs
  sub_pheno <- pheno[, c(env, geno)]
  for (i in 1:nrow(env_combos)) {
    e1_df <- sub_pheno[sub_pheno[[env]] == env_combos$env1[i], ]
    e2_df <- sub_pheno[sub_pheno[[env]] == env_combos$env2[i], ]
    shared_genos <- intersect(e1_df[[geno]], e2_df[[geno]])
    denom <- ifelse(proportion, length(union(e1_df[[geno]], e2_df[[geno]])), 1)
    if (length(shared_genos) > 0) {
      env_combos$n_genos[i] <- length(shared_genos) / denom
    }
  }
  
  ## Return dataframe or sparse matrix
  if (output == "dataframe") {
      return(env_combos)
  } else {
      return(with(env_combos, 
                  Matrix::sparseMatrix(i = as.numeric(env1),
                                       j = as.numeric(env2),
                                       x = n_genos,
                                       dimnames = list(levels(env1), levels(env2)))))
  }
}