#' Simulate a multi-environment breeding trial
#'
#' @param X A Marker matrix in minor-allele dosage format, with individuals in
#' rows and SNPs in columns
#' @param snps Either an integer n, in which case n SNPs are randomly sampled out
#' of X, or a vector of SNP indices, or a vector of SNP names. Setting to "all"
#' will use all SNPs
#' @param inds Either an integer n, in which case n individuals are randomly sampled out
#' of X, or a vector of individual indices, or a vector of individual names. Setting
#' to "all" will use all individuals
#' @param qtls Either an integer n, in which case n evenly-spaced SNPs will be
#' chosen as QTLs, or a vector of SNP indices to assign as QTLs, or vector of 
#' SNP names to assign as QTLs
#' @param n_envs Integer - number of environments to simulate
#' @param n_reps Integer - number of replicates within each environment to simulate
#' @param var_phen Positive float - Phenotypic variance within each environment
#' @param h2 Proper fraction - Within-environment narrow-sense heritability value
#' @param var_env Positive float - Variance between environments
#' @param var_rep Positive float - Variance between replications within same environment
#' @param beta_ab Improper fraction - Beta distribution shape parameter to control
#' GxE effect. See details.
#' @param return_scaled Logical - Indicates whether to scale the genetic signal and
#'   the error by the phenotypic standard deviation before returning the outputted
#'   plots data
#' @return A list containing the following elements:
#' * plots_data - Dataframe containing data for individual plots in each environment
#' * qtl_effects - Dataframe containing QTL effects in each environment
#' @details QTL effects are allowed to vary across environments by sampling out
#' of a symmetric beta distribution. This implies that genotype-by-environment
#' (GxE) interaction decreases as the α and β shape parameters of the
#' distribution increase. At the limits, setting α = β = 1 makes the
#' beta distribution equivalent to a uniform distribution - QTL effects may vary
#' without any central tendency. Alternatively, setting α = β = `Inf`
#' will set QTL effects constant across environments.
#' Setting n_reps to 1 will simulate within-environment means. In this
#' case, the var_rep value will have no effect.
#' @md 
#' @export
simulate_plots <- function(X, snps = 1000, inds = 100, qtls = 20, n_envs = 4,
                           n_reps = 2, var_phen = 10, h2 = 0.75, var_env = 5,
                           var_rep = 0.5, beta_ab = 5, return_scaled = FALSE) {
  
  ## Sanity checks on input
  stopifnot(n_envs >= 1)
  stopifnot(n_reps >= 1)
  stopifnot(var_phen > 0)
  stopifnot(h2 > 0 & h2 <= 1)
  stopifnot(var_env > 0)
  stopifnot(var_rep >= 0)
  stopifnot(beta_ab >= 1)
  
  ## Subset and scale input SNP matrix
  if (snps != "all") {
    if (length(snps) == 1) {
      X <- X[, sample(1:ncol(X), snps)]
    } else if (length(snps) > 1) {
      X <- X[, snps]
    }
  }
  
  if (inds == "all") {
    n_inds <- nrow(X)
  } else {
    if (length(inds) == 1) {
      n_inds <- inds
      X <- X[sample(1:nrow(X), inds), ]
    } else if (length(inds) > 1) {
      n_inds <- length(inds)
      X <- X[inds, ]
    }
  }
  
  X <- scale(X)
  
  ## Assign QTLs
  if (length(qtls) == 1) {
    n_qtls <- qtls
    QTL <- floor(seq(from = 50, to = ncol(X), length = n_qtls))
    qtl_names <- colnames(X)[QTL]
  } else {
    n_qtls <- length(qtls)
    if (is.numeric(qtls)) {
      QTL <- qtls
      qtl_names <- colnames(X)[QTL]
    } else if (is.character(qtls)) {
      QTL <- which(colnames(X) %in% qtls)
      qtl_names <- qtls
    }
  }
  
  ## Generate environment means
  env_means <- stats::rnorm(n_envs, sd = sqrt(var_env))
  names(env_means) <- paste0("env", 1:n_envs)
  
  ## Initialize lists for output data
  qtl_list <- vector("list", length = n_envs)
  re_combos <- expand.grid(1:n_envs, 1:n_reps)
  sim_list <- vector("list", length = n_envs * n_reps)
  names(sim_list) <- paste(re_combos$Var1, re_combos$Var2, sep = "_")
  
  ## Loop through envs and reps
  for (e in 1:n_envs) {
    
    ## Vary QTL effects across envs. unless the beta shape parameter is set to infinite
    if (is.infinite(beta_ab)) {
      b <- rep(1, n_qtls) * sqrt(h2 / n_qtls)
    } else {
      ge_dev <- stats::rbeta(n_qtls, beta_ab, beta_ab)
      b <- ge_dev * n_qtls / sum(ge_dev) * sqrt(h2 / n_qtls)
    }
    qtl_list[[e]] <- data.frame("ENV" = paste0("env", e),
                                "QTL" = qtl_names,
                                "EFFECT" = b)
    
    ## Create genetic component of phenotype
    signal <- X[, QTL] %*% b
    
    ## Generate phenotype
    for (r in 1:n_reps) {
      index <- paste(e, r, sep = "_")
      rep_mean <- stats::rnorm(1, sd = sqrt(var_rep))
      error <- stats::rnorm(n_inds, sd = sqrt(1 - h2))
      sim_list[[index]] <- data.frame("IID" = rownames(X),
                                      "ENV" = paste0("env", e),
                                      "REP" = paste0("rep", r),
                                      "GEN" = signal,
                                      "ERROR" = error,
                                      "PHEN" = sqrt(var_phen) * (signal + error) + env_means[e] + rep_mean)
    }
  }
  plots_data <- do.call("rbind", sim_list)
  qtl_effects <- do.call("rbind", qtl_list)
  
  ## Optionally scale genetic signal and error prior to outputting
  if (return_scaled) {
    plots_data$GEN <- plots_data$GEN * sqrt(var_phen)
    plots_data$ERROR <- plots_data$ERROR * sqrt(var_phen)
  }
  
  return(list("plots_data" = plots_data, "qtl_effects" = qtl_effects))  
}