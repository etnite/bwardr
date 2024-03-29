% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_plots.R
\name{simulate_plots}
\alias{simulate_plots}
\title{Simulate a multi-environment breeding trial}
\usage{
simulate_plots(
  X,
  snps = 1000,
  inds = 100,
  qtls = 20,
  n_envs = 4,
  n_reps = 2,
  var_phen = 10,
  h2 = 0.75,
  var_env = 5,
  var_rep = 0.5,
  beta_ab = 5,
  return_scaled = FALSE
)
}
\arguments{
\item{X}{A Marker matrix in minor-allele dosage format, with individuals in
rows and SNPs in columns}

\item{snps}{Either an integer n, in which case n SNPs are randomly sampled out
of X, or a vector of SNP indices, or a vector of SNP names. Setting to "all"
will use all SNPs}

\item{inds}{Either an integer n, in which case n individuals are randomly sampled out
of X, or a vector of individual indices, or a vector of individual names. Setting
to "all" will use all individuals}

\item{qtls}{Either an integer n, in which case n evenly-spaced SNPs will be
chosen as QTLs, or a vector of SNP indices to assign as QTLs, or vector of
SNP names to assign as QTLs}

\item{n_envs}{Integer - number of environments to simulate}

\item{n_reps}{Integer - number of replicates within each environment to simulate}

\item{var_phen}{Positive float - Phenotypic variance within each environment}

\item{h2}{Proper fraction - Within-environment narrow-sense heritability value}

\item{var_env}{Positive float - Variance between environments}

\item{var_rep}{Positive float - Variance between replications within same environment}

\item{beta_ab}{Improper fraction - Beta distribution shape parameter to control
GxE effect. See details.}

\item{return_scaled}{Logical - Indicates whether to scale the genetic signal and
the error by the phenotypic standard deviation before returning the outputted
plots data}
}
\value{
A list containing the following elements:
\itemize{
\item plots_data - Dataframe containing data for individual plots in each environment
\item qtl_effects - Dataframe containing QTL effects in each environment
}
}
\description{
Simulate a multi-environment breeding trial
}
\details{
QTL effects are allowed to vary across environments by sampling out
of a symmetric beta distribution. This implies that genotype-by-environment
(GxE) interaction decreases as the α and β shape parameters of the
distribution increase. At the limits, setting α = β = 1 makes the
beta distribution equivalent to a uniform distribution - QTL effects may vary
without any central tendency. Alternatively, setting α = β = \code{Inf}
will set QTL effects constant across environments.
Setting n_reps to 1 will simulate within-environment means. In this
case, the var_rep value will have no effect.
}
