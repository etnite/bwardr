% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mc_train_val_multienv.R
\name{mc_train_val_multienv}
\alias{mc_train_val_multienv}
\title{Monte Carlo training/validation splits for multi-environment data}
\usage{
mc_train_val_multienv(pheno, prop_val = 0.2, cv_scheme = "CV1")
}
\arguments{
\item{pheno}{Dataframe containing phenotypic data. The dataframe must 
contain genotype IDs in a column named "IID", and environment designators 
in a column named "ENV". The data can be replicated within environment, and
can have any additional number of columns.}

\item{prop_val}{Number between 0 and 1 indicating the proportion of phenotypic
data that should be set to the validation set.}

\item{cv_scheme}{Character string consisting of either "CV1" or "CV2". "CV1"
assigns validation data by genotype (i.e. simulates introducing new genotypes
which have not been tested in any environment). CV2 assigns validation data
by genotype-environment combination (i.e. simulates introducing genotypes into
new environments).}
}
\value{
A dataframe identical to the input phenotypic dataframe, except with
  a VAL_SET column added (if not previously present) or else re-randomized 
  (if the column was already present) to indicate lines assigned to the
  validation set
}
\description{
Monte Carlo training/validation splits for multi-environment data
}
\details{
This function only performs Monte Carlo (i.e. random subsampling)
  cross-validation training/validation assignment. Note that k-fold cross
  validation becomes difficult to perform for a CV2 scheme. 
  If the CV2 cross-validation scheme is selected, the function will
  attempt to choose genotype-environment combinations to assign to the
  validation set such that any given genotype will only be assigned to the
  validation set in at most one environment. However, if the user sets a
  prop_val value higher than (1 / # envs), then some genotypes will be assigned to
  the validation set in multiple environments. Note that the value of
  prop_val at which this occurs will be lower for data that is unbalanced
  across environments. As prop_val is increased, the CV2 cross-validation
  scheme will begin to more closely resemble CV1, as a higher proportion of 
  genotypes will be assigned to the validation set across all environments.
}
