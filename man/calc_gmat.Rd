% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_gmat.R
\name{calc_gmat}
\alias{calc_gmat}
\title{Calculate Realized Genomic Relationship Matrices}
\usage{
calc_gmat(data, method = "GOF")
}
\arguments{
\item{data}{Dataframe or matrix with individuals in rows and markers in columns.
Markers should be encoded in minor-allele dosage format
(i.e. 0 = homozygous for major allele, 1 = heterozygous, 2 = homozygous
for minor allele).}

\item{method}{String to select the method to use for calculating the 
relationship matrix. Available methods are:
\enumerate{
  \item "GOF" - observed allele frequencies (VanRaden, 2008)
  \item "GD" - weighted markers by recipricals of expected variance (Forni et al., 2011)
  \item "G05" - allele frequencies fixed at 0.5 (Forni et al., 2011)
  \item "GM" - allele frequencies fixed at mean for each locus (Forni et al., 2011)
}}
}
\value{
The full realized genetic relationship matrix (G)
}
\description{
Calculate Realized Genomic Relationship Matrices
}
\details{
This function is adapted from an example from Chapter 11 of
  Genetic Data Analysis for Plant and Animal Breeding, by Isik, Holland and
  Maltecca, 2017 (DOI 10.1007/978-3-319-55177-7). However, their code is
  more complex, as it includes one additional computation method that uses
  observed pedigree data, and it outputs both realized genetic and additive
  relationships. This function only returns the realized genetic relationship
  matrix.
}
