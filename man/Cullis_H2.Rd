% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cullis_H2.R
\name{Cullis_H2}
\alias{Cullis_H2}
\title{Calculate Generalized Heritability from lme4 Model}
\usage{
Cullis_H2(model, geno_label = "GENO")
}
\arguments{
\item{model}{A lme4 model object}

\item{geno_label}{A string denoting the label of the random genotypic effect in the
supplied lme4 model object}
}
\value{
A list containing the following elements:
\itemize{
\item avsed The average standard error of differences between adjusted means estimates
\item H2 The generalized heritability estimate
}
}
\description{
Calculate Generalized Heritability from lme4 Model
}
\details{
This function calculates generalized heritability using the method of
Cullis et al., 2006 (\url{https://doi.org/10.1198/108571106X154443}).
Specifically, their formula is H2 = 1 - (vblup / (2 * var_g)). Where the
generalized heritability (H2) is a function of the reliability of the BLUPs
(vblup - the average standard error of differences between BLUPs squared),
and the genotypic variance (var_g). This method can be used in unbalanced
applications where the traditional entry-mean heritability calculation will
give biased estimates. The method of doing this using lme4 is detailed by
Ben Bolker at \url{https://stackoverflow.com/questions/38697477/mean-variance-of-a-difference-of-blues-or-blups-in-lme4}.
This method yields values that are slightly different (I have observed up to
0.75\%) from ASReml-R's results. Another solution I came across at
\url{https://shantel-martinez.github.io/resources.html} seems to produce
results that are more divergent from ASReml-R's.
}
